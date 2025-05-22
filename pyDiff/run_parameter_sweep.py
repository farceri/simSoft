import yaml
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import time  # For adding delays if needed, and for unique temp files

# --- Configuration for the Parameter Sweep ---
BASE_CONFIG_FILE = "config_template.yaml"
PYTHON_EXECUTABLE = "python"  # Or "python3"
SIMULATION_SCRIPT = "RandomWalk-optimized.py"  # Your main simulation script
SWEEP_RESULTS_CSV = "truncated_gaussian_sweep_results.csv"  # Updated name for clarity

# Define the range of disorder parameters for 'truncated_gaussian_landscape_per_trial_vectorized'
# These are the parameters of the *underlying* Gaussian before truncation to [0,1]
mean_rest_prob_values = [-0.9,-0,8,-0.7]  # Example: Mean of the underlying Gaussian
std_dev_rest_prob_values = [0.01,0.05,0,75,0.01]  # Example: Std dev of the underlying Gaussian

# Parameters for t_diff and D extraction (tune these as needed)
SLOPE_TOLERANCE = 1e-4
INITIAL_TRANSIENT_TIME = 1.0  # seconds


# --- Function to Extract t_diff and D (Same as before) ---
def extract_diffusion_parameters(data_filename, slope_tolerance, initial_transient_threshold):
    """
    Loads data from the simulation output file, calculates the slope of MSD/Time,
    estimates t_diff, and then calculates the diffusion coefficient D.
    Returns a dictionary {'t_diff': value, 'D': value, 'avg_msd_over_time': value}.
    Returns np.nan for values if they cannot be determined.
    """
    try:
        data = np.loadtxt(data_filename, delimiter=',', skiprows=1, usecols=(0, 2))
        time_values = data[:, 0]
        msd_over_time_values = data[:, 1]

        if len(time_values) < 5:
            print(f"  Not enough data points in {data_filename} to analyze.")
            return {'t_diff': np.nan, 'D': np.nan, 'avg_msd_over_time': np.nan}

        slope_of_msd_over_time = np.gradient(msd_over_time_values, time_values)

        t_diff_val = np.nan
        d_val = np.nan
        avg_msd_over_time_val = np.nan

        indices_after_transient = np.where(time_values >= initial_transient_threshold)[0]

        if indices_after_transient.size > 0:
            first_search_idx = indices_after_transient[0]
            potential_tdiff_indices_relative = np.where(
                np.abs(slope_of_msd_over_time[first_search_idx:]) < slope_tolerance
            )[0]

            if potential_tdiff_indices_relative.size > 0:
                t_diff_idx = first_search_idx + potential_tdiff_indices_relative[0]
                t_diff_val = time_values[t_diff_idx]

                if t_diff_idx < len(msd_over_time_values) - 1:
                    msd_over_time_linear_regime = msd_over_time_values[t_diff_idx:]
                    if msd_over_time_linear_regime.size > 1:
                        avg_msd_over_time_val = np.mean(msd_over_time_linear_regime)
                        d_val = avg_msd_over_time_val / 4.0  # Assuming 2D

        return {'t_diff': t_diff_val, 'D': d_val, 'avg_msd_over_time': avg_msd_over_time_val}

    except FileNotFoundError:
        print(f"  Error: Data file {data_filename} not found for analysis.")
        return {'t_diff': np.nan, 'D': np.nan, 'avg_msd_over_time': np.nan}
    except Exception as e:
        print(f"  Error analyzing {data_filename}: {e}")
        return {'t_diff': np.nan, 'D': np.nan, 'avg_msd_over_time': np.nan}


# --- Main Loop for Parameter Sweep ---
collected_results = []

try:
    with open(BASE_CONFIG_FILE, 'r') as f:
        config_template = yaml.safe_load(f)
except FileNotFoundError:
    print(f"Error: Base config file '{BASE_CONFIG_FILE}' not found. Exiting.")
    exit()

output_dir_name = "run_data_outputs_truncated_gaussian"  # Or adjust based on disorder type
os.makedirs(output_dir_name, exist_ok=True)

for mean_rp in mean_rest_prob_values:
    for std_dev_rp in std_dev_rest_prob_values:
        variance_rp = std_dev_rp ** 2

        print(
            f"\nProcessing: mean_rest_prob = {mean_rp:.2f}, std_dev_rest_prob = {std_dev_rp:.2f} (variance = {variance_rp:.4f})")

        # Define unique temporary config filename for this iteration
        # Using a more robust way to ensure uniqueness and placement if script is moved
        temp_config_filename = os.path.join(os.getcwd(), f"temp_config_{int(time.time() * 100000 + os.getpid())}.yaml")

        run_data_filename = ""  # Initialize

        try:
            # 1. Modify Configuration
            current_config = config_template.copy()

            current_config['disorder']['type'] = 'truncated_gaussian'  # Or your chosen type
            current_config['disorder']['params']['mean_rest_prob'] = float(mean_rp)
            current_config['disorder']['params']['std_dev_rest_prob'] = float(std_dev_rp)
            current_config['disorder']['params'].pop('max_allowed_rest_prob', None)
            current_config['disorder']['params'].pop('min_rest_prob', None)
            current_config['disorder']['params'].pop('rest_level', None)

            mean_str = str(mean_rp).replace('.', 'p')
            std_str = str(std_dev_rp).replace('.', 'p')
            run_data_filename = os.path.join(output_dir_name,
                                             f"results_type_{current_config['disorder']['type']}_mean{mean_str}_std{std_str}.txt")
            current_config['save_data']['filename'] = run_data_filename
            current_config['save_data']['enabled'] = True

            if 'plotting' not in current_config:  # Ensure plotting section exists
                current_config['plotting'] = {}
            current_config['plotting']['show_plots_at_end'] = False  # Suppress plots

            if 'animation' not in current_config: current_config['animation'] = {}
            current_config['animation']['enabled'] = False
            if 'histograms' not in current_config: current_config['histograms'] = {}
            current_config['histograms']['enabled'] = False

            with open(temp_config_filename, 'w') as f:
                yaml.dump(current_config, f)

            # 2. Run Simulation
            print(f"  Running simulation with {temp_config_filename}...")
            extracted_data = {'t_diff': np.nan, 'D': np.nan, 'avg_msd_over_time': np.nan}  # Default for this iteration

            try:
                process_result = subprocess.run(
                    [PYTHON_EXECUTABLE, SIMULATION_SCRIPT, temp_config_filename],
                    capture_output=True, text=True, check=True, timeout=3600
                )
                print(f"  --- Subprocess STDOUT for mean={mean_rp:.2f}, std={std_dev_rp:.2f} ---")
                print(process_result.stdout)
                print(f"  --- Subprocess STDERR for mean={mean_rp:.2f}, std={std_dev_rp:.2f} ---")
                if process_result.stderr:
                    print(f"  WARNING: STDERR is not empty!")
                print(f"  --- End Subprocess Output ---")
                print(
                    f"  Simulation for mean={mean_rp:.2f}, std={std_dev_rp:.2f} completed (Return Code: {process_result.returncode}).")

                # Now check for the file
                if os.path.exists(run_data_filename):

                # 3. Extract t_diff and D
                    print(f"  Analyzing output file: {run_data_filename}")
                    extracted_data = extract_diffusion_parameters(run_data_filename, SLOPE_TOLERANCE,
                                                                  INITIAL_TRANSIENT_TIME)
                    print(
                        f"  Extracted t_diff: {extracted_data.get('t_diff', float('nan')):.4f}, D: {extracted_data.get('D', float('nan')):.4e}, AvgMSD/t: {extracted_data.get('avg_msd_over_time', float('nan')):.4e}")
                else:
                    print(f"  Output file {run_data_filename} not found after simulation. Cannot extract parameters.")

            except subprocess.CalledProcessError as e:
                print(f"  ERROR during simulation run for mean={mean_rp:.2f}, std={std_dev_rp:.2f}:")
                print(f"  STDOUT: {e.stdout}")
                print(f"  STDERR: {e.stderr}")
            except subprocess.TimeoutExpired:
                print(f"  TIMEOUT during simulation run for mean={mean_rp:.2f}, std={std_dev_rp:.2f}.")

            # 4. Store Results
            collected_results.append({
                'mean_input_rest_prob': mean_rp,
                'std_dev_input_rest_prob': std_dev_rp,
                'variance_input_rest_prob': variance_rp,
                't_diff': extracted_data['t_diff'],
                'D': extracted_data['D'],
                'avg_msd_over_time_at_D': extracted_data['avg_msd_over_time']
            })

        finally:
            # *** Ensure temporary config file is always deleted ***
            if os.path.exists(temp_config_filename):
                try:
                    os.remove(temp_config_filename)
                    print(f"  Cleaned up temporary config: {temp_config_filename}")
                except Exception as e_rem:
                    print(f"  Warning: Could not remove temporary config {temp_config_filename}: {e_rem}")
            # *****************************************************
# 5. Save Aggregated Results
results_df = pd.DataFrame(collected_results)
results_df.to_csv(SWEEP_RESULTS_CSV, index=False, na_rep='NaN')
print(f"\nParameter sweep complete. All results saved to {SWEEP_RESULTS_CSV}")

# 6. Plot Final Results
if not results_df.empty:
    # Ensure columns exist before plotting and handle NaNs for plotting
    # For D vs mean (color by variance)
    plot_df1 = results_df.dropna(subset=['mean_input_rest_prob', 'D', 'variance_input_rest_prob'])
    # For t_diff vs variance (color by mean)
    plot_df2 = results_df.dropna(subset=['variance_input_rest_prob', 't_diff', 'mean_input_rest_prob'])

    plt.figure(figsize=(14, 6)) # Keep figure size or adjust

    # --- Plot 1: D vs. mean_input_rest_prob (color by variance_input_rest_prob) ---
    plt.subplot(1, 2, 1)
    if not plot_df1.empty:
        scatter1 = plt.scatter(plot_df1['mean_input_rest_prob'], plot_df1['D'],
                               c=plot_df1['variance_input_rest_prob'], cmap='viridis', # Changed X and C
                               alpha=0.7, edgecolors='k', linewidths=0.5)
        plt.xlabel('Mean of Underlying Rest Prob Distribution') # Changed X label
        plt.ylabel('Estimated Diffusion Coefficient (D)')
        plt.title('D vs. Mean of Rest Prob (Color by Variance)') # Changed Title
        cbar1 = plt.colorbar(scatter1, label='Variance of Underlying Rest Prob') # Changed colorbar label
        plt.grid(True, linestyle=':')
    else:
        plt.text(0.5, 0.5, 'No valid data for D vs Mean plot', horizontalalignment='center', verticalalignment='center')


    # --- Plot 2: t_diff vs. variance_input_rest_prob (color by mean_input_rest_prob) ---
    plt.subplot(1, 2, 2)
    if not plot_df2.empty:
        scatter2 = plt.scatter(plot_df2['variance_input_rest_prob'], plot_df2['t_diff'],
                               c=plot_df2['mean_input_rest_prob'], cmap='plasma', # Changed X and C
                               alpha=0.7, edgecolors='k', linewidths=0.5)
        plt.xlabel('Variance of Underlying Rest Prob Distribution (std_dev^2)') # Changed X label
        plt.ylabel('Estimated Convergence Time t_diff (s)')
        plt.title('t_diff vs. Variance of Rest Prob (Color by Mean)') # Changed Title
        cbar2 = plt.colorbar(scatter2, label='Mean of Underlying Rest Prob') # Changed colorbar label
        plt.grid(True, linestyle=':')
    else:
        plt.text(0.5, 0.5, 'No valid data for t_diff vs Variance plot', horizontalalignment='center', verticalalignment='center')

    plt.tight_layout()
    plt.savefig("sweep_summary_plots_MEDIABASSA_DEVST09.png") # Changed save filename slightly
    plt.show()
else:
    print("No results collected to plot.")