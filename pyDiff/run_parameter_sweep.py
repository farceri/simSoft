import yaml
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import time  # For adding delays if needed, and for unique temp files
import re  # For parsing output

# --- Configuration for the Parameter Sweep ---
BASE_CONFIG_FILE = "config_template.yaml"
PYTHON_EXECUTABLE = "python"  # Or "python3"
SIMULATION_SCRIPT = "RandomWalk-optimized.py"  # Your main simulation script
ANALYSIS_MODE = "var_D"
# Updated name to reflect the inclusion of actual mean/variance
# Update CSV name based on analysis mode or make it more general
if ANALYSIS_MODE == "var_D":
    SWEEP_RESULTS_CSV = "var_D_sweep_results_PROVA.csv"
else:
    SWEEP_RESULTS_CSV = "t_diff_D_sweep_results.csv"

# Define the range of disorder parameters for 'truncated_gaussian_landscape_per_trial_vectorized'
mean_rest_prob_values = [0]
#std_dev_rest_prob_values = [0.01,0.05,0.07,0.1]
std_dev_rest_prob_values = [0.3,0.5,1]
#std_dev_rest_prob_values = [0.3,0.5,0.6,0.9,1]

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

output_dir_name = f"run_data_outputs_{ANALYSIS_MODE}" # Directory name based on mode
os.makedirs(output_dir_name, exist_ok=True)

for mean_rp in mean_rest_prob_values:
    for std_dev_rp in std_dev_rest_prob_values:
        variance_rp = std_dev_rp ** 2

        print(
            f"\nProcessing: mean_rest_prob = {mean_rp:.2f}, std_dev_rest_prob = {std_dev_rp:.2f} (variance = {variance_rp:.4f})")

        temp_config_filename = os.path.join(os.getcwd(), f"temp_config_{int(time.time() * 100000 + os.getpid())}.yaml")
        run_data_filename = ""

        # Initialize all potential values for this iteration
        actual_mean_prest_val = np.nan  # E[ E[rho|landscape_i] ] (average of trial means)
        avg_actual_var_prest_val = np.nan  # E[ Var(rho|landscape_i) ] (average of trial variances)
        var_of_means_prest_val = np.nan  # Var( E[rho|landscape_i] ) (variance of trial means)

        extracted_t_diff = np.nan
        extracted_D_avg = np.nan
        extracted_avg_msd_over_time = np.nan
        parsed_var_D_late_time = np.nan

        try:
            # 1. Modify Configuration
            current_config = config_template.copy()

            # This should remain to ensure RandomWalk-optimized.py calculates actual mean/var
            current_config['disorder']['type'] = 'truncated_gaussian'
            current_config['disorder']['params']['mean_rest_prob'] = float(mean_rp)
            current_config['disorder']['params']['std_dev_rest_prob'] = float(std_dev_rp)
            # ... (remove other disorder params as in your original script) ...
            current_config['disorder']['params'].pop('max_allowed_rest_prob', None)
            current_config['disorder']['params'].pop('min_rest_prob', None)
            current_config['disorder']['params'].pop('rest_level', None)


            mean_str = str(mean_rp).replace('.', 'p')
            std_str = str(std_dev_rp).replace('.', 'p')
            # Ensure run_data_filename is defined before use, e.g.
            run_data_filename = os.path.join(output_dir_name,
                                             f"results_type_{current_config['disorder']['type']}_mode{ANALYSIS_MODE}_mean{mean_str}_std{std_str}.txt")
            current_config['save_data']['filename'] = run_data_filename
            current_config['save_data']['enabled'] = True # Raw data always saved by RW-opt.py

            # Suppress plots/animation in simulation script
            if 'plotting' not in current_config: current_config['plotting'] = {}
            current_config['plotting']['show_plots_at_end'] = False
            if 'animation' not in current_config: current_config['animation'] = {}
            current_config['animation']['enabled'] = False
            if 'histograms' not in current_config: current_config['histograms'] = {}
            current_config['histograms']['enabled'] = False

            with open(temp_config_filename, 'w') as f:
                yaml.dump(current_config, f)

            # Initialize all possible result fields that might be extracted
            extracted_t_diff = np.nan
            extracted_D_avg = np.nan
            extracted_avg_msd_over_time = np.nan
            parsed_var_D_late_time = np.nan
            # actual_mean_prest_val and actual_var_prest_val are initialized outside this try,
            # but re-parsing them from stdout for each run is good practice.

            # Build the command for subprocess
            cmd = [PYTHON_EXECUTABLE, SIMULATION_SCRIPT, temp_config_filename]
            if ANALYSIS_MODE == "var_D":
                cmd.append("--analysis_mode")
                cmd.append("var_D")
            # else it defaults to 't_diff_D' in RandomWalk-optimized.py

            print(f"  Running simulation with command: {' '.join(cmd)}")
            # The try block for subprocess.run and parsing:
            try:
                process_result = subprocess.run(
                    cmd, # Use the constructed command
                    capture_output=True, text=True, check=True, timeout=7200 # Increased timeout slightly
                )
                print(f"  --- Subprocess STDOUT for mean={mean_rp:.2f}, std={std_dev_rp:.2f} ---")
                # print(process_result.stdout) # Can be very verbose, print selectively or upon error
                print(f"  --- Subprocess STDERR for mean={mean_rp:.2f}, std={std_dev_rp:.2f} ---")
                if process_result.stderr:
                    print(process_result.stderr) # Print if not empty
                    # print(f"  WARNING: STDERR is not empty!")
                print(f"  --- End Subprocess Output ---")
                print(f"  Simulation completed (Return Code: {process_result.returncode}).")

                # Parse stdout
                if process_result.stdout:
                    for line in process_result.stdout.splitlines():
                        if "SWEEP_DATA_ACTUAL_MEAN_PREST:" in line:  # This is avg_actual_mean_prest
                            try:
                                actual_mean_prest_val = float(line.split(":")[1])
                            except:
                                print(f"  Warning: Could not parse ACTUAL_MEAN_PREST: {line}")
                        # MODIFIED PARSING KEY:
                        if "SWEEP_DATA_AVG_ACTUAL_VAR_PREST:" in line:  # This is avg_actual_var_prest
                            try:
                                avg_actual_var_prest_val = float(line.split(":")[1])
                            except:
                                print(f"  Warning: Could not parse AVG_ACTUAL_VAR_PREST: {line}")
                        # NEW PARSING KEY:
                        if "SWEEP_DATA_VAR_OF_MEANS_PREST:" in line:  # This is var_of_trial_means_prest
                            try:
                                var_of_means_prest_val = float(line.split(":")[1])
                            except:
                                print(f"  Warning: Could not parse VAR_OF_MEANS_PREST: {line}")

                        if ANALYSIS_MODE == "var_D":
                            if "SWEEP_DATA_VAR_D_LATE_TIME:" in line:
                                try:
                                    parsed_var_D_late_time = float(line.split(":")[1])
                                except:
                                    print(f"  Warning: Could not parse VAR_D_LATE_TIME: {line}")

                print(f"  Parsed actual_mean_prest: {actual_mean_prest_val}")
                print(f"  Parsed avg_actual_var_prest: {avg_actual_var_prest_val}")
                print(f"  Parsed var_of_means_prest: {var_of_means_prest_val}")


                if ANALYSIS_MODE == "var_D":
                    print(f"  Analysis Mode: var_D. Parsing Var(D) from stdout.")
                    if process_result.stdout:
                        for line in process_result.stdout.splitlines():
                            if "SWEEP_DATA_VAR_D_LATE_TIME:" in line:
                                try:
                                    parsed_var_D_late_time = float(line.split(":")[1])
                                except (ValueError, IndexError):
                                    print(f"  Warning: Could not parse Var(D) from line: {line}")
                    print(f"  Parsed Var(D)_late_time: {parsed_var_D_late_time}")
                    if not os.path.exists(run_data_filename): # Check if the raw data file was created
                         print(f"  Warning: Simulation output file {run_data_filename} not found (though not used for t_diff in this mode).")

                else: # Default "t_diff_D" mode
                    print(f"  Analysis Mode: t_diff_D. Running t_diff/D extraction.")
                    if os.path.exists(run_data_filename):
                        extracted_params = extract_diffusion_parameters(run_data_filename, SLOPE_TOLERANCE, INITIAL_TRANSIENT_TIME)
                        extracted_t_diff = extracted_params['t_diff']
                        extracted_D_avg = extracted_params['D']
                        extracted_avg_msd_over_time = extracted_params['avg_msd_over_time']
                        print(
                            f"  Extracted t_diff: {extracted_t_diff:.4f}, D_avg: {extracted_D_avg:.4e}, AvgMSD/t: {extracted_avg_msd_over_time:.4e}"
                        )
                    else:
                        print(f"  Output file {run_data_filename} not found. Cannot extract t_diff/D.")

            except subprocess.CalledProcessError as e:
                print(f"  ERROR during simulation run for mean={mean_rp:.2f}, std={std_dev_rp:.2f}:")
                print(f"  STDOUT: {e.stdout}")
                print(f"  STDERR: {e.stderr}")
            except subprocess.TimeoutExpired:
                print(f"  TIMEOUT during simulation run for mean={mean_rp:.2f}, std={std_dev_rp:.2f}.")
            except Exception as e_run: # Catch any other unexpected errors during run/parse
                 print(f"  UNEXPECTED ERROR during simulation run or parsing for mean={mean_rp:.2f}, std={std_dev_rp:.2f}: {e_run}")


            # Store results
            collected_results.append({
                'mean_input_rest_prob': mean_rp,
                'std_dev_input_rest_prob': std_dev_rp,
                'variance_input_rest_prob': variance_rp,
                'actual_mean_prest': actual_mean_prest_val,          # E[ E[rho|landscape_i] ]
                'avg_within_landscape_var_prest': avg_actual_var_prest_val, # E[ Var(rho|landscape_i) ]
                'var_between_landscape_mean_prest': var_of_means_prest_val, # Var( E[rho|landscape_i] )
                't_diff': extracted_t_diff,
                'D_avg': extracted_D_avg,
                'avg_msd_over_time_at_D': extracted_avg_msd_over_time,
                'var_D_late_time': parsed_var_D_late_time
            })
        finally:
            if os.path.exists(temp_config_filename):
                try:
                    os.remove(temp_config_filename)
                    print(f"  Cleaned up temporary config: {temp_config_filename}")
                except Exception as e_rem:
                    print(f"  Warning: Could not remove temporary config {temp_config_filename}: {e_rem}")

results_df = pd.DataFrame(collected_results)
results_df.to_csv(SWEEP_RESULTS_CSV, index=False, na_rep='NaN')
print(f"\nParameter sweep complete ({ANALYSIS_MODE} mode). All results saved to {SWEEP_RESULTS_CSV}")

# Plot Final Results (Conditional Plotting based on ANALYSIS_MODE)
if not results_df.empty:
    plt.figure(figsize=(14, 7)) # Single figure for conditional plots

    if ANALYSIS_MODE == "var_D":
        if 'avg_within_landscape_var_prest' in results_df.columns and \
                'var_between_landscape_mean_prest' in results_df.columns:
            results_df['total_actual_var_prest'] = results_df['avg_within_landscape_var_prest'].fillna(0) + \
                                                   results_df['var_between_landscape_mean_prest'].fillna(0)
        else:
            print(
                "Warning: Source columns for 'total_actual_var_prest' are missing in results_df. Plotting may fail or use fallbacks.")
            results_df['total_actual_var_prest'] = np.nan  # Ensure column exists to prevent other KeyErrors
        # Plot Var(D)_late_time vs. actual_variance_rest_prob
        plot_df_varD = results_df.dropna(subset=['var_D_late_time', 'total_actual_var_prest', 'actual_mean_prest'])

        if not plot_df_varD.empty:
            plt.scatter(plot_df_varD['total_actual_var_prest'], plot_df_varD['var_D_late_time'],
                        c=plot_df_varD['actual_mean_prest'], cmap='coolwarm',  # Use 'actual_mean_prest'
                        alpha=0.7, edgecolors='k', linewidths=0.5)
            plt.xlabel('Total Actual Variance of Rest Prob ($Var_{total}(\rho)$)')
            plt.ylabel('Variance of D (Late Time Regime)')
            plt.title('$Var(D)_{late}$ vs. $Var_{total}(\rho)$ (Color by Actual Mean)')
            plt.colorbar(label='Actual Mean of Rest Prob ($E[\\bar{\\rho}_i]$)')
            plt.grid(True, linestyle=':')
        else:
            plt.text(0.5, 0.5, 'No valid data for Var(D) plot after processing.', ha='center', va='center',
                     transform=plt.gca().transAxes)

    else: # "t_diff_D" mode plots
        # Plot 1: D_avg vs. actual_mean_rest_prob
        plt.subplot(1, 2, 1)
        plot_df1 = results_df.dropna(subset=['actual_mean_rest_prob', 'D_avg', 'actual_variance_rest_prob'])
        if not plot_df1.empty:
            scatter1 = plt.scatter(plot_df1['actual_mean_rest_prob'], plot_df1['D_avg'],
                                   c=plot_df1['actual_variance_rest_prob'], cmap='viridis',
                                   alpha=0.7, edgecolors='k', linewidths=0.5)
            plt.xlabel('Actual Mean of Rest Prob Distribution')
            plt.ylabel('Estimated Average Diffusion Coefficient (D_avg)')
            plt.title('D_avg vs. Actual Mean (Color by Actual Variance)')
            plt.colorbar(scatter1, label='Actual Variance of Rest Prob')
            plt.grid(True, linestyle=':')
        else:
            plt.text(0.5, 0.5, 'No valid data for D_avg vs Actual Mean plot', ha='center', va='center', transform=plt.gca().transAxes)

        # Plot 2: t_diff vs. actual_variance_rest_prob
        plt.subplot(1, 2, 2)
        plot_df2 = results_df.dropna(subset=['actual_variance_rest_prob', 't_diff', 'actual_mean_rest_prob'])
        if not plot_df2.empty:
            scatter2 = plt.scatter(plot_df2['actual_variance_rest_prob'], plot_df2['t_diff'],
                                   c=plot_df2['actual_mean_rest_prob'], cmap='plasma',
                                   alpha=0.7, edgecolors='k', linewidths=0.5)
            plt.xlabel('Actual Variance of Rest Prob Distribution')
            plt.ylabel('Estimated Convergence Time t_diff (s)')
            plt.title('t_diff vs. Actual Variance (Color by Actual Mean)')
            plt.colorbar(scatter2, label='Actual Mean of Rest Prob')
            plt.grid(True, linestyle=':')
        else:
            plt.text(0.5, 0.5, 'No valid data for t_diff vs Actual Variance plot', ha='center', va='center', transform=plt.gca().transAxes)

    plt.tight_layout()
    plt.savefig(f"sweep_summary_plots_{ANALYSIS_MODE}.png")
    plt.show()
else:
    print("No results collected to plot.")