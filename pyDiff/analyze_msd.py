import numpy as np
import matplotlib.pyplot as plt

# Define the filename
filename = "simulation_results_TEST.txt"  # Make sure this file is in the correct path

try:
    # Load data: Time (col 0), MSDoverTime (col 2)
    data = np.loadtxt(filename, delimiter=',', skiprows=1, usecols=(0, 2))
    time_values = data[:, 0]
    msd_over_time_values = data[:, 1]

    if len(time_values) < 2:
        raise ValueError("Not enough data points to calculate gradient.")

    # 1. Calculate the numerical derivative (slope) of MSDoverTime w.r.t. Time
    slope_of_msd_over_time = np.gradient(msd_over_time_values, time_values)

    # 2. Estimate t_diff
    # Parameters for t_diff estimation (these might need tuning based on your data)
    slope_tolerance = 1e-5  # How close to zero the slope should be to be considered "flat"
    # Inspect your slope plot to set a reasonable value.
    # This was chosen based on typical scales but might be too strict/loose.
    initial_transient_time_threshold = 1.2  # Ignore data before this time (s) for finding t_diff
    # This helps avoid premature t_diff due to initial noise/steepness.
    # From TEST.png, 2.0s looked like a good visual estimate.

    t_diff_determined_value = None
    t_diff_determined_index = -1
    average_msd_over_time_for_D = np.nan
    diffusion_coefficient_D = np.nan

    # Find indices of data points that are past the initial transient period
    indices_after_transient = np.where(time_values >= initial_transient_time_threshold)[0]

    if indices_after_transient.size > 0:
        first_index_to_search = indices_after_transient[0]
        # Look for the first point *after the transient period* where the slope is within tolerance
        # This is a simple heuristic: find the first time it becomes "flat".
        # A more robust method might check if it *stays* flat for a certain duration.
        potential_tdiff_indices_relative = np.where(
            np.abs(slope_of_msd_over_time[first_index_to_search:]) < slope_tolerance
        )[0]

        if potential_tdiff_indices_relative.size > 0:
            t_diff_determined_index = first_index_to_search + potential_tdiff_indices_relative[0]
            t_diff_determined_value = time_values[t_diff_determined_index]
            print(f"Estimated t_diff (convergence to linear regime) ≈ {t_diff_determined_value:.4f} s")

            # 3. Calculate D for t >= t_diff_determined_value
            # Consider data from this t_diff index onwards
            if t_diff_determined_index < len(msd_over_time_values) - 1:  # Ensure there's data after t_diff
                msd_over_time_in_linear_regime = msd_over_time_values[t_diff_determined_index:]
                time_in_linear_regime = time_values[t_diff_determined_index:]  # For reference or weighted average

                if msd_over_time_in_linear_regime.size > 1:  # Need some points to average
                    average_msd_over_time_for_D = np.mean(msd_over_time_in_linear_regime)
                    print(
                        f"Average MSD/Time in linear regime (for t >= {t_diff_determined_value:.4f}s): {average_msd_over_time_for_D:.4e}")

                    # For 2D diffusion, MSD/Time = 4D
                    diffusion_coefficient_D = average_msd_over_time_for_D / 4.0
                    print(f"Estimated Diffusion Coefficient (D) from linear regime: {diffusion_coefficient_D:.4e}")
                else:
                    print("Not enough data points found in the determined linear regime to calculate D robustly.")
            else:
                print("t_diff is at or near the end of the data; cannot calculate D from subsequent points.")
        else:
            print(
                f"Could not determine t_diff: Slope did not consistently fall below tolerance {slope_tolerance} after t={initial_transient_time_threshold}s.")
            print("Consider adjusting 'slope_tolerance' or 'initial_transient_time_threshold'.")
            print("You may need to visually inspect the slope plot to choose a t_diff manually or refine criteria.")
    else:
        print(
            f"No data points available after the initial transient time threshold of {initial_transient_time_threshold}s to determine t_diff.")

    # --- Plotting ---
    fig, ax1 = plt.subplots(figsize=(12, 7))

    color = 'tab:blue'
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('MSD/Time', color=color)
    ax1.plot(time_values, msd_over_time_values, color=color, label='MSD/Time', alpha=0.8)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, linestyle=':', axis='y', which='both')

    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('Slope of MSD/Time (Derivative)', color=color)
    ax2.plot(time_values, slope_of_msd_over_time, color=color, linestyle='--', label='Slope of MSD/Time', alpha=0.8)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.axhline(0, color='gray', lw=0.7, linestyle=':')
    ax2.axhline(slope_tolerance, color='gray', lw=0.5, linestyle=':', label=f'Tolerance (+/-{slope_tolerance:.0e})')
    ax2.axhline(-slope_tolerance, color='gray', lw=0.5, linestyle=':')

    # Annotate t_diff and average MSD/Time for D on the plot if found
    if t_diff_determined_value is not None:
        ax1.axvline(t_diff_determined_value, color='green', linestyle='-.', lw=1.5,
                    label=f'Est. t_diff ≈ {t_diff_determined_value:.2f}s')
        if not np.isnan(average_msd_over_time_for_D):
            # Plot the average line only over the region it was calculated for
            time_for_avg_line = time_values[t_diff_determined_index:]
            ax1.plot(time_for_avg_line, np.full_like(time_for_avg_line, average_msd_over_time_for_D),
                     color='purple', linestyle=':', lw=2,
                     label=f'Avg. MSD/Time (for D) ≈ {average_msd_over_time_for_D:.2e}')

    fig.suptitle('MSD/Time Analysis for Diffusion Coefficient', fontsize=16)

    # Collect all legend handles and labels
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper right')

    plt.show()

except FileNotFoundError:
    print(f"Error: The file '{filename}' was not found. Make sure it's in the correct path.")
except ValueError as ve:
    print(f"ValueError: {ve}. Please check data integrity or script logic.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")