import pandas as pd
import matplotlib.pyplot as plt
import numpy as np  # For handling potential NaNs
import os

# --- Define the input CSV files ---
# Add all the CSV files you want to merge into this list
csv_files_to_merge = [
    "truncated_gaussian_sweep_results_MEDIABASSADEVSTBASSA.csv",  #
    "truncated_gaussian_sweep_results_MEDIABASSADEVSTALTA.csv",  #
    # Add more file names here if you have them, e.g.:
     "truncated_gaussian_sweep_results_MEDIABASSADEVST09.csv",
     "truncated_gaussian_sweep_results_MEDIABASSADEVSTALTA2.csv"
]

# Output plot filename
output_plot_filename = "all_combined_sweep_summary_plots.png"

# List to hold all DataFrames
list_of_dataframes = []
all_files_found = True

# Load each CSV file
for f_name in csv_files_to_merge:
    if os.path.exists(f_name):
        try:
            df = pd.read_csv(f_name)
            if not df.empty:
                list_of_dataframes.append(df)
                print(f"Successfully loaded: {f_name} ({len(df)} rows)")
            else:
                print(f"Warning: File {f_name} is empty and will be skipped.")
        except pd.errors.EmptyDataError:
            print(f"Warning: File {f_name} is empty or invalid and will be skipped.")
        except Exception as e:
            print(f"Error loading {f_name}: {e}")
            all_files_found = False  # Mark that at least one file had issues
    else:
        print(f"Warning: File {f_name} not found and will be skipped.")
        all_files_found = False  # Mark that at least one file was missing

if not list_of_dataframes:
    print("No data loaded. Exiting.")
    exit()

# Merge all loaded DataFrames
combined_df = pd.concat(list_of_dataframes, ignore_index=True)
print(f"\nSuccessfully combined data from {len(list_of_dataframes)} file(s). Total data points: {len(combined_df)}")

# --- Plotting Final Results (Same plotting logic as before) ---
if not combined_df.empty:
    # Prepare data for Plot 1: D vs. mean_input_rest_prob (color by variance_input_rest_prob)
    plot_df1 = combined_df.dropna(subset=['mean_input_rest_prob', 'D', 'variance_input_rest_prob'])

    # Prepare data for Plot 2: t_diff vs. variance_input_rest_prob (color by mean_input_rest_prob)
    plot_df2 = combined_df.dropna(subset=['variance_input_rest_prob', 't_diff', 'mean_input_rest_prob'])

    plt.figure(figsize=(14, 6))

    # Plot 1: D vs. mean_input_rest_prob (color by variance_input_rest_prob)
    plt.subplot(1, 2, 1)
    if not plot_df1.empty:
        scatter1 = plt.scatter(plot_df1['mean_input_rest_prob'], plot_df1['D'],
                               c=plot_df1['variance_input_rest_prob'], cmap='viridis',
                               alpha=0.7, edgecolors='k', linewidths=0.5)
        plt.xlabel('Mean of Underlying Rest Prob Distribution')
        plt.ylabel('Estimated Diffusion Coefficient (D)')
        plt.title('D vs. Mean of Rest Prob (Color by Variance)')
        cbar1 = plt.colorbar(scatter1, label='Variance of Underlying Rest Prob')
        plt.grid(True, linestyle=':')
    else:
        plt.text(0.5, 0.5, 'No valid data for D vs Mean plot',
                 horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)

    # Plot 2: t_diff vs. variance_input_rest_prob (color by mean_input_rest_prob)
    plt.subplot(1, 2, 2)
    if not plot_df2.empty:
        scatter2 = plt.scatter(plot_df2['variance_input_rest_prob'], plot_df2['t_diff'],
                               c=plot_df2['mean_input_rest_prob'], cmap='plasma',
                               alpha=0.7, edgecolors='k', linewidths=0.5)
        plt.xlabel('Variance of Underlying Rest Prob Distribution (std_dev^2)')
        plt.ylabel('Estimated Convergence Time t_diff (s)')
        plt.title('t_diff vs. Variance of Rest Prob (Color by Mean)')
        cbar2 = plt.colorbar(scatter2, label='Mean of Underlying Rest Prob')
        plt.grid(True, linestyle=':')
    else:
        plt.text(0.5, 0.5, 'No valid data for t_diff vs Variance plot',
                 horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)

    plt.tight_layout()
    plt.savefig(output_plot_filename)
    print(f"\nCombined plot saved as {output_plot_filename}")
    plt.show()
else:
    print("No data collected (combined DataFrame is empty) to plot.")

if not all_files_found:
    print(
        "\nNote: Some CSV files were not found or were empty. The plot generated is based on successfully loaded files.")