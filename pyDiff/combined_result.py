import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import os
import statsmodels.api as sm

# --- CONTROL FLAGS ---
PERFORM_PARAMETRIC_FITTING = False  # Set to True for per-family exponential/linear fits
PERFORM_LOWESS_SMOOTHING = False  # Set to True for per-family LOWESS trend lines

FIT_PLOT1_GLOBAL_EXCL_CRIT1 = True  # Fit 1 for Plot 1
FIT_PLOT1_GLOBAL_EXCL_CRIT2 = True  # Fit 2 for Plot 1
FIT_GLOBAL_LINEAR_PLOT2 = False  # Control for global linear fit on Plot 2 (middle plot)
FIT_GLOBAL_LINEAR_PLOT3 =True  # Control for global linear fit on Plot 3 (new rightmost plot)
# --------------------

PERFORM_BINNED_V_FITTING_PLOT1 = False  # <<< New flag for this feature
NUM_V_BINS_PLOT1 = 2                   # <<< Number of bins for V




# --- NORMALIZATION PARAMETER ---
D_STAR_NORMALIZATION = 2.5 / 1000  # User-modifiable D* value for normalization.
# --------------------

# --- Exclusion Thresholds for Plot 1 Global Fits ---
PLOT1_EXCL_CRIT1_MEAN_GT = 0.3
PLOT1_EXCL_CRIT1_VAR_GT = 0.04

PLOT1_EXCL_CRIT2_MEAN_GT = 0.3
PLOT1_EXCL_CRIT2_VAR_LT = 0.02
# --------------------

if PERFORM_PARAMETRIC_FITTING:
    def model_exponential_decay(x, A, B): return A * np.exp(-B * x)


    def model_linear(x, A, B): return A + B * x


    def model_quadratic(x, A, B, C): return A + B * x + C * x ** 2

# --- Configuration ---
csv_files_to_merge = [
    "truncated_gaussian_sweep_results_ACTUAL_MEDIAALTADEVSTALTA1.csv",
    "truncated_gaussian_sweep_results_ACTUAL_MEDIAALTADEVSTALTA2.csv",
    "truncated_gaussian_sweep_results_ACTUAL_MEDIAALTADEVSTALTA3.csv",
    "truncated_gaussian_sweep_results_ACTUAL_MEDIAALTADEVSTALTA4.csv",
    "truncated_gaussian_sweep_results_ACTUAL_MEDIAALTADEVSTBASSA1.csv",
    "truncated_gaussian_sweep_results_ACTUAL_MEDIAALTADEVSTBASSA2.csv",
    "truncated_gaussian_sweep_results_ACTUAL_MEDIABASSADEVSTALTA1.csv",
    "var_D_sweep_results_PROVA.csv",  # This file contains var_D_late_time and the new variance components
    "var_D_sweep_results_PROVA2.csv"
]
output_plot_filename = "all_combined_summary_plots_corrected_v13.png"

# --- Load and Merge Data ---
list_of_dataframes = []
all_files_found_and_valid = True
print("Attempting to load and merge the following CSV files:")
for f_name in csv_files_to_merge:
    if not os.path.exists(f_name):
        print(f"  Warning: File '{f_name}' not found and will be skipped.")
        all_files_found_and_valid = False
        continue
    try:
        df = pd.read_csv(f_name)
        if not df.empty:
            # Standardize D_avg column
            if 'D' in df.columns and 'D_avg' not in df.columns:
                df.rename(columns={'D': 'D_avg'}, inplace=True)

            # Standardize actual_mean_prest column
            if 'actual_mean_rest_prob' in df.columns and 'actual_mean_prest' not in df.columns:
                df.rename(columns={'actual_mean_rest_prob': 'actual_mean_prest'}, inplace=True)

            # Standardize avg_within_landscape_var_prest (from older actual_variance_rest_prob if new one not present)
            if 'actual_variance_rest_prob' in df.columns and 'avg_within_landscape_var_prest' not in df.columns:
                df.rename(columns={'actual_variance_rest_prob': 'avg_within_landscape_var_prest'}, inplace=True)

            # Ensure all potentially needed columns exist, filling with NaN if necessary
            ensure_cols = ['actual_mean_prest', 'avg_within_landscape_var_prest',
                           'var_between_landscape_mean_prest', 'D_avg', 't_diff', 'var_D_late_time']
            for col in ensure_cols:
                if col not in df.columns:
                    df[col] = np.nan

            list_of_dataframes.append(df)
        else:
            print(f"  Warning: File {f_name} is empty.")
    except Exception as e:
        print(f"  Error loading {f_name}: {e}");
        all_files_found_and_valid = False

if not list_of_dataframes: print("\nNo data. Exiting."); exit()
combined_df = pd.concat(list_of_dataframes, ignore_index=True)
print(f"\nCombined data: {len(combined_df)} points.")
if not all_files_found_and_valid: print("Note: Some files missing/empty or had errors.")

# --- Create 'total_actual_var_prest' column in combined_df ---
# This will be the primary variance measure for x-axes where appropriate
if 'avg_within_landscape_var_prest' in combined_df.columns and \
        'var_between_landscape_mean_prest' in combined_df.columns:
    # Use fillna(0) for components if you want a sum even if one part is missing (treats missing as 0)
    # Or, let NaNs propagate by removing fillna(0):
    # combined_df['total_actual_var_prest'] = combined_df['avg_within_landscape_var_prest'] + combined_df['var_between_landscape_mean_prest']
    combined_df['total_actual_var_prest'] = combined_df['avg_within_landscape_var_prest'].fillna(0) + \
                                            combined_df['var_between_landscape_mean_prest'].fillna(0)
    print(
        "Calculated 'total_actual_var_prest' using sum of 'avg_within_landscape_var_prest' and 'var_between_landscape_mean_prest'.")
else:
    print(
        "Warning: Source columns for 'total_actual_var_prest' ('avg_within_landscape_var_prest', 'var_between_landscape_mean_prest') not found.")
    # Fallback to using 'avg_within_landscape_var_prest' if 'var_between_landscape_mean_prest' is missing
    # (This covers cases where CSVs might only have the old 'actual_variance_rest_prob' renamed to 'avg_within_landscape_var_prest')
    if 'avg_within_landscape_var_prest' in combined_df.columns:
        combined_df['total_actual_var_prest'] = combined_df['avg_within_landscape_var_prest']
        print("Using 'avg_within_landscape_var_prest' as 'total_actual_var_prest'.")
    else:
        combined_df['total_actual_var_prest'] = np.nan
        print("'total_actual_var_prest' set to NaN as essential variance components are missing.")

# --- Plotting and Trend Extrapolation ---
if combined_df.empty: print("No data to plot."); exit()
plt.figure(figsize=(27, 7))  # Width increased for 3 plots

effective_d_star = D_STAR_NORMALIZATION if D_STAR_NORMALIZATION > 0 else 1.0







# --- Plot 1: D_avg/D* vs. actual_mean_prest ---
ax1 = plt.subplot(1, 3, 1)
plot_df1_base = combined_df.dropna(subset=['actual_mean_prest', 'D_avg', 'total_actual_var_prest'])

if not plot_df1_base.empty:
    plot_df1_base = plot_df1_base.copy()
    plot_df1_base.loc[:, 'D_normalized'] = plot_df1_base['D_avg'] / effective_d_star
else:
    # Ensure column exists even if df is empty to prevent potential downstream errors
    # if code expects the column on an empty DataFrame.
    plot_df1_base['D_normalized'] = pd.Series(dtype=float)


    # ++++ START: DEVIATION CALCULATIONS FOR PLOT 1 ++++
    if 'actual_mean_prest' in plot_df1_base.columns and not plot_df1_base['D_normalized'].isnull().all():
        plot_df1_base.loc[:, 'D_theory_plot1'] = 1 - plot_df1_base['actual_mean_prest']
        plot_df1_base.loc[:, 'residuals_D_plot1'] = plot_df1_base['D_normalized'] - plot_df1_base['D_theory_plot1']

        # Ensure residuals are not all NaN before calculating stats
        valid_residuals_p1 = plot_df1_base['residuals_D_plot1'].dropna()
        if not valid_residuals_p1.empty:
            mae_D1 = valid_residuals_p1.abs().mean()
            mse_D1 = (valid_residuals_p1**2).mean()
            rmse_D1 = np.sqrt(mse_D1)
            mean_residual_D1 = valid_residuals_p1.mean()
            std_residual_D1 = valid_residuals_p1.std()

            print("\n--- Deviations from D_avg/D* = 1 - actual_mean_prest (Plot 1 Data) ---")
            print(f"  Number of data points for D_avg analysis: {len(valid_residuals_p1)}")
            print(f"  Mean Absolute Error (MAE): {mae_D1:.4e}")
            print(f"  Mean Squared Error (MSE): {mse_D1:.4e}")
            print(f"  Root Mean Squared Error (RMSE): {rmse_D1:.4e}")
            print(f"  Mean of Residuals: {mean_residual_D1:.4e} (ideally close to 0)")
            print(f"  Std Dev of Residuals: {std_residual_D1:.4e}")
        else:
            print("\n--- Deviations for Plot 1 (D_avg/D*): No valid residual data for statistics. ---")
    else:
        print("\n--- Deviations for Plot 1 (D_avg/D*): Not enough data or D_normalized is all NaN. ---")
    # ++++ END: DEVIATION CALCULATIONS FOR PLOT 1 ++++


ax1.set_xlabel('Actual Mean of Rest Prob ($E[\\bar{\\rho}_i]$)')
if effective_d_star == 1.0:
    ax1.set_ylabel('Avg. Diffusion Coefficient ($D_{avg}$)')
else:
    ax1.set_ylabel(f'Normalized Avg. Diff. Coeff. ($D / D^*$) [$D^*={effective_d_star:.1e}$]')
ax1.grid(True, linestyle=':')
ax1.set_title('$D_{avg}/D^*$ vs. Actual Mean')

legend_handles_plot1 = []
legend_labels_plot1 = []

if not plot_df1_base.empty:
    vmin_var_p1 = plot_df1_base['total_actual_var_prest'].min()
    vmax_var_p1 = plot_df1_base['total_actual_var_prest'].max()
    main_scatter1 = ax1.scatter(plot_df1_base['actual_mean_prest'], plot_df1_base['D_normalized'],
                                c=plot_df1_base['total_actual_var_prest'], cmap='viridis',
                                vmin=vmin_var_p1, vmax=vmax_var_p1,
                                alpha=0.9, edgecolors='grey', linewidths=0.5, s=50, zorder=3, label="Data")
    if main_scatter1 not in legend_handles_plot1: legend_handles_plot1.append(main_scatter1)
    if "Data" not in legend_labels_plot1: legend_labels_plot1.append("Data")
    cbar1 = plt.colorbar(main_scatter1, ax=ax1, label='Total Actual Var Rest Prob' )

    # Fit 1 for Plot 1
    if FIT_PLOT1_GLOBAL_EXCL_CRIT1 and len(plot_df1_base) >= 2:
        df_for_fit = plot_df1_base.copy()
        exclusion_condition = (df_for_fit['actual_mean_prest'] > PLOT1_EXCL_CRIT1_MEAN_GT) & \
                              (df_for_fit['total_actual_var_prest'] > PLOT1_EXCL_CRIT1_VAR_GT)
        num_excluded = exclusion_condition.sum()
        df_for_fit = df_for_fit[~exclusion_condition]
        if len(df_for_fit) >= 2:
            try:
                x_fit_data = df_for_fit['actual_mean_prest'].values
                y_fit_data = df_for_fit['D_normalized'].values
                slope, intercept = np.polyfit(x_fit_data, y_fit_data, 1)
                x_line = np.array([x_fit_data.min(), x_fit_data.max()])
                y_line = slope * x_line + intercept
                label = (f'Fit 1 (Excl. M>{PLOT1_EXCL_CRIT1_MEAN_GT}, V>{PLOT1_EXCL_CRIT1_VAR_GT})\n'
                         f'y={slope:.2e}x+{intercept:.2e}')
                line, = ax1.plot(x_line, y_line, color='magenta', linestyle='-.', linewidth=2.5, zorder=2)
                if line not in legend_handles_plot1: legend_handles_plot1.append(line)
                if label not in legend_labels_plot1: legend_labels_plot1.append(label)
            except Exception as e:
                print(f"Plot 1 Global Fit (Crit 1) failed: {e}")
        elif num_excluded > 0:
            print(f"Not enough data for Plot 1 Global Fit (Crit 1) after {num_excluded} exclusions.")

    # Fit 2 for Plot 1
    if FIT_PLOT1_GLOBAL_EXCL_CRIT2 and len(plot_df1_base) >= 2:
        df_for_fit = plot_df1_base.copy()
        exclusion_condition = (df_for_fit['actual_mean_prest'] > PLOT1_EXCL_CRIT2_MEAN_GT) & \
                              (df_for_fit['total_actual_var_prest'] < PLOT1_EXCL_CRIT2_VAR_LT)
        num_excluded = exclusion_condition.sum()
        df_for_fit = df_for_fit[~exclusion_condition]
        if len(df_for_fit) >= 2:
            try:
                x_fit_data = df_for_fit['actual_mean_prest'].values
                y_fit_data = df_for_fit['D_normalized'].values
                slope, intercept = np.polyfit(x_fit_data, y_fit_data, 1)
                x_line = np.array([x_fit_data.min(), x_fit_data.max()])
                y_line = slope * x_line + intercept
                label = (f'Fit 2 (Excl. M>{PLOT1_EXCL_CRIT2_MEAN_GT}, V<{PLOT1_EXCL_CRIT2_VAR_LT})\n'
                         f'y={slope:.2e}x+{intercept:.2e}')
                line, = ax1.plot(x_line, y_line, color='red', linestyle='--', linewidth=2.5, zorder=2)
                if line not in legend_handles_plot1: legend_handles_plot1.append(line)
                if label not in legend_labels_plot1: legend_labels_plot1.append(label)
            except Exception as e:
                print(f"Plot 1 Global Fit (Crit 2) failed: {e}")
        elif num_excluded > 0:
            print(f"Not enough data for Plot 1 Global Fit (Crit 2) after {num_excluded} exclusions.")

    if PERFORM_LOWESS_SMOOTHING or PERFORM_PARAMETRIC_FITTING:  # Per-family fits
        unique_variances = sorted(plot_df1_base['total_actual_var_prest'].dropna().round(4).unique())
        added_lowess_legend1 = False
        added_param_legend1 = False
        for var_idx, var_val in enumerate(unique_variances):
            family_df = plot_df1_base[np.isclose(plot_df1_base['total_actual_var_prest'].round(4), var_val)]
            if len(family_df) < 3: continue
            x_data = family_df['actual_mean_prest'].values
            y_data = family_df['D_normalized'].values
            sort_indices = np.argsort(x_data)
            x_data_sorted, y_data_sorted = x_data[sort_indices], y_data[sort_indices]
            # fit_line_color = cmap_fits1(norm_fits1(var_val)) # Requires cmap_fits1, norm_fits1 from main_scatter1
            # ... (rest of LOWESS/Parametric fitting as in your full script) ...

    if legend_handles_plot1:
        ax1.legend(handles=legend_handles_plot1, labels=legend_labels_plot1, title="Trend Types", fontsize='small',
                   loc='best')

        # --- Binned Fitting by V (Total Actual Variance of Rest Prob) for Plot 1 ---
        if PERFORM_BINNED_V_FITTING_PLOT1 and len(plot_df1_base) > NUM_V_BINS_PLOT1:  # Ensure enough data for binning
            from scipy.optimize import curve_fit


            def model_for_binned_fit(x_val, K_eff_bin_val):
                """Model function: y = (1-x) - K_eff_bin / (1-x)"""
                # Ensure 1-x_val is not zero; x_val should be < 1
                # This should generally be true if x_val is a probability < 1
                denominator = 1.0 - np.asarray(x_val)
                # Avoid division by zero or very small numbers by adding a small epsilon
                # or by ensuring x_val is strictly less than 1 before calling.
                # For now, assuming x_val < 1 from context.
                return (1.0 - np.asarray(x_val)) - K_eff_bin_val / denominator


            print(
                "\n--- Performing Binned Fitting for Plot 1 (D_avg/D* vs. Mean Rest Prob, binned by Var Rest Prob) ---")
            print(f"Target function: y = (1-x) - K_eff_bin / (1-x)")

            fitted_K_eff_results = []  # To store (V_bin_mean, K_eff_fitted)

            try:
                # Create bins for 'total_actual_var_prest' (V)
                # Using qcut for (roughly) equal number of points per bin
                plot_df1_base.loc[:, 'V_bin_category'] = pd.qcut(
                    plot_df1_base['total_actual_var_prest'],
                    q=NUM_V_BINS_PLOT1,
                    labels=False,  # Get integer labels for bins
                    duplicates='drop'  # If not enough unique values to form q bins, drop some
                )

                bin_colors = plt.cm.get_cmap('cool', NUM_V_BINS_PLOT1 + 2)  # Get some distinct colors

                for bin_idx in sorted(plot_df1_base['V_bin_category'].unique()):
                    bin_df = plot_df1_base[plot_df1_base['V_bin_category'] == bin_idx]

                    if len(bin_df) < 2:  # Need at least 2 points for a fit (ideally more for K_eff)
                        print(f"  Skipping V bin {bin_idx}: only {len(bin_df)} data point(s).")
                        continue

                    V_bin_mean = bin_df['total_actual_var_prest'].mean()
                    V_bin_min = bin_df['total_actual_var_prest'].min()
                    V_bin_max = bin_df['total_actual_var_prest'].max()

                    x_bin_data = bin_df['actual_mean_prest'].values
                    y_bin_data = bin_df['D_normalized'].values

                    # Ensure x_bin_data values are less than 1 for the model
                    if np.any(x_bin_data >= 1.0):
                        print(f"  Warning: Skipping V bin {bin_idx} (V_mean={V_bin_mean:.3e}) due to x_data >= 1.")
                        continue

                    try:
                        # Initial guess for K_eff_bin: K_eff ~ K_global * V_bin_mean.
                        # If we assume K_global is small, say 0.1 to 1
                        initial_K_eff_guess = 0.1 * V_bin_mean if V_bin_mean > 0 else 0.01
                        popt, pcov = curve_fit(model_for_binned_fit, x_bin_data, y_bin_data, p0=[initial_K_eff_guess],
                                               maxfev=5000)
                        K_eff_fitted = popt[0]

                        fitted_K_eff_results.append(
                            {'V_bin_mean': V_bin_mean, 'K_eff_fitted': K_eff_fitted, 'num_points': len(bin_df)})
                        print(
                            f"  Bin {bin_idx} (V range [{V_bin_min:.2e}, {V_bin_max:.2e}], V_mean={V_bin_mean:.3e}, N={len(bin_df)}): Fitted K_eff_bin = {K_eff_fitted:.3e}")

                        # Plot the fitted curve for this bin
                        x_line_bin = np.linspace(min(0, x_bin_data.min()), x_bin_data.max(), 100)
                        # Ensure x_line_bin does not hit 1.0 exactly for plotting
                        x_line_bin = np.clip(x_line_bin, None, 0.9999)
                        y_line_bin = model_for_binned_fit(x_line_bin, K_eff_fitted)

                        bin_label_text = f'Fit (V_mean≈{V_bin_mean:.2e}, K_eff≈{K_eff_fitted:.1e})'
                        line_bin, = ax1.plot(x_line_bin, y_line_bin, linestyle=':', linewidth=2.0,
                                             color=bin_colors(bin_idx), label=bin_label_text, zorder=1)
                        if line_bin not in legend_handles_plot1: legend_handles_plot1.append(line_bin)
                        if bin_label_text not in legend_labels_plot1: legend_labels_plot1.append(bin_label_text)

                    except RuntimeError:
                        print(f"  Warning: Fit did not converge for V bin {bin_idx} (V_mean={V_bin_mean:.3e}).")
                    except Exception as e:
                        print(f"  Error fitting V bin {bin_idx} (V_mean={V_bin_mean:.3e}): {e}")

                # Re-apply legend for plot 1 if new items were added
                if legend_handles_plot1:  # Check if list is not empty
                    ax1.legend(handles=legend_handles_plot1, labels=legend_labels_plot1, title="Trend Types",
                               fontsize='small', loc='best')


            except ValueError as ve:
                print(
                    f"Warning: Could not create quantile bins for V ('total_actual_var_prest'). Min/Max: {plot_df1_base['total_actual_var_prest'].min()}/{plot_df1_base['total_actual_var_prest'].max()}. Error: {ve}")
                print("Skipping binned V fitting for Plot 1.")
            except ImportError:
                print("Warning: scipy.optimize.curve_fit not available. Skipping binned V fitting.")

            # --- Analysis of K_eff_bin vs V_bin_mean ---
            if fitted_K_eff_results:
                print("\n--- Analysis of K_eff_bin vs. V_bin_mean ---")
                V_means_for_K_fit = np.array([res['V_bin_mean'] for res in fitted_K_eff_results])
                K_effs_for_K_fit = np.array([res['K_eff_fitted'] for res in fitted_K_eff_results])

                if len(V_means_for_K_fit) >= 2:
                    # Fit K_eff = K * V (linear through origin)
                    # Forcing through origin: y = a*x. Reshape V_means to be (N,1)
                    K_global_slope_origin, residuals, rank, singular_values = np.linalg.lstsq(
                        V_means_for_K_fit[:, np.newaxis], K_effs_for_K_fit, rcond=None)
                    K_global_origin = K_global_slope_origin[0]
                    print(f"  Fit K_eff_bin = K * V_bin_mean (through origin): Global K ≈ {K_global_origin:.3f}")

                    # Fit K_eff = K * V + K0 (with intercept)
                    K_global_slope_intercept, K_global_intercept = np.polyfit(V_means_for_K_fit, K_effs_for_K_fit, 1)
                    print(
                        f"  Fit K_eff_bin = K * V_bin_mean + K0:          Global K ≈ {K_global_slope_intercept:.3f}, Intercept K0 ≈ {K_global_intercept:.3f}")

                    # You can create a new plot for K_eff_bin vs V_bin_mean if desired
                    # plt.figure(figsize=(7,5))
                    # plt.plot(V_means_for_K_fit, K_effs_for_K_fit, 'o', label='Fitted $K_{eff,bin}$ values')
                    # plt.plot(V_means_for_K_fit, K_global_origin * V_means_for_K_fit, '--', label=f'$K_{{eff}} = {K_global_origin:.2f} \\times V_{{mean}}$')
                    # plt.plot(V_means_for_K_fit, K_global_slope_intercept * V_means_for_K_fit + K_global_intercept, ':', label=f'$K_{{eff}} = {K_global_slope_intercept:.2f} \\times V_{{mean}} + {K_global_intercept:.2f}$')
                    # plt.xlabel('Mean V in Bin ($V_{bin\_mean}$)')
                    # plt.ylabel('Fitted $K_{eff\_bin}$')
                    # plt.title('$K_{eff\_bin}$ vs. $V_{bin\_mean}$ and Global K Estimation')
                    # plt.legend()
                    # plt.grid(True)
                    # if show_plots_main_run: # Using your existing flag
                    #     plt.show()

                else:
                    print("  Not enough successfully fitted bins to analyze K_eff_bin vs. V_bin_mean.")
        # --- End of Binned Fitting Section ---


else:
    ax1.text(0.5, 0.5, 'No valid data for $D_{avg}/D^*$ vs Actual Mean plot', transform=ax1.transAxes, ha='center',
             va='center')







# --- Plot 2: t_diff vs. total_actual_var_prest ---
ax2 = plt.subplot(1, 3, 2)
plot_df2_base = combined_df.dropna(subset=['total_actual_var_prest', 't_diff', 'actual_mean_prest'])
ax2.set_xlabel('Total Actual Variance of Rest Prob ')
ax2.set_ylabel('Estimated Convergence Time $t_{diff}$ (s)')
ax2.grid(True, linestyle=':')
ax2.set_title('$t_{diff}$ vs. ')

legend_handles_plot2 = []
legend_labels_plot2 = []

if not plot_df2_base.empty:
    vmin_mean_p2 = plot_df2_base['actual_mean_prest'].min()
    vmax_mean_p2 = plot_df2_base['actual_mean_prest'].max()
    main_scatter2 = ax2.scatter(plot_df2_base['total_actual_var_prest'], plot_df2_base['t_diff'],
                                c=plot_df2_base['actual_mean_prest'], cmap='viridis',
                                vmin=vmin_mean_p2, vmax=vmax_mean_p2,
                                alpha=0.9, edgecolors='grey', linewidths=0.5, s=50, zorder=3, label="Data")
    if main_scatter2 not in legend_handles_plot2: legend_handles_plot2.append(main_scatter2)
    if "Data" not in legend_labels_plot2: legend_labels_plot2.append("Data")
    cbar2 = plt.colorbar(main_scatter2, ax=ax2, label='Actual Mean of Rest Prob ($E[\\bar{\\rho}_i]$)')

    if FIT_GLOBAL_LINEAR_PLOT2 and len(plot_df2_base) >= 2:
        try:
            x_all_p2 = plot_df2_base['total_actual_var_prest'].values
            y_all_p2 = plot_df2_base['t_diff'].values
            slope, intercept = np.polyfit(x_all_p2, y_all_p2, 1)
            x_line_p2 = np.array([x_all_p2.min(), x_all_p2.max()])
            y_line_p2 = slope * x_line_p2 + intercept
            line_global2, = ax2.plot(x_line_p2, y_line_p2, color='blue', linestyle=':', linewidth=2.5, zorder=2)
            if line_global2 not in legend_handles_plot2: legend_handles_plot2.append(line_global2)
            label_fit2 = f'Global Linear Fit (y={slope:.2e}x+{intercept:.2e})'
            if label_fit2 not in legend_labels_plot2: legend_labels_plot2.append(label_fit2)
            print(f"Plot 2 Global Fit: t_diff = {slope:.3e}*TotalVar + {intercept:.3e}")
        except Exception as e:
            print(f"Plot 2 Global Fit failed: {e}")

    if legend_handles_plot2:
        ax2.legend(handles=legend_handles_plot2, labels=legend_labels_plot2, title="Trend Types", fontsize='small',
                   loc='best')
else:
    ax2.text(0.5, 0.5, 'No valid data for $t_{diff}$ vs Total Actual Variance plot', transform=ax2.transAxes,
             ha='center', va='center')

# --- Plot 3: Var(D)_late_time / (D*)^2 vs. total_actual_var_prest ---
ax3 = plt.subplot(1, 3, 3)
plot_df3_base = combined_df.dropna(subset=['var_D_late_time', 'total_actual_var_prest', 'actual_mean_prest'])

if not plot_df3_base.empty:
    plot_df3_base = plot_df3_base.copy()
    plot_df3_base.loc[:, 'var_D_normalized_sq'] = plot_df3_base['var_D_late_time'] / (effective_d_star ** 2)
else:
    if 'var_D_normalized_sq' not in plot_df3_base.columns: plot_df3_base['var_D_normalized_sq'] = pd.Series(dtype=float)

    # ++++ START: DEVIATION CALCULATIONS FOR PLOT 3 ++++
    if 'total_actual_var_prest' in plot_df3_base.columns and not plot_df3_base['var_D_normalized_sq'].isnull().all():
        plot_df3_base.loc[:, 'VarD_theory_plot3'] = plot_df3_base['total_actual_var_prest']  # Theory: Y = X
        plot_df3_base.loc[:, 'residuals_VarD_plot3'] = plot_df3_base['var_D_normalized_sq'] - plot_df3_base[
            'VarD_theory_plot3']

        valid_residuals_p3 = plot_df3_base['residuals_VarD_plot3'].dropna()
        if not valid_residuals_p3.empty:
            mae_VarD3 = valid_residuals_p3.abs().mean()
            mse_VarD3 = (valid_residuals_p3 ** 2).mean()
            rmse_VarD3 = np.sqrt(mse_VarD3)
            mean_residual_VarD3 = valid_residuals_p3.mean()
            std_residual_VarD3 = valid_residuals_p3.std()

            print("\n--- Deviations from Var(D)/(D*)^2 = total_actual_var_prest (Plot 3 Data) ---")
            print(f"  Number of data points for Var(D) analysis: {len(valid_residuals_p3)}")
            print(f"  Mean Absolute Error (MAE): {mae_VarD3:.4e}")
            print(f"  Mean Squared Error (MSE): {mse_VarD3:.4e}")
            print(f"  Root Mean Squared Error (RMSE): {rmse_VarD3:.4e}")
            print(f"  Mean of Residuals: {mean_residual_VarD3:.4e} (measures avg offset from Y=X line)")
            print(f"  Std Dev of Residuals: {std_residual_VarD3:.4e} (measures scatter around Y=X line after offset)")
        else:
            print("\n--- Deviations for Plot 3 (Var(D)): No valid residual data for statistics. ---")
    else:
        print("\n--- Deviations for Plot 3 (Var(D)): Not enough data or var_D_normalized_sq is all NaN. ---")
    # ++++ END: DEVIATION CALCULATIONS FOR PLOT 3 ++++






ax3.set_xlabel('Total Actual Variance of Rest Prob')
ax3.grid(True, linestyle=':')
ax3.set_title('$Var(D)_{late}/(D^*)^2$')
if effective_d_star == 1.0:
    ax3.set_ylabel('$Var(D)_{late}$')
else:
    ax3.set_ylabel(f'Normalized $Var(D)$ ($Var(D) / (D^*)^2$) [$D^*={effective_d_star:.1e}$]')

legend_handles_plot3 = []
legend_labels_plot3 = []

if not plot_df3_base.empty:
    vmin_mean_p3 = plot_df3_base['actual_mean_prest'].min()
    vmax_mean_p3 = plot_df3_base['actual_mean_prest'].max()

    main_scatter3 = ax3.scatter(plot_df3_base['total_actual_var_prest'], plot_df3_base['var_D_normalized_sq'],
                                c=plot_df3_base['actual_mean_prest'], cmap='coolwarm',
                                vmin=vmin_mean_p3, vmax=vmax_mean_p3,
                                alpha=0.9, edgecolors='black', linewidths=0.5, s=50, zorder=3, label="Data Points")
    if main_scatter3 not in legend_handles_plot3: legend_handles_plot3.append(main_scatter3)
    if "Data Points" not in legend_labels_plot3: legend_labels_plot3.append("Data Points")
    cbar3 = plt.colorbar(main_scatter3, ax=ax3, label='Actual Mean of Rest Prob ($E[\\bar{\\rho}_i]$)')

    if FIT_GLOBAL_LINEAR_PLOT3 and len(plot_df3_base) >= 2:
        try:
            x_all_p3 = plot_df3_base['total_actual_var_prest'].values
            y_all_p3 = plot_df3_base['var_D_normalized_sq'].values

            valid_fit_indices_p3 = ~np.isnan(x_all_p3) & ~np.isnan(y_all_p3)
            if np.sum(valid_fit_indices_p3) < 2:
                print("Not enough valid (non-NaN) points for Plot 3 Global Fit after final check.")
            else:
                x_fit_plot_p3 = x_all_p3[valid_fit_indices_p3]
                y_fit_plot_p3 = y_all_p3[valid_fit_indices_p3]  # Use only valid y for fitting too
                slope, intercept = np.polyfit(x_fit_plot_p3, y_fit_plot_p3, 1)

                x_line_p3 = np.array([x_fit_plot_p3.min(), x_fit_plot_p3.max()])
                y_line_p3 = slope * x_line_p3 + intercept

                line_global3, = ax3.plot(x_line_p3, y_line_p3, color='green', linestyle=':', linewidth=2.5, zorder=2)
                if line_global3 not in legend_handles_plot3: legend_handles_plot3.append(line_global3)
                label_fit3 = f'Global Linear Fit (y={slope:.2e}x+{intercept:.2e})'
                if label_fit3 not in legend_labels_plot3: legend_labels_plot3.append(label_fit3)

                y_axis_label_fit_p3 = "Var(D)/(D*^2)" if effective_d_star != 1.0 else "Var(D)"
                print(f"Plot 3 Global Fit: {y_axis_label_fit_p3} = {slope:.3e}*TotalActualVar + {intercept:.3e}")
        except Exception as e:
            print(f"Plot 3 Global Fit failed: {e}")

    if legend_handles_plot3:
        ax3.legend(handles=legend_handles_plot3, labels=legend_labels_plot3, title="Trend Types", fontsize='small',
                   loc='best')
else:
    ax3.text(0.5, 0.5, 'No data for $Var(D)_{late}/(D^*)^2$ vs Total Actual Var plot', transform=ax3.transAxes,
             ha='center', va='center')

plt.tight_layout()
plt.savefig(output_plot_filename)
print(f"\nCombined plot saved as {output_plot_filename}")
plt.show()