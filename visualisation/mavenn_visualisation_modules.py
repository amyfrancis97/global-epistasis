# Import modules
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.interpolate import make_interp_spline
import re 
import glob
import os
import math
from matplotlib.gridspec import GridSpec

# Create a dictionary for colors (optional)
colour_palette = sns.color_palette("colorblind")
def plot_affinity_distribution(trainval_df, out_dir, colour_palette, wild_type_line=False, test_df=pd.DataFrame()):
    # Create a 2x2 grid of plots
    fig, axs = plt.subplots(2, 2, figsize=(14, 12))

    # Ensure the colour_palette has enough colors
    if len(colour_palette) < 4:
        raise ValueError("The colour_palette must contain at least 4 colors.")

    # Plotting helper function
    def plot_single_distribution(ax, df, title, color, wild_type_line=False):
        # Plot the KDE curve on the provided axes
        sns.kdeplot(df['y'], ax=ax, color=color, fill=True, alpha=0.1, linewidth=2)
        ax.set_xlabel('Absolut! Binding Affinity', fontsize=18)
        ax.set_ylabel('Density', fontsize=18)
        ax.set_xlim(-105, -55)
        ax.set_ylim(0, 0.06)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
   
        # Optional: add wild type line
        if wild_type_line:
            ax.axvline(x=-86.05, color='black', linestyle='--', linewidth=1.5)
            ax.text(-102, ax.get_ylim()[1] * 0.9, 'Wild Type\n-86.05', color='black', 
                    verticalalignment='top', fontsize=20, ha='left')

        # Set the title for each subplot
        ax.set_title(title, fontsize=18, weight='bold')

    # Combine trainval_df and test_df for "All Mutants" if test_df is not empty
    df_all = pd.concat([trainval_df, test_df]).reset_index(drop=True) if not test_df.empty else trainval_df

    # Subplot 1: All Mutants
    plot_single_distribution(axs[0, 0], df_all, "All Mutants", colour_palette[0], wild_type_line)

    # Subplot 2: Single Mutants
    df_single = pd.concat([trainval_df[trainval_df['dist'] == 1], test_df[test_df['dist'] == 1]]).reset_index(drop=True) if not test_df.empty else trainval_df[trainval_df['dist'] == 1]
    plot_single_distribution(axs[0, 1], df_single, "Single Mutants", colour_palette[1], wild_type_line)

    # Subplot 3: Double Mutants
    df_double = pd.concat([trainval_df[trainval_df['dist'] == 2], test_df[test_df['dist'] == 2]]).reset_index(drop=True) if not test_df.empty else trainval_df[trainval_df['dist'] == 2]
    plot_single_distribution(axs[1, 0], df_double, "Double Mutants", colour_palette[2], wild_type_line)

    # Subplot 4: Triple Mutants
    df_triple = pd.concat([trainval_df[trainval_df['dist'] == 3], test_df[test_df['dist'] == 3]]).reset_index(drop=True) if not test_df.empty else trainval_df[trainval_df['dist'] == 3]
    plot_single_distribution(axs[1, 1], df_triple, "Triple Mutants", colour_palette[3], wild_type_line)

    # Adjust layout for better spacing
    plt.tight_layout(pad=3.0)

    # Save the combined plot
    plt.savefig(f"{out_dir}/affinity_distributions_grid_pretty_palette.png", dpi=300)

    # Show the plot
    plt.show()

def plot_measurement_process(model_dict, test_df):

    model_names_filtered = [f'1FBI_X_models_additive',
                            f'1FBI_X_models_neighbor',
                            f'1FBI_X_models_pairwise',
                            f'1FBI_X_models_blackbox']
    
    # Create figure and axes for plotting
    fig, axs = plt.subplots(2, 2, figsize=[10, 10])
    axs = axs.ravel()

    # Loop over models
    for ax, name in zip(axs, model_names_filtered):
        # Get model
        model = model_dict[name]
        
        y_test = test_df['y']
        x_test = test_df['x']    

        # Compute phi on test data
        phi_test = model.x_to_phi(x_test)

        # Set phi limits and create a grid in phi space
        phi_lim = [min(phi_test) - .5, max(phi_test) + .5]
        phi_grid = np.linspace(phi_lim[0], phi_lim[1], 1000)

        # Compute yhat for each phi gridpoint
        yhat_grid = model.phi_to_yhat(phi_grid)

        # Compute 95% CI for each yhat
        q = [0.025, 0.975]
        yqs_grid = model.yhat_to_yq(yhat_grid, q=q)

        # Plot 95% confidence interval
        ax.fill_between(phi_grid, yqs_grid[:, 0], yqs_grid[:, 1],
                        alpha=0.2, color='C1', lw=0, label='95% CI')

        # Plot GE nonlinearity
        ax.plot(phi_grid, yhat_grid,
                linewidth=3, color='C1', label='nonlinearity')

        # Plot scatter of phi and y values
        ax.scatter(phi_test, y_test,
                    color='C0', s=10, alpha=.3, label='test data', zorder=+100)

        # Style plot
        ax.set_xlim(phi_lim)
        ax.set_xlabel('latent phenotype ($\phi$)')
        ax.set_ylabel('measurement ($y$)')
        ax.set_title(f'name', fontsize=15)
        ax.legend(loc='lower right')

    # Adjust layout to prevent overlap
    fig.tight_layout()

    # Show the plot
    plt.show()

def plot_I_var_pred(info_df, name):
    print(name)

    # Filter the dataframe based on the full name pattern
    filtered_df = info_df[info_df['name'] == name].reset_index(drop=True)
    print(filtered_df)

    # Create figure
    fig, ax = plt.subplots(figsize=[8, 4])

    # Plot bars with seaborn
    barplot = sns.barplot(ax=ax,
                          data=filtered_df,
                          hue='metric',
                          x='gpmap',
                          y='I',
                          ci=None)  # Disable internal confidence intervals since we'll use custom error bars

    # Add error bars to the plot
    for i, bar in enumerate(barplot.patches):
        if i < len(filtered_df):
            # Get the center of the bar
            x = bar.get_x() + bar.get_width() / 2
            y = bar.get_height()
            
            # Add error bars
            ax.errorbar(x=x, 
                        y=y, 
                        yerr=filtered_df['dI'].iloc[i],  # Match the correct error with the bar
                        fmt='none', 
                        color='black', 
                        capsize=3, 
                        elinewidth=1, 
                        capthick=1)

    # Set axis labels and limits
    ax.set_ylabel('Information (bits)')
    ax.set_xlabel('')
    ax.set_ylim([-1, 2])

    # Place the legend
    ax.legend(loc='upper left')

    # Show the plot
    plt.show()

def plot_metrics_double_triple(metric_df, metric_type, colour_palette):
    # Set up the figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Define bar width and positions
    bar_width = 0.2
    index = np.arange(len(metric_df['gpmap_type'].unique()))  # One bar per unique 'gpmap_type'

    # Define the unique model types and test data types
    model_types = metric_df['gpmap_type'].unique()
    test_data_types = metric_df['test_data_type'].unique()

    # Plot the bars within each model type
    for i, test_data_type in enumerate(test_data_types):
        # Extract data for each test_data_type
        subset = metric_df[metric_df['test_data_type'] == test_data_type]
        
        # Align bars using index + i * bar_width to create grouped bars
        ax.bar(index + i * bar_width, subset[metric_type], bar_width, label=test_data_type, color = colour_palette[i], edgecolor = "black")

    # Set labels, title, and tick positions
    ax.set_xlabel('Model Type')
    ax.set_ylabel(metric_type)
    ax.set_title(f'{metric_type.capitalize()} for different models and test data types')
    ax.set_xticks(index + bar_width)
    ax.set_xticklabels(model_types)
    ax.legend(title="Test Data Type")

    # Display the plot
    plt.tight_layout()
    plt.show()

def plot_sample_size_metrics(results, out_dir, sample = False, palette = colour_palette):
    # Define the list of models to compare and metrics to plot
    models = ['additive', 'neighbor', 'pairwise', 'blackbox']
    metrics = ['r2', 'spearman']

    # Create a 2x2 grid of subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Flatten the axes array for easy iteration
    axes = axes.flatten()

    # Loop over each metric to create a subplot
    for j, metric in enumerate(metrics):
        ax = axes[j]
        
        # Loop over each model to plot them on the same graph
        for i, model in enumerate(models):
            if sample:
                model_data = results.loc[results['gpmap_type'] == model, metric].sample(sample)
            else:
                model_data = results.loc[results['gpmap_type'] == model, metric]
            sample_size = results.loc[results['gpmap_type'] == model, 'sample_size'][model_data.index]

            # Sort the data by sample size
            sorted_indices = np.argsort(sample_size)
            sample_size_sorted = sample_size.iloc[sorted_indices]
            model_data_sorted = model_data.iloc[sorted_indices]

            # Smooth the curve using SciPy's spline interpolation
            xnew = np.linspace(sample_size_sorted.min(), sample_size_sorted.max(), 300)  # Create more points for smoothness
            spl = make_interp_spline(sample_size_sorted, model_data_sorted, k=3)  # k=3 gives a cubic spline
            y_smooth = spl(xnew)

            # Plot the data for the current model with the Seaborn color palette
            ax.plot(sample_size_sorted, model_data_sorted, label=f'{model.capitalize()} Model', color=palette[i])


        # Set labels and title for each subplot
        ax.set_xlabel('Sample Size', fontsize=12)
        #ax.set_xlabel('Number of Doubles', fontsize=10)
        plt.xticks(fontsize=12)  
        plt.yticks(fontsize=12)
        ax.set_ylabel(metric.capitalize(), fontsize=12)
        ax.set_title(f'{metric.capitalize()} vs Sample Size', fontsize=14)
    # ax.set_title(f'{metric.capitalize()} vs Number of Doubles', fontsize=12)
        ax.axvline(x=1000, color='grey', linestyle='--')
        ax.text(1000+300, ax.get_ylim()[1], '1000', color='grey', verticalalignment='top', rotation=90, fontsize=10)
        ax.axvline(x=2500, color='grey', linestyle='--')
        ax.text(2500+300, ax.get_ylim()[1], '2500', color='grey', verticalalignment='top', rotation=90, fontsize=10)
        ax.axvline(x=5000, color='grey', linestyle='--')
        ax.text(5000+300, ax.get_ylim()[1], '5000', color='grey', verticalalignment='top', rotation=90, fontsize=10)

    # Adjust the layout and add a legend
    plt.tight_layout()
    plt.legend(models)
    #fig.legend(models, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=len(models))
    plt.savefig(f"{out_dir}/all_sampled_mutants_random.png", dpi=300)
    # Show the plot
    plt.show()

def double_triple_heatmap(df, model, metric, out_dir, max_double=1500, max_triple=1500):
    # Filter the DataFrame to include only rows where the gpmap_type is 'blackbox'
    new_data = df[(df['gpmap_type'] == model) & (df['triple_sample_size'] < max_triple) & (df['double_sample_size'] < max_double)]

    # Create a pivot table for the heatmap with 'double_sample_size' on the x-axis and 'triple_sample_size' on the y-axis
    # Use the 'spearman' metric as the values for the heatmap
    spearman_pivot = new_data.pivot(index='triple_sample_size', columns='double_sample_size', values=metric)

    # Plotting the heatmap
    plt.figure(figsize=(16, 16))  # Adjust the figure size for better readability
    sns.heatmap(spearman_pivot, annot=True, cmap="coolwarm", cbar_kws=None, 
                xticklabels=True, yticklabels=True, cbar=False, annot_kws={"size": 13.5})

    # Reverse the y-axis so that the values increase from bottom to top
    plt.gca().invert_yaxis()

    # Set plot labels and title
    plt.title(f'{model.capitalize()} {metric.capitalize()} Score for Different Double and Triple Sample Sizes', fontsize=20)
    plt.xlabel('Number of Doubles', fontsize=16)
    plt.ylabel('Number of Triples', fontsize=16)

    # Rotate x and y axis ticks and increase their size
    plt.xticks(rotation=45, ha='right', fontsize=16)  # Rotate x-axis ticks 45 degrees
    plt.yticks(rotation=45, ha='right', fontsize=16)  # Rotate y-axis ticks 45 degrees

    plt.tight_layout()
    plt.savefig(f"{out_dir}/double_triple_{metric}_{model}_heatmap.png", dpi=300)

    # Show the plot
    plt.show()

def plot_active_learning(df, out_dir, threshold, selection_index):
    # Calculate the difference in Spearman's Rho and R² between consecutive points
    df['delta_spearman_rho'] = df['spearman_rho'].diff()
    df['delta_r2'] = df['r2'].diff()

    # Identify points of diminishing returns
    diminishing_returns_spearman = df[df['delta_spearman_rho'].abs() < threshold]

    # Plot Spearman's Rho and R² vs. training size
    plt.figure(figsize=(10, 6))

    # Dynamically setting y-ticks based on data range
    y_min = round(min(df['spearman_rho'].min(), df['r2'].min()))
    y_max = round(max(df['spearman_rho'].max(), df['r2'].max()))
    
    # Set y-ticks within the range of y_min and y_max with a reasonable step size (e.g., 0.2)
    plt.yticks(np.arange(y_min-1, y_max + 0.2, 0.2))
    
    plt.plot(df['training_size'], df['spearman_rho'], marker='o', label="Spearman's Rho")
    plt.plot(df['training_size'], df['r2'], marker='x', label="R²")

    # Add labels, title, and legend
    plt.title('Spearman\'s Rho and R² vs Training Size')
    plt.xlabel('Training Size')
    plt.ylabel('Metric Value')

    # Highlight the point where diminishing returns start for Spearman's Rho
    if not diminishing_returns_spearman.empty:
        plt.axvline(x=diminishing_returns_spearman['training_size'].iloc[selection_index], 
                    linestyle='--', color='grey', label='Diminishing Returns (Spearman)')
    
    plt.legend()
    
    # Save and show the plot
    plt.savefig(f"{out_dir}/active_learning.png", dpi=300)
    plt.show()

    return diminishing_returns_spearman



def plot_grouped_bar_chart_per_model(metrics, out_dir,colour_palette):
    # Set up the figure and axes with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)  # 3 columns, shared y-axis
    test_data_types = ['all', 'doubles_only', 'triples_only']
    titles = ['Test Data: Full Dataset', 'Test Data: Doubles Only', 'Test Data: Triples Only']

    # Loop over test_data_types and create a plot for each
    i = 0
    for ax, test_data_type, title in zip(axes, test_data_types, titles):
        sns.barplot(
            data=metrics[metrics['test_data_type'] == test_data_type], 
            x='gpmap_type', 
            y='spearman', 
            hue='train_data_type',
            ax=ax,
            palette=colour_palette,  # Use full palette
            edgecolor='black'
        )
        
        # Set title and labels
        ax.set_title(title)
        ax.set_xlabel('Model Type')
        
        if i == 0:  # Only set the y-axis label for the first plot
            ax.set_ylabel("Spearman's Rho")
        else:
            ax.set_ylabel("")  # Clear y label for other plots
        
        ax.set_ylim(-0.4, 1)  # Adjust y-axis range
        ax.legend(title='Training Data Type', loc='lower right')

        i += 1

    # Set a common title for the whole figure
    plt.suptitle(f"Training on all singles, 1500 doubles, and 1500 triples", fontsize=16)
    
    # Adjust layout and save the plot
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Leave space for suptitle
    plt.savefig(f"{out_dir}/spearman_rho_per_model_comparison.png", dpi=300)
    plt.show()


def plot_amino_acid_distribution(df, out_dir, reference_sequence="LYYYGTSYGVL"):
    # Function to calculate amino acid counts by mutation type (single, double, triple, and overall)
    def get_amino_acid_counts(filtered_df):
        sequence_length = len(reference_sequence)
        amino_acids = 'ACDEFGHIKLMNPQRSTVWY'  # Standard amino acids
        position_counts = {i: {aa: 0 for aa in amino_acids} for i in range(sequence_length)}

        for seq in filtered_df['x']:
            for i, aa in enumerate(seq):
                position_counts[i][aa] += 1

        return pd.DataFrame(position_counts).T

    # Split dataset into single, double, triple mutants, and overall
    df['num_mutations'] = df['x'].apply(lambda seq: sum(1 for i, aa in enumerate(seq) if aa != reference_sequence[i]))

    single_mutants = df[df['num_mutations'] == 1]
    double_mutants = df[df['num_mutations'] == 2]
    triple_mutants = df[df['num_mutations'] == 3]
    
    # Create a grid of four subplots (overall, single, double, triple mutants)
    fig, axes = plt.subplots(2, 2, figsize=(18, 12), sharex=True, sharey=True)

    # List of data and titles for each subplot
    datasets = [
        (df, "Overall Distribution"),
        (single_mutants, "Single Mutants"),
        (double_mutants, "Double Mutants"),
        (triple_mutants, "Triple Mutants")
    ]
    
    # Generate the heatmaps without the side bar (color bar)
    for ax, (mutant_df, title) in zip(axes.flatten(), datasets):
        heatmap_df = get_amino_acid_counts(mutant_df)
        sns.heatmap(heatmap_df.T, cmap="coolwarm", annot=True, cbar=False, fmt="d", ax=ax, annot_kws={"size": 14})  # No color bar
        ax.set_title(title, fontsize = 18)
        ax.set_xticklabels(range(1, len(reference_sequence) + 1), fontsize = 16)
        ax.set_xlabel("Position in Sequence", fontsize = 16)
        ax.set_ylabel("Amino Acid", fontsize = 16)
        ax.tick_params(axis='y', labelsize=16)

    # Adjust layout to create more space and avoid overlap
    plt.tight_layout()
    plt.savefig(f"{out_dir}/amino_acid_distribution_grid.png", dpi=300)
    plt.xticks(fontsize = 16)


def plot_mutation_diversity(df, reference_sequence="LYYYGTSYGVL"):
    # Calculate the diversity of amino acids at each position
    sequence_length = len(reference_sequence)
    amino_acid_sets = {i: set() for i in range(sequence_length)}

    for seq in df['x']:
        for i, aa in enumerate(seq):
            amino_acid_sets[i].add(aa)

    # Count how many unique amino acids are found at each position
    diversity_counts = {i: len(amino_acid_sets[i]) for i in range(sequence_length)}

    # Plot the diversity of mutations at each position
    plt.figure(figsize=(10, 6))
    plt.bar(diversity_counts.keys(), diversity_counts.values(), color='green')
    plt.xlabel("Position in Sequence")
    plt.ylabel("Number of Unique Amino Acids")
    plt.title("Diversity of Mutations at Each Position in 11-mer Sequences")
    plt.xticks(range(0, sequence_length))
    plt.show()
    
def plot_mutation_trends_over_iterations(df, reference_sequence="LYYYGTSYGVL"):
    # Count mutations at each position per iteration
    sequence_length = len(reference_sequence)
    df['iteration'] = df['iteration']  # Assuming 'iteration' is a column in the DataFrame

    iteration_mutations = {i: [] for i in range(sequence_length)}

    for _, row in df.iterrows():
        seq = row['x']
        iteration = row['iteration']
        for i, aa in enumerate(seq):
            if aa != reference_sequence[i]:
                iteration_mutations[i].append(iteration)

    # Count how often each position is mutated per iteration
    for position, iterations in iteration_mutations.items():
        plt.hist(iterations, bins=range(1, df['iteration'].max() + 2), alpha=0.5, label=f"Position {position}")

    plt.xlabel("Iteration")
    plt.ylabel("Mutation Count")
    plt.title("Mutation Trends Over Iterations")
    plt.legend()
    plt.show()

def plot_spearman_comparison(results_all_epitopes, results_single_epitope, 
                            gpmap_type='pairwise', max_sample_size=2000, save_path=None):
    """
    Plots a comparison of Spearman's Rho for pairwise models with and without epitope switching.

    Parameters:
    - results_all_epitopes (DataFrame): Data for all epitopes.
    - results_single_epitope (DataFrame): Data filtered for single epitopes.
    - gpmap_type (str): Type of gpmap model to filter by (default: 'pairwise').
    - max_sample_size (int): Maximum sample size threshold for filtering (default: 2000).
    - save_path (str): Path to save the plot image (optional).

    Returns:
    - None: Displays and optionally saves the plot.
    """

    # Filter results for 'all epitopes'
    filtered_all_epitopes = results_all_epitopes[
        (results_all_epitopes['double_sample_size'] == results_all_epitopes['triple_sample_size']) & 
        (results_all_epitopes['gpmap_type'] == gpmap_type) & 
        (results_all_epitopes['triple_sample_size'] < max_sample_size)
    ]

    # Filter results for 'single epitope'
    filtered_single_epitope = results_single_epitope[
        (results_single_epitope['double_sample_size'] == results_single_epitope['triple_sample_size']) & 
        (results_single_epitope['gpmap_type'] == gpmap_type)
    ]

    # Sort data for consistent plotting
    filtered_all_epitopes = filtered_all_epitopes.sort_values(by='double_sample_size')
    filtered_single_epitope = filtered_single_epitope.sort_values(by='double_sample_size')

    # Plotting
    plt.figure(figsize=(10, 8))

    # Line for 'all epitopes'
    plt.plot(
        filtered_all_epitopes['double_sample_size'] + filtered_all_epitopes['triple_sample_size'], 
        filtered_all_epitopes['spearman'], 
        label="All Epitopes", marker='o'
    )

    # Line for 'single epitope'
    plt.plot(
        filtered_single_epitope['double_sample_size'] + filtered_single_epitope['triple_sample_size'], 
        filtered_single_epitope['spearman'], 
        label="Single Epitope", marker='o'
    )

    # Adding labels, title, and legend
    plt.title("Comparison of Spearman's Rho for Pairwise Models with and without Epitope Switching")
    plt.xlabel("Sample Size")
    plt.ylabel("Spearman's Rho")
    plt.legend()

    # Save the plot if a save path is provided
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved to {save_path}")

    # Display the plot
    plt.show()

def plot_r2_vs_sample_size(df, title="R2 vs Sample Size"):
    """
    Plots 'r2' vs 'sample_size' from a given DataFrame with an elegant design.
    
    Parameters:
    - df: DataFrame with columns 'sample_size' and 'r2'.
    - title: Title for the plot.
    """
    # Sort data for consistent plotting
    df = df.sort_values(by='sample_size')

    # Set a Seaborn theme
    sns.set_theme(style="whitegrid", font_scale=1.2)

    # Create the figure and axes
    plt.figure(figsize=(12, 4))
    
    # Smooth lineplot with error shading (if confidence intervals exist)
    sns.lineplot(
        data=df,
        x='sample_size',
        y='r2',
        marker='o',
        color='#007ACC',
        linewidth=2.5,
        markersize=10,
        alpha=0.8,
    )

    # Customize axis labels and title
    plt.xlabel('Sample Size', fontsize=18, labelpad=15)
    plt.ylabel('R2', fontsize=18, labelpad=15)
    plt.title(title, fontsize=18, pad=20, weight='bold')
    plt.axhline(y=0.5, linestyle='--', linewidth=2)
    plt.axvline(x=375, linestyle='--', linewidth=2)
    


    # Customize ticks
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    # Remove unnecessary spines for a clean look
    sns.despine()

    # Add gridlines for readability
    plt.grid(visible=True, linestyle='--', alpha=0.7)

    # Save and show
    plt.tight_layout()
    plt.savefig("/Users/uw20204/Desktop/1fbi_x_r2_vs_sample_size_plot.png", dpi=300)
    plt.show()

def plot_heatmap_single_antigen(results_list, datasets, out_dir, gpmap):
    """
    Plots a heatmap where:
    - X-axis: Sample sizes
    - Y-axis: Dataset names
    - Color: Spearman's rho for the blackbox model
    """
    # Initialize a dictionary to store R2 for each dataset and sample size
    heatmap_data = {}

    # Iterate over the single result
    for idx, result in enumerate(results_list):
        dataset_name = datasets[idx]

        # Filter the data for the 'gpmap' model
        blackbox_data = result.loc[(result['gpmap_type'] == gpmap)]

        # Group by sample size and calculate the mean R2
        spearman_rho = blackbox_data.groupby('sample_size')['r2'].mean()
        heatmap_data[dataset_name] = spearman_rho

    # Convert the dictionary to a DataFrame for the heatmap
    heatmap_df = pd.DataFrame(heatmap_data)

    # Plot the heatmap
    plt.figure(figsize=(18, 2))
    ax = sns.heatmap(
        heatmap_df.T, annot=True, cmap='coolwarm', fmt=".2f", cbar_kws={'label': "R2"}
    )
    for text in ax.texts:  # Iterate over the text elements
        text.set_rotation(90)  # Set the rotation
        text.set_size(10)
    # Customize the color bar
    cbar = ax.collections[0].colorbar
    cbar.set_label("R2", fontsize=14)
    cbar.ax.tick_params(labelsize=12)

    # Customize labels
    plt.xlabel('Total number of doubles & triples', fontsize=14)
    plt.ylabel('Antigen', fontsize=14)
    plt.title(f'R2 for {gpmap} Model', fontsize=16)

    # Save the plot
    plt.tight_layout()
    plt.savefig(f"{out_dir}/r2_heatmap_{gpmap}.png", dpi=300)
    plt.show()

def plot_heatmap_of_spearman(results_list, datasets, out_dir, gpmap):
    """
    Plots a heatmap where:
    - X-axis: Sample sizes
    - Y-axis: Dataset names (sorted by numeric prefix, then alphabetically)
    - Color: Spearman's rho for the blackbox model
    """
    # Initialize a dictionary to store Spearman's rho for each dataset and sample size
    heatmap_data = {}

    # Iterate over each result (dataset)
    for idx, result in enumerate(results_list):
        dataset_name = datasets[idx]

        # Filter the data for the 'blackbox' model and the 'spearman' metric
        blackbox_data = result.loc[(result['gpmap_type'] == gpmap)]

        # Get unique sample sizes
        sample_sizes = blackbox_data['sample_size'].unique()

        # Store Spearman's rho values in the dictionary with the dataset name
        spearman_rho = blackbox_data.groupby('sample_size')['r2'].mean()
        heatmap_data[dataset_name] = spearman_rho

    # Convert the dictionary to a DataFrame for the heatmap
    heatmap_df = pd.DataFrame(heatmap_data)

    # Extract numeric prefix and the remaining part, then sort accordingly
    sorted_datasets = sorted(
        heatmap_df.columns,
        key=lambda x: (int(re.match(r'\d+', x).group()), x)  # First by number, then alphabetically
    )
    heatmap_df = heatmap_df[sorted_datasets]

    # Plot the heatmap
    plt.figure(figsize=(20, 20))
    ax = sns.heatmap(heatmap_df.T, annot=False, cmap = 'coolwarm', fmt=".2f", cbar_kws={'label': "R2"})
    
    # Increase the font size of the color bar label
    cbar = ax.collections[0].colorbar
    cbar.set_label("R2", fontsize=16)  # Set label and font size
    cbar.ax.tick_params(labelsize=14)  # Set font size for color bar ticks

    plt.xticks(ha='right', fontsize=14)  # Rotate x-axis ticks 45 degrees
    plt.yticks(ticks=np.arange(len(heatmap_df.columns)) + 0.5, labels=heatmap_df.columns, fontsize=14)

    # Set labels and title
    plt.xlabel('Total number of doubles & triples', fontsize=16)
    plt.ylabel('Antigens', fontsize=16)
    plt.title(f'R2 for {gpmap} Model Across Antigens Trained on All 209 Single Mutants and Different Numbers of Doubles & Triples', fontsize=20)

    # Save the plot
    plt.tight_layout()
    plt.savefig(f"{out_dir}/r2_heatmap_{gpmap}.png", dpi=300)
    plt.show()



def plot_num_epitopes_vs_r2(num_epitopes_r2_df, output_path=None):
    """
    Plots the relationship between the number of epitopes and R² values.

    Parameters:
    - num_epitopes_r2_df (DataFrame): Data containing 'num_epitopes' and 'r2' columns.
    - output_path (str): Path to save the plot image (optional).

    Returns:
    - None: Displays or saves the plot.
    """
    sns.set(style="whitegrid", context="talk")  # Seaborn style
    plt.figure(figsize=(10, 6))

    # Scatter plot
    plt.scatter(num_epitopes_r2_df['num_epitopes'], num_epitopes_r2_df['r2'], 
                color='royalblue', edgecolor='black', alpha=0.7, s=100, label="Antigens")

    # Add regression line
    sns.regplot(x='num_epitopes', y='r2', data=num_epitopes_r2_df, scatter=False, color='black', ci=None)

    # Add labels and title
    plt.title("Relationship Between Number of Epitopes and R²", fontsize=18, pad=15)
    plt.xlabel("Number of Epitopes", fontsize=14, labelpad=10)
    plt.ylabel("R² Value", fontsize=14, labelpad=10)
    plt.legend(loc='best', fontsize=12)

    # Grid lines and customisation
    plt.grid(visible=True, linestyle="--", alpha=0.5)
    plt.tight_layout()

    # Save or display the plot
    if output_path:
        plt.savefig(output_path, dpi=300)
        print(f"Plot saved to {output_path}")
    plt.show()


def plot_min_doubles_triples(results_list, datasets, gpmap, threshold=0.6):
    """
    Plots the minimum number of doubles/triples required to reach >0.6 R2 for each antigen.
    - X-axis: Antigens (sorted by number of doubles/triples required)
    - Y-axis: Minimum number of doubles/triples required
    """
    min_doubles_triples = {}

    # Iterate over each result (dataset)
    for idx, result in enumerate(results_list):
        dataset_name = datasets[idx]

        # Filter the data for the specified gpmap and calculate the minimum sample size to achieve the threshold R2
        filtered_data = result[result['gpmap_type'] == gpmap]
        above_threshold = filtered_data[filtered_data['r2'] > threshold]
        
        if not above_threshold.empty:
            min_sample_size = above_threshold['sample_size'].min()
        else:
            min_sample_size = float('inf')  # If no R2 exceeds the threshold, set to infinity
        
        min_doubles_triples[dataset_name] = min_sample_size

    # Convert to DataFrame for plotting
    min_doubles_triples_df = pd.DataFrame.from_dict(min_doubles_triples, orient='index', columns=['MinDoublesTriples'])
    min_doubles_triples_df = min_doubles_triples_df.reset_index()
    min_doubles_triples_df.columns = ['Antigen', 'MinDoublesTriples']

    # Sort by the minimum number of doubles/triples
    min_doubles_triples_df = min_doubles_triples_df.sort_values(by='MinDoublesTriples')

    # Plot the bar chart
    plt.figure(figsize=(14, 8))
    plt.bar(min_doubles_triples_df['Antigen'], min_doubles_triples_df['MinDoublesTriples'], color='skyblue')
    plt.xticks(rotation=90, ha='right', fontsize=10)
    plt.xlabel('Antigens', fontsize=12)
    plt.ylabel('Total Number of Doubles/Triples', fontsize=12)
    plt.yscale('log')
    plt.title(f'Minimum Doubles/Triples for R2 > {threshold} ({gpmap})', fontsize=14)
    plt.tight_layout()

    # Save the plot
    plt.savefig(f"min_doubles_triples_{gpmap}.png", dpi=300)
    plt.show()

    return min_doubles_triples_df

import matplotlib.pyplot as plt
import seaborn as sns

def plot_max_spearman_additive(df, output_path=None, palette="viridis"):
    """
    Plots the Max Spearman's Rho for additive models sorted by antigen prefix.

    Parameters:
    - df (DataFrame): A DataFrame containing 'antigen' and 'max_spearman_additive' columns.
    - output_path (str): Path to save the figure (optional). Defaults to None.
    - palette (str): Colour palette for the barplot.

    Returns:
    - None: Displays and optionally saves the plot.
    """

    # Helper function to extract numeric prefixes for sorting
    def extract_prefix(antigen):
        """Extract numeric prefix for sorting."""
        try:
            return int(antigen.split('_')[0])  # Assume prefix is numeric before '_'
        except ValueError:
            return float('inf')  # Handle non-numeric prefixes, placing them last

    # Add a column for sorting
    df['prefix'] = df['antigen'].apply(extract_prefix)

    # Sort the DataFrame by the new prefix column
    df_sorted = df.sort_values(by=['prefix', 'antigen'])

    # Plotting
    plt.figure(figsize=(14, 6))
    sns.barplot(
        x='antigen',
        y='max_spearman_additive',
        data=df_sorted,
        order=df_sorted['antigen'],  # Explicit x-axis order
        palette=palette
    )

    # Customise labels and title
    plt.xlabel("Antigen", fontsize=14)
    plt.ylabel("Max Spearman's Rho", fontsize=14)
    plt.title("Max Spearman's Rho for Additive Models for Each Antigen", fontsize=16)

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=90, fontsize=10)

    # Adjust layout
    plt.tight_layout()

    # Save the plot if output_path is specified
    if output_path:
        plt.savefig(output_path, dpi=300)
        print(f"Plot saved to {output_path}")

    # Show the plot
    plt.show()

colour_palette = sns.color_palette("pastel", 4)

def generate_panel_of_plots(data_dir, output_path, gpmap, max_plots_per_page=30, colour_palette = colour_palette):
    print(data_dir)
    # Debug: List all files in the data directory
    all_files = os.listdir(data_dir)
    print("All files in directory:", all_files)

    # Filter prediction, grid, and test data files
    prediction_files = [f for f in all_files if "prediction_data" in f and 'nonlinear' in f]
    prediction_files = [f for f in prediction_files if gpmap in f]
    grid_files = [f for f in all_files if "grid_data" in f and 'nonlinear' in f]
    grid_files = [f for f in grid_files if gpmap in f]
    test_files = [f for f in all_files if "test_data" in f and 'linear' in f]
    test_files = [f for f in test_files if gpmap in f]

    print("Prediction files found:", prediction_files)
    print("Grid files found:", grid_files)
    print("Test files found:", test_files)

    # Extract and sort prefixes
    def safe_sort_key(prefix):
        """Extract the numeric part of the prefix for sorting. Non-numeric prefixes are sorted after numeric ones."""
        try:
            # Extract the first part and convert to an integer
            return int(prefix.split('_')[0]), prefix
        except ValueError:
            # Non-numeric prefixes will be sorted as strings
            return float('inf'), prefix

    prefixes = sorted(
        {f.split('_')[0] + "_" + f.split('_')[1] for f in prediction_files},
        key=safe_sort_key
    )
    print("Sorted prefixes found:", prefixes)

    # Initialise lists for plot data
    prediction_data = []
    ge_data = []

    # Process each matching set of files
    for prefix in prefixes:
        # Match files with the same prefix
        pred_file = next((f for f in prediction_files if f.startswith(prefix)), None)
        grid_file = next((f for f in grid_files if f.startswith(prefix)), None)
        test_file = next((f for f in test_files if f.startswith(prefix)), None)

        # Ensure all files are found
        if not all([pred_file, grid_file, test_file]):
            print(f"Warning: Missing files for prefix {prefix}. Skipping.")
            continue

        # Load prediction, grid, and test data
        pred_data = pd.read_csv(os.path.join(data_dir, pred_file))
        grid_data = pd.read_csv(os.path.join(data_dir, grid_file))
        test_data = pd.read_csv(os.path.join(data_dir, test_file))

        # Extract data
        phi_test = test_data['phi_test'].values
        y_test = pred_data['y_test'].values

        # Check and filter mismatched sizes
        if len(phi_test) != len(y_test):
            print(f"Warning: Mismatched sizes for prefix {prefix}. Adjusting lengths.")
            min_length = min(len(phi_test), len(y_test))
            phi_test = phi_test[:min_length]
            y_test = y_test[:min_length]

        # Append data for later plotting
        prediction_data.append((prefix, pred_data))
        ge_data.append((prefix, grid_data, phi_test, y_test))

    # Split data into batches
    def chunk_data(data, chunk_size):
        """Split data into chunks of max chunk_size."""
        for i in range(0, len(data), chunk_size):
            yield data[i:i + chunk_size]

    # Create and save multiple panels for prediction plots
    for batch_idx, batch in enumerate(chunk_data(prediction_data, max_plots_per_page), start=1):
        num_plots = len(batch)
        num_cols = 5
        num_rows = math.ceil(num_plots / num_cols)
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(5 * num_cols, 5 * num_rows))
        axes = axes.flatten()

        for i, (prefix, pred_data) in enumerate(batch):
            ax = axes[i]
            # Scatter plot

            ax.scatter(pred_data['yhat_test'], pred_data['y_test'], alpha=0.9, s=5, color=colour_palette[0])
            # Diagonal line
            ax.plot([min(pred_data['yhat_test']), max(pred_data['yhat_test'])],
                    [min(pred_data['yhat_test']), max(pred_data['yhat_test'])], '--k', color = 'black', linewidth = 3)
            # Calculate R^2
            yhat = pred_data['yhat_test']
            y = pred_data['y_test']
            r2 = np.corrcoef(yhat, y)[0, 1] ** 2
            # Add annotation for R^2
            ax.text(0.05, 0.95, f"$R^2$ = {r2:.2f}", fontsize=14, transform=ax.transAxes, 
                    verticalalignment='top', horizontalalignment='left')
            # Set title and labels
            ax.set_title(f'{prefix}', fontsize=22)
            ax.set_xlabel('$\hat{y}$', fontsize=18)
            ax.set_ylabel('$y$', fontsize=18)
            ax.tick_params(axis='both', labelsize=18)

        # Remove unused axes
        for i in range(num_plots, len(axes)):
            fig.delaxes(axes[i])

        plt.tight_layout()
        plt.savefig(os.path.join(output_path, f"prediction_panell_{gpmap}_{batch_idx}_latest.png"))
        plt.close()

    # Create and save multiple panels for GE measurement process plots
    for batch_idx, batch in enumerate(chunk_data(ge_data, max_plots_per_page), start=1):
        num_plots = len(batch)
        num_cols = 5
        num_rows = math.ceil(num_plots / num_cols)
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(5 * num_cols, 5 * num_rows))
        axes = axes.flatten()

        for i, (prefix, grid_data, phi_test, y_test) in enumerate(batch):
            ax = axes[i]
            ax.fill_between(grid_data['phi_grid'], grid_data['yqs_lower'], grid_data['yqs_upper'], alpha=0.2, color=colour_palette[1])
            ax.plot(grid_data['phi_grid'], grid_data['yhat_grid'], color='black', linewidth = 3)
            ax.scatter(phi_test, y_test, alpha=0.9, s=5, color=colour_palette[2])
            ax.set_title(f'{prefix}', fontsize=22)
            ax.set_xlabel('$\phi$', fontsize=18)
            ax.set_ylabel('$y$', fontsize=18)
            ax.tick_params(axis='both', labelsize=18)

        # Remove unused axes
        for i in range(num_plots, len(axes)):
            fig.delaxes(axes[i])

        plt.tight_layout()
        plt.savefig(os.path.join(output_path, f"ge_process_panel_{gpmap}_{batch_idx}_latest.png"), dpi=300)
        plt.close()


def generate_comparison_panels(data_dir, output_path, gpmap_types, max_plots=5):
    # Debug: List all files in the data directory
    all_files = os.listdir(data_dir)

    # Filter files for each GP map type
    gpmap_files = {gpmap: [f for f in all_files if gpmap in f] for gpmap in gpmap_types}

    # Extract and sort prefixes
    def safe_sort_key(prefix):
        """Extract the numeric part of the prefix for sorting. Non-numeric prefixes are sorted after numeric ones."""
        try:
            # Extract the first part and convert to an integer
            return int(prefix.split('_')[0]), prefix
        except ValueError:
            # Non-numeric prefixes will be sorted as strings
            return float('inf'), prefix

    sorted_prefixes = {
        gpmap: sorted(
            {f.split('_')[0] + "_" + f.split('_')[1] for f in files if "prediction_data" in f},
            key=safe_sort_key
        )[:max_plots]
        for gpmap, files in gpmap_files.items()
    }

    # Prepare data for prediction and measurement process plots
    prediction_data = {}
    ge_data = {}
    for gpmap, prefixes in sorted_prefixes.items():
        prediction_data[gpmap] = []
        ge_data[gpmap] = []
        for prefix in prefixes:
            # Find corresponding files
            pred_file = next((f for f in gpmap_files[gpmap] if f.startswith(prefix) and "prediction_data" in f), None)
            grid_file = next((f for f in gpmap_files[gpmap] if f.startswith(prefix) and "grid_data" in f), None)
            test_file = next((f for f in gpmap_files[gpmap] if f.startswith(prefix) and "test_data" in f), None)
            if pred_file:
                pred_data = pd.read_csv(os.path.join(data_dir, pred_file))
                prediction_data[gpmap].append((prefix, pred_data))
            if grid_file and test_file:
                grid_data = pd.read_csv(os.path.join(data_dir, grid_file))
                test_data = pd.read_csv(os.path.join(data_dir, test_file))
                phi_test = test_data['phi_test'].values
                y_test = test_data['y_test'].values
                ge_data[gpmap].append((prefix, grid_data, phi_test, y_test))

    # Generate prediction plot panel
    def generate_panel(data, file_name, x_label, y_label, title_key):
        num_gpmap_types = len(gpmap_types)
        # Reduce title row space using smaller height ratio
        fig = plt.figure(figsize=(5 * max_plots, 5 * num_gpmap_types))
        grid = GridSpec(num_gpmap_types * 2, max_plots, figure=fig, height_ratios=[0.2, 1] * num_gpmap_types)

        for row_idx, gpmap in enumerate(gpmap_types):
            # Add a title spanning all columns for this GP map type
            title_ax = fig.add_subplot(grid[row_idx * 2, :])
            title_ax.axis("off")
            title_ax.text(0.5, 0.5, gpmap.upper(), fontsize=18, weight='bold', ha='center', va='center')

            for col_idx in range(max_plots):
                plot_idx = row_idx * 2 + 1
                ax = fig.add_subplot(grid[plot_idx, col_idx])
                if col_idx < len(data[gpmap]):
                    entry = data[gpmap][col_idx]
                    prefix = entry[0]
                    if title_key == "prediction":
                        pred_data = entry[1]
                        ax.scatter(pred_data['yhat_test'], pred_data['y_test'], alpha=0.3)
                        ax.plot([min(pred_data['yhat_test']), max(pred_data['yhat_test'])],
                                [min(pred_data['yhat_test']), max(pred_data['yhat_test'])], '--k')
                    elif title_key == "measurement":
                        grid_data, phi_test, y_test = entry[1], entry[2], entry[3]
                        ax.fill_between(grid_data['phi_grid'], grid_data['yqs_lower'], grid_data['yqs_upper'], alpha=0.2, color='C1')
                        ax.plot(grid_data['phi_grid'], grid_data['yhat_grid'], color='C1')
                        ax.scatter(phi_test, y_test, alpha=0.3, color='C0', s=10)
                    ax.set_title(f'{prefix}', fontsize=12)
                    ax.set_xlabel(x_label, fontsize=10)
                    ax.set_ylabel(y_label, fontsize=10)
                    ax.tick_params(axis='both', labelsize=8)
                else:
                    ax.axis('off')  # Hide unused axes

        # Adjust layout and save the figure
        plt.tight_layout(h_pad=1.5)
        plt.savefig(os.path.join(output_path, file_name), dpi=300)
        plt.close()

    # Create prediction plot panel
    generate_panel(
        prediction_data,
        file_name="prediction_comparison_panel.png",
        x_label="$\hat{y}$",
        y_label="$y$",
        title_key="prediction"
    )

    # Create GE measurement process plot panel
    generate_panel(
        ge_data,
        file_name="ge_process_comparison_panel.png",
        x_label="$\phi$",
        y_label="$y$",
        title_key="measurement"
    )