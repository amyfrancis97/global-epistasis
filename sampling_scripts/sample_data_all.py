import os
import argparse
import pandas as pd
import gc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from utils import *  # Assuming necessary functions like fit_model, calculate_metrics, etc.
#from config import *  # Assuming necessary configuration values

# Random sampling implementation
def random_sample_experiment(name, trainval_file, test_file, num_samples_list, gpmap_types, output_path):
    # Read in trainval_df and test_df
    trainval_df = pd.read_csv(trainval_file, sep="\t")
    test_df = pd.read_csv(test_file, sep="\t")
    out_prefix = trainval_file.split("_trainval_df")[0]
    out_prefix = out_prefix.split("data/")[1]
    results = []

    # Loop over different random sample sizes
    for num_samples in num_samples_list:
        print(f"Processing num_samples={num_samples}")

        # Filter rows where dist == 1
        single_rows = trainval_df[trainval_df['dist'] == 1]

        # Sample rows where dist != 1
        sampled_rows = trainval_df[trainval_df['dist'] != 1].sample(n=num_samples, random_state=42)

        # Concatenate the two DataFrames
        sampled_trainval_df = pd.concat([single_rows, sampled_rows]).reset_index(drop=True)


        # Dictionary to store metrics for different G-P map types
        metrics_for_this_run = []

        # Train and evaluate each type of G-P map
        for gpmap_type in gpmap_types:
            print(f"Training model for gpmap_type={gpmap_type}")

            # Fit the model
            model = fit_model(sampled_trainval_df, gpmap_type)

            # Compute metrics on the test dataset
            mse, r2, spearman = calculate_metrics(model, test_df)
            metrics_for_this_run.append({
                'gpmap_type': gpmap_type,
                'sample_size': num_samples,
                'mse': mse,
                'r2': r2,
                'spearman': spearman
            })

        # Append results for this random sample size
        for metric in metrics_for_this_run:
            print(metric)
            results.append(metric)

        # Clean up memory by deleting models and collecting garbage
        del model
        gc.collect()

    # Convert results to a DataFrame and return
    results_df = pd.DataFrame(results)
    
    # Save results to a file
    results_df.to_csv(os.path.join(output_path, f"{out_prefix}_{name}_random_downsampling_results.tsv"), sep="\t", index=False)

    return results_df

# Command-line argument parsing
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Random Sampling Experiment')

    # Add arguments for the trainval_df, test_df, and output_path
    parser.add_argument('--name', type=str, required=True, help='Prefix name of files')
    parser.add_argument('--trainval_file', type=str, required=True, help='Path to the trainval_df file')
    parser.add_argument('--test_file', type=str, required=True, help='Path to the test_df file')
    parser.add_argument('--output_path', type=str, required=True, help='Directory to save output files')
    parser.add_argument('--num_samples_list', type=str, default="50-3000-25,3500-45000-250", help='Range of sample sizes (e.g., "50-3000-25,3500-45000-250")')
    parser.add_argument('--gpmap_types', type=str, default="additive,neighbor,pairwise,blackbox", help='Comma-separated list of GP map types')

    args = parser.parse_args()

    # Parse num_samples_list into a list of integers
    num_samples_list = []
    for range_str in args.num_samples_list.split(','):
        start, end, step = map(int, range_str.split('-'))
        num_samples_list.extend(range(start, end, step))

    # Parse gpmap_types into a list
    gpmap_types = args.gpmap_types.split(',')

    # Run the experiment
    random_results = random_sample_experiment(
        name=args.name,
        trainval_file=args.trainval_file,
        test_file=args.test_file,
        num_samples_list=num_samples_list,
        gpmap_types=gpmap_types,
        output_path=args.output_path
    )

    # Plot the results
    plot_sample_size_metrics(random_results, args.output_path)

