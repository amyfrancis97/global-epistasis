# This script samples all singles, and different combinations of specified double and triple mutants.

import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gc
from utils import *
#from config import *

def run_sampling_experiment(name, trainval_file, test_file, output_path, gpmap_types, range_to_sample):

    # Read in trainval_df and test_df from the specified directory
    trainval_df = pd.read_csv(trainval_file, sep="\t")
    test_df = pd.read_csv(test_file, sep="\t")

    print(trainval_file.split("_trainval_df"))
    out_prefix = trainval_file.split("_trainval_df")[0]
    out_prefix = out_prefix.split("data/")[1]
    print(out_prefix)

    # Precompute single training data and validation data once to avoid repetitive computation
    single_training_data = trainval_df[trainval_df['dist'] == 1]

    # Available counts of double and triple mutants
    available_doubles = trainval_df[trainval_df['dist'] == 2].shape[0]
    available_triples = trainval_df[trainval_df['dist'] == 3].shape[0]

    # Initialise results storage
    results_2 = []

    # Main loop over sample sizes
    for num_doubles in range_to_sample:
        if num_doubles > available_doubles:
            print(f"Skipping num_doubles={num_doubles} (only {available_doubles} available)")
            continue  # Skip if the sample size exceeds available doubles

        for num_triples in range_to_sample:
            if num_triples > available_triples:
                print(f"Skipping num_triples={num_triples} (only {available_triples} available)")
                continue  # Skip if the sample size exceeds available triples

            print(f"Processing num_doubles={num_doubles}, num_triples={num_triples}")

            # Create a new dataset by concatenating sampled data
            new_double_training_data = trainval_df[trainval_df['dist'] == 2].sample(num_doubles)
            new_triple_training_data = trainval_df[trainval_df['dist'] == 3].sample(num_triples)
            new_dataset = pd.concat([single_training_data, new_double_training_data, new_triple_training_data]).reset_index(drop=True)

            # Split the new dataset into training and validation (5% validation)
            train_data, val_data = train_test_split(new_dataset, test_size=0.05, random_state=42)

            # Mark the validation data
            train_data['validation'] = False
            val_data['validation'] = True

            # Concatenate the training and validation data back together
            new_dataset = pd.concat([train_data, val_data])

            # Train and evaluate each type of G-P map
            metrics_for_this_run = []
            for gpmap_type in gpmap_types:
                print(f"Training model for gpmap_type={gpmap_type}")

                # Fit the model by passing new_dataset as trainval_df
                model = fit_model(new_dataset, gpmap_type)

                # Compute metrics on the test dataset
                mse, r2, spearman = calculate_metrics(model, test_df)
                metrics_for_this_run.append({
                    'gpmap_type': gpmap_type,
                    'mse': mse,
                    'r2': r2,
                    'spearman': spearman
                })

            # Append results for this combination of doubles and triples
            for metric in metrics_for_this_run:
                metric.update({
                    'double_sample_size': num_doubles,
                    'triple_sample_size': num_triples
                })
                print(metric)
                results_2.append(metric)

            # Clean up memory by deleting models and collecting garbage
            del model
            gc.collect()

    # Convert results to a DataFrame and save
    results_2_df = pd.DataFrame(results_2)
    results_2_df.to_csv(os.path.join(output_path, f"{out_prefix}_{name}_sampling_results_double_triple_final.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    # Command-line argument parsing
    parser = argparse.ArgumentParser(description='Double and Triple Sampling Experiment')

    # Add arguments for trainval_df, test_df, output_path, name_prefix, and other parameters
    parser.add_argument('--name', type=str, required=True, help='Prefix for output file names')
    parser.add_argument('--trainval_file', type=str, required=True, help='Path to the trainval_df file')
    parser.add_argument('--test_file', type=str, required=True, help='Path to the test_df file')
    parser.add_argument('--output_path', type=str, required=True, help='Directory to save output files')
    parser.add_argument('--gpmap_types', type=str, default="additive,neighbour,pairwise,blackbox", help='Comma-separated list of GP map types')
    parser.add_argument('--range_to_sample', type=str, default="1,50-2000-50,4000-40000-2000", help='Comma-separated ranges for double and triple sample sizes')

    # Parse arguments
    args = parser.parse_args()

    # Parse range_to_sample into a list of integers
    range_to_sample = []
    for range_str in args.range_to_sample.split(','):
        if '-' in range_str:
            start, end, step = map(int, range_str.split('-'))
            range_to_sample.extend(range(start, end, step))
        else:
            range_to_sample.append(int(range_str))

    # Parse gpmap_types into a list
    gpmap_types = args.gpmap_types.split(',')

    # Run the experiment
    run_sampling_experiment(
        name=args.name,
        trainval_file=args.trainval_file,
        test_file=args.test_file,
        output_path=args.output_path,
        gpmap_types=gpmap_types,
        range_to_sample=range_to_sample
    )

