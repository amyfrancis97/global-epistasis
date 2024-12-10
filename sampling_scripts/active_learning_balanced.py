import os
import argparse
import pandas as pd
from utils import *

def active_learning_experiment(name, gpmap_type, trainval_file, test_file, output_path, n_queries=10, query_size=50, n_initial_doubles=10, n_initial_triples=5, performance_threshold=None):

    # Read in trainval_df and test_df from the specified directory
    trainval_df = pd.read_csv(trainval_file, sep="\t")
    test_df = pd.read_csv(test_file, sep="\t")

    # Start with an initial training set that includes all 'dist == 1' samples (single mutants)
    single_mutants = trainval_df[trainval_df['dist'] == 1]

    # Sample double and triple mutants
    double_mutants = trainval_df[trainval_df['dist'] == 2].sample(n=n_initial_doubles, random_state=42)
    triple_mutants = trainval_df[trainval_df['dist'] == 3].sample(n=n_initial_triples, random_state=42)

    # Combine the initial training set
    trainval_df_filtered = pd.concat([single_mutants, double_mutants, triple_mutants]).reset_index(drop=True)

    # Define pool: only training data where 'dist != 1' (remaining unlabeled training data)
    remaining_pool = trainval_df[trainval_df['dist'] != 1]
    remaining_pool = remaining_pool[~remaining_pool['x'].isin(trainval_df_filtered['x'])]

    # Downsample 'dist == 3' to balance with 'dist == 2'
    dist_2_pool = remaining_pool[remaining_pool['dist'] == 2]
    dist_3_pool = remaining_pool[remaining_pool['dist'] == 3].sample(n=len(dist_2_pool), random_state=42)

    # Combine the downsampled pools back together
    remaining_pool_balanced = pd.concat([dist_2_pool, dist_3_pool]).reset_index(drop=True)

    # Store performance metrics over time
    metrics_over_time = []

    # Train initial model on the initial training set
    train_data, val_data = train_test_split(trainval_df_filtered, test_size=0.05, random_state=42)
    train_data['validation'] = False
    val_data['validation'] = True
    trainval_df_filtered = pd.concat([train_data, val_data]).reset_index(drop=True)

    # Train the model on the initial training set
    model = fit_model(trainval_df_filtered, gpmap_type)

    selected_sequences = []
    
    # Active Learning Loop
    for i in range(n_queries):
        print(f"Active Learning Iteration {i+1}/{n_queries}")

        # Check if remaining pool is empty
        if remaining_pool_balanced.empty:
            print("Remaining pool is empty, stopping active learning.")
            break

        # Compute uncertainty for the remaining pool using your preferred metric
        uncertainties = compute_uncertainty(model, remaining_pool_balanced)

        # Select the top 'query_size' uncertain samples
        if len(remaining_pool_balanced) < query_size:
            print(f"Only {len(remaining_pool_balanced)} samples left in the pool, selecting all.")
            query_size = len(remaining_pool_balanced)  # Adjust query size to avoid exceeding pool size

        uncertain_indices = uncertainties.argsort()[-query_size:]
        new_samples = remaining_pool_balanced.iloc[uncertain_indices]

        selected_sequences.append(new_samples['x'].tolist()) 

        # Add the selected uncertain samples to the training set
        trainval_df_filtered = pd.concat([trainval_df_filtered, new_samples]).reset_index(drop=True)

        # Remove the selected samples from the pool
        remaining_pool_balanced = remaining_pool_balanced.drop(new_samples.index)

        # Re-split the dataset into training (95%) and validation (5%)
        train_data, val_data = train_test_split(trainval_df_filtered, test_size=0.05, random_state=42)
        train_data['validation'] = False
        val_data['validation'] = True
        trainval_df_filtered = pd.concat([train_data, val_data]).reset_index(drop=True)

        # Retrain the model with the updated training set
        model = fit_model(trainval_df_filtered, gpmap_type)

        # Evaluate the model on the test set
        mse, r2, spearman_rho = calculate_metrics(model, test_df)

        # Track distribution of singles, doubles, and triples
        num_singles = (trainval_df_filtered['dist'] == 1).sum()
        num_doubles = (trainval_df_filtered['dist'] == 2).sum()
        num_triples = (trainval_df_filtered['dist'] == 3).sum()

        # Log the performance metrics and distribution of singles, doubles, and triples
        metrics_over_time.append({
            'iteration': i+1,
            'training_size': len(trainval_df_filtered),
            'mse': mse,
            'r2': r2,
            'spearman_rho': spearman_rho,
            'num_singles': num_singles,
            'num_doubles': num_doubles,
            'num_triples': num_triples
        })

        print(f"Iteration {i+1} - MSE: {mse}, R2: {r2}, Spearman Rho: {spearman_rho}, Training Size: {len(trainval_df_filtered)}")
        print(f"Distribution - Singles: {num_singles}, Doubles: {num_doubles}, Triples: {num_triples}")

        # Check for point of diminishing returns based on performance change
        if i > 0 and abs(metrics_over_time[-1]['spearman_rho'] - metrics_over_time[-2]['spearman_rho']) < 0.001:
            print(f"Point of diminishing returns at iteration {i+1}. Saving training and test datasets.")
            # Save the training and test datasets
            trainval_save_path = os.path.join(output_path, f"{name}_trainval_df_iteration_relaxed_{gpmap_type}_{i+1}.csv")
            test_save_path = os.path.join(output_path, f"{name}_test_df_iteration_relaxed_{gpmap_type}_{i+1}.csv")
            trainval_df_filtered.to_csv(trainval_save_path, index=False)
            test_df.to_csv(test_save_path, index=False)

        # Early stopping based on performance threshold
        if performance_threshold and r2 >= performance_threshold:
            print(f"Performance threshold reached (R2: {r2} >= {performance_threshold}). Stopping early.")
            break

    selected_seq_save_path = os.path.join(output_path, f"{name}_selected_seq.csv")
    pd.DataFrame(selected_sequences).to_csv(selected_seq_save_path, index=False)

    # Convert the metrics log into a DataFrame
    metrics_df = pd.DataFrame(metrics_over_time)

    # Save the metrics to the output directory
    metrics_save_path = os.path.join(output_path, f"{name}_active_learning_results_relaxed_{gpmap_type}.tsv")
    metrics_df.to_csv(metrics_save_path, sep="\t", index=False)

    return metrics_df


if __name__ == "__main__":
    # Command-line argument parsing
    parser = argparse.ArgumentParser(description='Active Learning Experiment')

    # Add arguments for the trainval_df, test_df, and output_path
    parser.add_argument('--name', type=str, required=True, help='Prefix for output files')
    parser.add_argument('--trainval_file', type=str, required=True, help='Path to the trainval_df file')
    parser.add_argument('--test_file', type=str, required=True, help='Path to the test_df file')
    parser.add_argument('--output_path', type=str, required=True, help='Directory to save output files')
    parser.add_argument('--gpmap_type', type=str, default='blackbox', help='GP map type to use')
    parser.add_argument('--n_queries', type=int, default=400, help='Number of queries to perform')
    parser.add_argument('--query_size', type=int, default=10, help='Number of samples to query in each iteration')
    parser.add_argument('--n_initial_doubles', type=int, default=10, help='Initial number of double mutants')
    parser.add_argument('--n_initial_triples', type=int, default=5, help='Initial number of triple mutants')
    parser.add_argument('--performance_threshold', type=float, default=None, help='Performance threshold for early stopping')

    # Parse the arguments
    args = parser.parse_args()

    # Run the experiment
    metrics_df = active_learning_experiment(
        name=args.name,
        gpmap_type=args.gpmap_type,
        trainval_file=args.trainval_file,
        test_file=args.test_file,
        output_path=args.output_path,
        n_queries=args.n_queries,
        query_size=args.query_size,
        n_initial_doubles=args.n_initial_doubles,
        n_initial_triples=args.n_initial_triples,
        performance_threshold=args.performance_threshold
    )

    # Print final metrics
    print(metrics_df)

