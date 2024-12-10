import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import gc
from utils import *

def run_sampling_experiment(name, trainval_file, test_file, output_path, gpmap_types, num_doubles, num_triples):
    # Read in trainval_df and test_df
    trainval_df = pd.read_csv(trainval_file, sep="\t")
    test_df = pd.read_csv(test_file, sep="\t")

    print(f"Running sampling for {name} with {num_doubles} doubles and {num_triples} triples.")

    # Precompute single training data
    single_training_data = trainval_df[trainval_df['dist'] == 1]

    # Validate available doubles and triples
    available_doubles = trainval_df[trainval_df['dist'] == 2].shape[0]
    available_triples = trainval_df[trainval_df['dist'] == 3].shape[0]

    if num_doubles > available_doubles or num_triples > available_triples:
        print(f"Insufficient doubles or triples. Skipping.")
        return

    # Sample doubles and triples
    new_double_training_data = trainval_df[trainval_df['dist'] == 2].sample(num_doubles)
    new_triple_training_data = trainval_df[trainval_df['dist'] == 3].sample(num_triples)
    new_dataset = pd.concat([single_training_data, new_double_training_data, new_triple_training_data]).reset_index(drop=True)

    # Split dataset into training and validation
    train_data, val_data = train_test_split(new_dataset, test_size=0.05, random_state=42)
    train_data['validation'] = False
    val_data['validation'] = True
    new_dataset = pd.concat([train_data, val_data])

    # Initialise results storage
    results_2 = []

    # Train and evaluate each G-P map model
    for gpmap_type in gpmap_types:
        print(f"Training model for gpmap_type={gpmap_type}")

        # Fit the model
        model = fit_model(new_dataset, gpmap_type)

        # Compute metrics
        mse, r2, spearman = calculate_metrics(model, test_df)

        # Generate and save plots
        save_model_plots(model, test_df, gpmap_type, output_path, name)

        # Store results
        results_2.append({
            'gpmap_type': gpmap_type,
            'mse': mse,
            'r2': r2,
            'spearman': spearman,
            'double_sample_size': num_doubles,
            'triple_sample_size': num_triples
        })

        # Clean up
        del model
        gc.collect()

    # Save results to a file
    results_2_df = pd.DataFrame(results_2)
    results_2_df.to_csv(os.path.join(output_path, f"{name}_nonlinear_sampling_results.tsv"), sep="\t", index=False)

def save_model_plots(model, test_df, gpmap_type, output_path, name):
    """Generates and saves model plots and data."""
    # Create directories for saving plots and data
    plot_dir = os.path.join(output_path, "plots")
    data_dir = os.path.join(output_path, "plot_data")
    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)

    # Prediction plot ($\hat{y}$ vs $y$)
    fig, ax = plt.subplots(figsize=[5, 5])
    y_test = test_df['y']
    yhat_test = model.x_to_yhat(test_df['x'])
    Rsq = np.corrcoef(yhat_test.ravel(), y_test)[0, 1] ** 2

    ax.scatter(yhat_test, y_test, color='C0', s=10, alpha=0.3, label='test data')
    xlim = [min(yhat_test), max(yhat_test)]
    ax.plot(xlim, xlim, '--', color='k', label='diagonal', zorder=100)
    ax.set_xlabel('model prediction ($\hat{y}$)')
    ax.set_ylabel('measurement ($y$)')
    ax.set_title(f'Standard metric of model performance:\n$R^2$={Rsq:.3f}')
    ax.legend()
    plt.tight_layout()
    plot_path = os.path.join(plot_dir, f"{name}_{gpmap_type}_nonlinear_prediction_plot.png")
    plt.savefig(plot_path)
    plt.close()

    # Save prediction data
    prediction_data = pd.DataFrame({"y_test": y_test, "yhat_test": yhat_test.ravel()})
    prediction_data.to_csv(os.path.join(data_dir, f"{name}_{gpmap_type}_nonlinear_prediction_data.csv"), index=False)

    # GE measurement process plot
    fig, ax = plt.subplots(figsize=[5, 5])
    phi_test = model.x_to_phi(test_df['x'])
    phi_lim = [min(phi_test) - 0.5, max(phi_test) + 0.5]
    phi_grid = np.linspace(phi_lim[0], phi_lim[1], 1000)
    yhat_grid = model.phi_to_yhat(phi_grid)
    yqs_grid = model.yhat_to_yq(yhat_grid, q=[0.025, 0.975])

    ax.fill_between(phi_grid, yqs_grid[:, 0], yqs_grid[:, 1], alpha=0.2, color='C1', lw=0, label='95% CI')
    ax.plot(phi_grid, yhat_grid, linewidth=3, color='C1', label='nonlinearity')
    ax.scatter(phi_test, y_test, color='C0', s=10, alpha=0.3, label='test data', zorder=100)
    ax.set_xlim(phi_lim)
    ax.set_xlabel('latent phenotype ($\phi$)')
    ax.set_ylabel('measurement ($y$)')
    ax.set_title('GE measurement process')
    ax.legend()
    plt.tight_layout()
    plot_path = os.path.join(plot_dir, f"{name}_{gpmap_type}_nonlinear_GE_process_plot.png")
    plt.savefig(plot_path)
    plt.close()

    # Save test data (phi_test) separately
    test_data = pd.DataFrame({
        "phi_test": phi_test,
        "y_test": y_test
    })
    test_data.to_csv(os.path.join(data_dir, f"{name}_{gpmap_type}_linear_test_data.csv"), index=False)

    # Save grid data (phi_grid, yhat_grid, yqs_grid) separately
    grid_data = pd.DataFrame({
        "phi_grid": phi_grid,
        "yhat_grid": yhat_grid,
        "yqs_lower": yqs_grid[:, 0],
        "yqs_upper": yqs_grid[:, 1],
    })
    grid_data.to_csv(os.path.join(data_dir, f"{name}_{gpmap_type}_nonlinear_grid_data.csv"), index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Double and Triple Sampling Experiment')

    # Add arguments
    parser.add_argument('--name', type=str, required=True, help='Prefix for output file names')
    parser.add_argument('--trainval_file', type=str, required=True, help='Path to the trainval_df file')
    parser.add_argument('--test_file', type=str, required=True, help='Path to the test_df file')
    parser.add_argument('--output_path', type=str, required=True, help='Directory to save output files')
    parser.add_argument('--gpmap_types', type=str, default="additive,neighbor,pairwise,blackbox", help='Comma-separated list of GP map types')
    parser.add_argument('--num_doubles', type=int, required=True, help='Number of doubles to sample')
    parser.add_argument('--num_triples', type=int, required=True, help='Number of triples to sample')

    args = parser.parse_args()

    gpmap_types = args.gpmap_types.split(',')

    run_sampling_experiment(
        name=args.name,
        trainval_file=args.trainval_file,
        test_file=args.test_file,
        output_path=args.output_path,
        gpmap_types=gpmap_types,
        num_doubles=args.num_doubles,
        num_triples=args.num_triples
    )

