# utils.py

import gc
import pandas as pd
import numpy as np
import mavenn
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.model_selection import train_test_split


# Function to calculate evaluation metrics
def calculate_metrics(model, test_dataset):
    y_test = test_dataset['y']
    x_test = test_dataset['x']    

    # Compute yhat on test data
    yhat_test = model.x_to_yhat(x_test)

    # Compute R^2 between yhat_test and y_test
    r2 = np.corrcoef(yhat_test.ravel(), y_test)[0, 1]**2
    mse = mean_squared_error(y_test, yhat_test.ravel())
    # Compute Spearman's Rho
    spearman_rho, _ = spearmanr(yhat_test.ravel(), y_test)

#    mse = mean_squared_error(y_test, y_pred)
    #r2 = r2_score(y_test, y_pred)

    #spearman_rho, _ = spearmanr(y_test, y_pred)

    return mse, r2, spearman_rho

# Function to fit the model
def fit_model(trainval_df, gpmap_type):
    index_first = trainval_df.index.tolist()[0]
    L = len(trainval_df.loc[index_first, 'x'])

    # Initialize the model
    model = mavenn.Model(
        L=L,
        alphabet='protein*',
        regression_type='GE',
        gpmap_type=gpmap_type,
        ge_nonlinearity_type='nonlinear',
        ge_nonlinearity_monotonic=False
    )

    # Set the data
    model.set_data(x=trainval_df['x'], y=trainval_df['y'], validation_flags=trainval_df['validation'])

    # Train the model and capture history
    model.fit(epochs=200, batch_size=32, learning_rate=0.0007, verbose=False, early_stopping=True, early_stopping_patience=2, linear_initialization=True)

    return model

# Function to compute uncertainty using Monte Carlo Dropout
def compute_uncertainty_mc_dropout(model, remaining_pool, num_passes=10):
    x_remaining = remaining_pool['x'].tolist()
    all_predictions = []
    for _ in range(num_passes):
        phi_remaining = model.x_to_phi(x_remaining)
        y_pred = model.phi_to_yhat(phi_remaining)
        all_predictions.append(y_pred)
    all_predictions = np.array(all_predictions)
    uncertainty = np.std(all_predictions, axis=0)
    return uncertainty

# Function to plot metrics vs sample size
def plot_sample_size_metrics(results, models=['additive', 'neighbor', 'pairwise', 'blackbox']):
    metrics = ['mse', 'r2', 'spearman']
    palette = sns.color_palette("husl", len(models))
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for j, metric in enumerate(metrics):
        ax = axes[j]
        for i, model in enumerate(models):
            model_data = results[results['gpmap_type'] == model][metric]
            sample_size = results[results['gpmap_type'] == model]['sample_size']
            sorted_indices = np.argsort(sample_size)
            sample_size_sorted = sample_size.iloc[sorted_indices]
            model_data_sorted = model_data.iloc[sorted_indices]
            ax.plot(sample_size_sorted, model_data_sorted, label=f'{model.capitalize()} Model', color=palette[i])
        ax.set_xlabel('Sample Size', fontsize=12)
        ax.set_ylabel(metric.capitalize(), fontsize=12)
        ax.set_title(f'{metric.capitalize()} vs Sample Size', fontsize=14)
        ax.axvline(x=1000, color='grey', linestyle='--')
        ax.axvline(x=2500, color='grey', linestyle='--')
        ax.axvline(x=5000, color='grey', linestyle='--')
    plt.tight_layout()
    plt.legend(models)
    plt.show()

# Function to compute uncertainty based on residuals
def compute_uncertainty(model, remaining_pool):
    x_remaining = remaining_pool['x'].tolist()
    y_true = remaining_pool['y'].values  # Ensure you have access to true labels in remaining_pool

    # Compute predictions
    y_pred = model.x_to_yhat(x_remaining)

    # Compute residuals
    residuals = np.abs(y_true - y_pred.ravel())
    
    return residuals  # Larger residuals indicate higher uncertainty

# Function to fit the model
def fit_model_linear_lp_mapping(trainval_df, gpmap_type, linear=False):
    """
    Fits a MAVE-NN model to the given training/validation dataset.
    
    Parameters:
        trainval_df (pd.DataFrame): Training and validation dataset.
        gpmap_type (str): The type of GP map to use ('additive', 'neighbor', 'pairwise', 'blackbox').
        linear (bool): Whether to enforce a linear mapping for the latent phenotype vs measurement.
    
    Returns:
        model: A trained MAVE-NN model.
    """
    index_first = trainval_df.index.tolist()[0]
    L = len(trainval_df.loc[index_first, 'x'])

    # Default model settings
    default_model_kwargs = {
        'L': L,
        'alphabet': 'protein*',
        'regression_type': 'GE',
        'ge_noise_model_type': 'SkewedT',
        'ge_heteroskedasticity_order': 2
    }

    # Default fitting settings
    default_fit_kwargs = {
        'learning_rate': .001,
        'epochs': 5,
        'batch_size': 200,
        'early_stopping': True,
        'early_stopping_patience': 30,
        'linear_initialization': False,
        'verbose': False
    }

    # Adjust settings based on gpmap_type
    if gpmap_type == 'pairwise':
        default_fit_kwargs['batch_size'] = 50
    elif gpmap_type == 'blackbox':
        default_model_kwargs['gpmap_kwargs'] = {'hidden_layer_sizes': [10] * 5, 'features': 'pairwise'}
        default_fit_kwargs['learning_rate'] = 0.0005
        default_fit_kwargs['batch_size'] = 50
        default_fit_kwargs['early_stopping_patience'] = 10

    # Adjust for linear mapping if specified
    if linear:
        default_model_kwargs['ge_nonlinearity_type'] = 'linear'
        default_model_kwargs['ge_nonlinearity_monotonic'] = False  # Linear maps are inherently monotonic
    else:
        default_model_kwargs['ge_nonlinearity_type'] = 'nonlinear'  # Default to sigmoid or other types

    # Initialize the model with the specified settings
    model = mavenn.Model(gpmap_type=gpmap_type, **default_model_kwargs)
    model.set_data(x=trainval_df['x'], y=trainval_df['y'], validation_flags=trainval_df['validation'])
    model.fit(**default_fit_kwargs)

    return model

