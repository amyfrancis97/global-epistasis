import os
import pandas as pd
import matplotlib.pyplot as plt

def generate_panel_of_plots(data_dir, output_path):
    prediction_files = [f for f in os.listdir(data_dir) if "prediction_data" in f]
    ge_files = [f for f in os.listdir(data_dir) if "GE_data" in f]

    for pred_file, ge_file in zip(prediction_files, ge_files):
        pred_data = pd.read_csv(os.path.join(data_dir, pred_file))
        ge_data = pd.read_csv(os.path.join(data_dir, ge_file))

        # Plot Prediction Data
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        ax[0].scatter(pred_data['yhat_test'], pred_data['y_test'], alpha=0.3)
        ax[0].plot([min(pred_data['yhat_test']), max(pred_data['yhat_test'])],
                   [min(pred_data['yhat_test']), max(pred_data['yhat_test'])], '--k')
        ax[0].set_title('Prediction Plot')
        ax[0].set_xlabel('$\hat{y}$')
        ax[0].set_ylabel('$y$')

        # Plot GE Measurement Data
        ax[1].fill_between(ge_data['phi_grid'], ge_data['yqs_lower'], ge_data['yqs_upper'], alpha=0.2, color='C1')
        ax[1].plot(ge_data['phi_grid'], ge_data['yhat_grid'], color='C1')
        ax[1].scatter(ge_data['phi_test'], pred_data['y_test'], alpha=0.3)
        ax[1].set_title('GE Measurement Process')
        ax[1].set_xlabel('$\phi$')
        ax[1].set_ylabel('$y$')

        # Save Combined Plot
        plot_name = pred_file.replace('_prediction_data.csv', '_combined_plot.png')
        plt.tight_layout()
        plt.savefig(os.path.join(output_path, plot_name))
        plt.close()

# Call function
generate_panel_of_plots(data_dir="path_to/plot_data", output_path="path_to/panel_plots")

