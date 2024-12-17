import torch
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
import seaborn as sns
from esm import pretrained

# Load a single model
def load_single_model(model_name):
    model, alphabet = pretrained.load_model_and_alphabet(model_name)
    model.eval()
    return model, alphabet

# Function to calculate log-probabilities for a model
def calculate_log_probs(sequence, model, alphabet):
    batch_converter = alphabet.get_batch_converter()
    data = [("sequence", sequence)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)

    with torch.no_grad():
        token_probs = torch.log_softmax(model(batch_tokens)["logits"], dim=-1)
    return token_probs.squeeze().cpu().numpy()

# Function to calculate relative log-probabilities for a subsequence
def calculate_relative_probs_in_subsequence(sequence, subsequence, model_log_probs, alphabet, amino_acids):
    start_pos = sequence.find(subsequence)
    if start_pos == -1:
        raise ValueError("Subsequence not found in the sequence.")
    end_pos = start_pos + len(subsequence)

    valid_indices = [alphabet.get_idx(aa) for aa in amino_acids]
    relative_log_probs = []

    for i in range(start_pos, end_pos):
        pos_relative_probs = []
        for j, mutant_aa in enumerate(amino_acids):
            mutant_prob = model_log_probs[i + 1, valid_indices[j]]
            pos_relative_probs.append(mutant_prob)

        relative_log_probs.append(pos_relative_probs)

    return np.array(relative_log_probs).T  # Transpose for heatmap compatibility

# Step 3: Generate single mutants with higher likelihood than the wild type
def generate_single_mutants(pseudo_wt_subsequence, relative_log_probs_pseudo_wt, amino_acids):
    single_mutants = []
    likelihood_differences = []  # Store likelihood differences for each single mutant
    
    for i, aa in enumerate(pseudo_wt_subsequence):
        pseudo_wt_prob = relative_log_probs_pseudo_wt[amino_acids.index(aa), i]
        for j, mutant_aa in enumerate(amino_acids):
            if relative_log_probs_pseudo_wt[j, i] > pseudo_wt_prob:  # Check if mutant AA has a higher likelihood
                mutant_seq = list(pseudo_wt_subsequence)
                mutant_seq[i] = mutant_aa
                mutant_seq_str = "".join(mutant_seq)
                
                # Calculate likelihood difference for this mutation
                likelihood_diff = relative_log_probs_pseudo_wt[j, i] - pseudo_wt_prob
                single_mutants.append(mutant_seq_str)
                likelihood_differences.append(likelihood_diff)
    
    return single_mutants, likelihood_differences

# Step 4: Generate double mutants
def generate_double_mutants(single_mutant, relative_log_probs_single_mutant, amino_acids, original_mutated_pos):
    double_mutants = []
    cumulative_likelihood_differences = []
    
    for i, aa in enumerate(single_mutant):
        if i == original_mutated_pos:  # Skip the already mutated position
            continue
        
        single_mutant_prob = relative_log_probs_single_mutant[amino_acids.index(aa), i]
        for j, mutant_aa in enumerate(amino_acids):
            if relative_log_probs_single_mutant[j, i] > single_mutant_prob:  # Check if mutant AA has a higher likelihood
                mutant_seq = list(single_mutant)
                mutant_seq[i] = mutant_aa
                mutant_seq_str = "".join(mutant_seq)
                
                # Calculate likelihood difference for this double mutation
                likelihood_diff = relative_log_probs_single_mutant[j, i] - single_mutant_prob
                double_mutants.append(mutant_seq_str)
                cumulative_likelihood_differences.append(likelihood_diff)
    
    return double_mutants, cumulative_likelihood_differences

# Step 5: Generate triple mutants
def generate_triple_mutants(double_mutant, relative_log_probs_double_mutant, amino_acids, mutated_positions):
    triple_mutants = []
    cumulative_likelihood_differences = []
    
    for i, aa in enumerate(double_mutant):
        if i in mutated_positions:  # Skip already mutated positions
            continue
        
        double_mutant_prob = relative_log_probs_double_mutant[amino_acids.index(aa), i]
        for j, mutant_aa in enumerate(amino_acids):
            if relative_log_probs_double_mutant[j, i] > double_mutant_prob:  # Check if mutant AA has a higher likelihood
                mutant_seq = list(double_mutant)
                mutant_seq[i] = mutant_aa
                mutant_seq_str = "".join(mutant_seq)
                
                # Calculate likelihood difference for this triple mutation
                likelihood_diff = relative_log_probs_double_mutant[j, i] - double_mutant_prob
                triple_mutants.append(mutant_seq_str)
                cumulative_likelihood_differences.append(likelihood_diff)
    
    return triple_mutants, cumulative_likelihood_differences

# Step 6: Create a DataFrame with pseudo-WT, single, double, and triple mutants
def create_dataframe(pseudo_wt_subsequence, single_mutants, double_mutants, triple_mutants, single_likelihood_diffs, double_likelihood_diffs, triple_likelihood_diffs):
    sequences = [pseudo_wt_subsequence] + single_mutants + double_mutants + triple_mutants
    distances = [0] + [1] * len(single_mutants) + [2] * len(double_mutants) + [3] * len(triple_mutants)
    
    likelihood_diffs = [0] + single_likelihood_diffs + double_likelihood_diffs + triple_likelihood_diffs  # 0 for pseudo-WT
    
    df = pd.DataFrame({'sequence': sequences, 'dist': distances, 'likelihood_diff': likelihood_diffs})
    return df

# Step 7: Replace the CDR sequence in the full sequence with the pseudo-WT
def replace_subsequence_in_full_sequence(full_sequence, original_subsequence, new_subsequence):
    return full_sequence.replace(original_subsequence, new_subsequence)

def main(model_names, iterations=1, least_likely=False):

    wild_type_sequence = "QVQLQQPGAELVKPGASVKLSCKASGYTFTSYWMHWVKQGPGQGLEWIGEIDPSDSYPNYNEKFKGKATLTVDKSSSTAYMQLSSLTSEDSAVYYCAS" + \
                         "LYYYGTSYGVL" + \
                         "DYWGQGTSVTVSSAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSPRPSETVT" + \
                         "CNVAHPASSTKVDKKIVP"
    wt_subsequence = "LYYYGTSYGVL"  # CDR for wild-type sequence

    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

    # List to store combined DataFrames for mutations across all models
    combined_dfs = []

    # Iterate over multiple models
    for model_name in model_names:
        print(f"Processing model: {model_name}")

        # Load the single ESM model
        model, alphabet = load_single_model(model_name)

        if least_likely:
            # Step 1: Calculate log-probabilities for the WT sequence
            log_probs_wt = calculate_log_probs(wild_type_sequence, model, alphabet)

            # Step 1: Calculate relative log-probabilities for the WT subsequence
            relative_log_probs_wt = calculate_relative_probs_in_subsequence(wild_type_sequence, wt_subsequence, log_probs_wt, alphabet, amino_acids)

            # Step 1: Find least likely mutations for the WT subsequence
            least_likely_mutations = find_least_likely_mutations(relative_log_probs_wt, amino_acids)

        # Generate pseudo-WT sequences and corresponding mutants
        for iteration in range(iterations):
            if least_likely:
                # Generate the pseudo-WT sequence by mutating four random least likely positions
                pseudo_wt_subsequence, selected_positions = generate_pseudo_wt_randomized(
                    wt_subsequence, relative_log_probs_wt, least_likely_mutations, amino_acids, num_mutations=4)

                # Replace the subsequence in the full sequence with the pseudo-WT sequence
                pseudo_wt_full_sequence = replace_subsequence_in_full_sequence(
                    wild_type_sequence, wt_subsequence, pseudo_wt_subsequence)

            else:
                pseudo_wt_full_sequence = wild_type_sequence
                pseudo_wt_subsequence = wt_subsequence

            # Step 3: Calculate log-probabilities for the pseudo-WT sequence
            log_probs_pseudo_wt = calculate_log_probs(pseudo_wt_full_sequence, model, alphabet)

            # Step 3: Calculate relative log-probabilities for the pseudo-WT subsequence
            relative_log_probs_pseudo_wt = calculate_relative_probs_in_subsequence(
                pseudo_wt_full_sequence, pseudo_wt_subsequence, log_probs_pseudo_wt, alphabet, amino_acids)

            # Step 4: Generate single mutants where the likelihood is higher than the pseudo-WT
            single_mutants, single_likelihood_diffs = generate_single_mutants(pseudo_wt_subsequence, relative_log_probs_pseudo_wt, amino_acids)

            # Step 5: Generate double mutants
            double_mutants, double_likelihood_diffs = [], []
            for single_mutant in single_mutants:
                original_mutated_pos = next(i for i, (a, b) in enumerate(zip(single_mutant, pseudo_wt_subsequence)) if a != b)
                single_mutant_full_sequence = replace_subsequence_in_full_sequence(
                    wild_type_sequence, wt_subsequence, single_mutant)
                log_probs_single_mutant = calculate_log_probs(single_mutant_full_sequence, model, alphabet)
                relative_log_probs_single_mutant = calculate_relative_probs_in_subsequence(
                    single_mutant_full_sequence, single_mutant, log_probs_single_mutant, alphabet, amino_acids)

                new_double_mutants, new_double_likelihood_diffs = generate_double_mutants(
                    single_mutant, relative_log_probs_single_mutant, amino_acids, original_mutated_pos)
                double_mutants.extend(new_double_mutants)
                double_likelihood_diffs.extend(new_double_likelihood_diffs)

            # Step 6: Generate triple mutants
            triple_mutants, triple_likelihood_diffs = [], []
            for double_mutant in double_mutants:
                mutated_positions = [i for i, (a, b) in enumerate(zip(double_mutant, pseudo_wt_subsequence)) if a != b]
                double_mutant_full_sequence = replace_subsequence_in_full_sequence(
                    wild_type_sequence, wt_subsequence, double_mutant)
                log_probs_double_mutant = calculate_log_probs(double_mutant_full_sequence, model, alphabet)
                relative_log_probs_double_mutant = calculate_relative_probs_in_subsequence(
                    double_mutant_full_sequence, double_mutant, log_probs_double_mutant, alphabet, amino_acids)

                new_triple_mutants, new_triple_likelihood_diffs = generate_triple_mutants(
                    double_mutant, relative_log_probs_double_mutant, amino_acids, mutated_positions)
                triple_mutants.extend(new_triple_mutants)
                triple_likelihood_diffs.extend(new_triple_likelihood_diffs)

            # Step 6: Create a DataFrame with the pseudo-WT, single, double, and triple mutants
            df = create_dataframe(pseudo_wt_subsequence, single_mutants, double_mutants, triple_mutants, single_likelihood_diffs, double_likelihood_diffs, triple_likelihood_diffs)

            # Append the results to the combined DataFrame for this model
            combined_dfs.append(df)

    # Combine all DataFrames from different models and iterations
    combined_df = pd.concat(combined_dfs).drop_duplicates(subset=['sequence']).reset_index(drop=True)

    # Save the combined DataFrame
    combined_df.to_csv("combined_pseudo_mutants.csv", sep="\t", header=True)
    combined_df.to_pickle("combined_pseudo_mutants.pkl")
    
    return combined_df

model_names = [
    'esm1b_t33_650M_UR50S',
    #'esm1v_t33_650M_UR90S_1',
    #'esm1v_t33_650M_UR90S_2',
    #'esm1v_t33_650M_UR90S_3',
    #'esm1v_t33_650M_UR90S_4',
    #'esm1v_t33_650M_UR90S_5',
]

result = main(model_names, iterations=1, least_likely=False)
print(result)
