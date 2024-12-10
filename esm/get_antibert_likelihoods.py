import numpy as np
import pandas as pd
import torch
from antiberty import AntiBERTyRunner
from torch.nn.functional import softmax

# Initialise the AntiBERTy model
antiberty = AntiBERTyRunner()

# Wild-type sequence and 11-mer subsequence
wild_type_sequence = ("QVQLQQPGAELVKPGASVKLSCKASGYTFTSYWMHWVKQGPGQGLEWIGEIDPSDSYPNYNEKFKGKATLTVDKSSSTAYMQLSSLTSEDSAVYYCAS" +
                      "LYYYGTSYGVL" +
                      "DYWGQGTSVTVSSAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSPRPSETVT" +
                      "CNVAHPASSTKVDKKIVP")

# Define the 11-mer subsequence and find its position in the full sequence
target_11mer = "LYYYGTSYGVL"
start_index = wild_type_sequence.index(target_11mer)  # Find start index of the 11-mer
end_index = start_index + len(target_11mer)  # End index of the 11-mer

# Define the 20 standard amino acids
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

# Initialise lists to store data for the DataFrame
mutated_sequences = []
raw_plls = []
relative_plls = []
normalised_likelihood_matrix = np.zeros((len(target_11mer), len(amino_acids)))

# Step 1: Calculate the pseudo-log-likelihood for the wild-type sequence
wild_type_log_likelihood = antiberty.pseudo_log_likelihood([wild_type_sequence], batch_size=16)

# Step 2: Iterate over positions within the 11-mer
for i in range(len(target_11mer)):
    position_in_full_sequence = start_index + i  # Map to full sequence position
    position_log_likelihoods = []

    # Step 3: For each position, iterate over all amino acids
    for j, aa in enumerate(amino_acids):
        # Create a mutant sequence by substituting the amino acid at position `start_index + i`
        mutant_sequence = (wild_type_sequence[:position_in_full_sequence] + aa +
                           wild_type_sequence[position_in_full_sequence + 1:])
        
        # Calculate the pseudo-log-likelihood for the mutant sequence
        pll_mutant = antiberty.pseudo_log_likelihood([mutant_sequence], batch_size=16).item()
        position_log_likelihoods.append(pll_mutant)  # Store the raw PLL for this mutation
        
        # Step 4: Calculate relative (normalised) log-likelihood
        relative_pll = pll_mutant - wild_type_log_likelihood.item()
        normalised_likelihood_matrix[i, j] = relative_pll  # Store in the matrix
        
        # Step 5: Store mutated 11-mer sequence, raw PLL, and relative PLL
        mutated_sequences.append(mutant_sequence[start_index:end_index])  # Just the 11-mer part
        raw_plls.append(pll_mutant)  # Raw PLL
        relative_plls.append(relative_pll)  # Relative PLL

# Step 6: Convert the relative log-likelihood matrix to probabilities using softmax
probability_matrix = np.exp(normalised_likelihood_matrix)

# Step 7: Create a DataFrame to store the results
df = pd.DataFrame({
    'Mutated Sequence': mutated_sequences,
    'Raw PLL': raw_plls,
    'Relative PLL': relative_plls
})

# Print the DataFrame and probability matrix
print(df)
print(probability_matrix)

# Save the DataFrame as a CSV file
df.to_csv('mutated_sequences_pll.csv', index=False)

print("DataFrame saved as 'mutated_sequences_pll.csv'")

# Save the probability matrix as an .npy file
np.save('probability_matrix.npy', probability_matrix)

print("Probability matrix saved as 'probability_matrix.npy'")

