import math
import numpy as np
from iglm import IgLM

iglm = IgLM()

sequence = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKEYYMHWVRQAPGKGLEWVGLIDPEQGNTIYDPKFQDRATISADNSKNTAYLQMNSLRAEDTAVYYCARDTAAYFDYWGQGTLVTVS"
chain_token = "[HEAVY]"
species_token = "[HUMAN]"

# Define the 20 standard amino acids
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

# Initialise an array to store log-likelihoods of mutated sequences
sequence_length = len(sequence)
likelihood_matrix = np.zeros((sequence_length, len(amino_acids)))

# Calculate the log-likelihood for the wild-type sequence
wild_type_log_likelihood = iglm.log_likelihood(sequence, chain_token, species_token)

# Iterate over each position in the sequence
for i in range(sequence_length):
    # Iterate over each possible amino acid mutation
    for j, aa in enumerate(amino_acids):
        # Create a mutant sequence by replacing the amino acid at position i
        mutant_sequence = sequence[:i] + aa + sequence[i+1:]
        
        # Calculate the log-likelihood for the mutant sequence
        mutant_log_likelihood = iglm.log_likelihood(mutant_sequence, chain_token, species_token)
        
        # Store the log-likelihood for this mutation in the matrix
        likelihood_matrix[i, j] = mutant_log_likelihood

# To convert log-likelihoods to probabilities, normalise the values
# Subtracting the log-likelihood of the wild-type sequence ensures probabilities are relative to it
probability_matrix = np.exp(likelihood_matrix - wild_type_log_likelihood)

# Print or process the probability matrix
print(probability_matrix)

