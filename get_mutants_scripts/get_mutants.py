# %%
import pandas as pd
import random
import math
import sys
import os
import argparse

# %%
# Create a list containing the AA alphabet
amino_acids = [
    'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
]

# %%
def create_single_mutants(WT_CDR, amino_acids):
    """
    Create all possible single mutants for a given wild-type CDR sequence.

    Parameters:
    WT_CDR (str): Wild-type CDR sequence.
    amino_acids (list): List of amino acids.

    Returns:
    list: List of single mutant sequences.
    int: Number of single mutant sequences.
    """
    mutants = []
    for i in range(len(WT_CDR)):
        for aa in amino_acids:
            if WT_CDR[i] != aa:
                mutant = WT_CDR[:i] + aa + WT_CDR[i+1:]
                mutants.append(mutant)
    
    number_mutants = len(mutants)
    print(f"Number of single mutant sequences: {number_mutants}")
    
    return mutants, number_mutants

# %%
def create_double_triple_mutants(WT_CDR, amino_acids, mutation_type='double'):
    """
    Create a single double or triple mutant for a given wild-type CDR sequence.

    Parameters:
    WT_CDR (str): Wild-type CDR sequence.
    amino_acids (list): List of amino acids.
    mutation_type (str): Type of mutation ('double' or 'triple').

    Returns:
    str: Mutant sequence.
    """
    num_mutations = 2 if mutation_type == 'double' else 3
    positions = random.sample(range(len(WT_CDR)), num_mutations)
    mutant = list(WT_CDR)
    for pos in positions:
        new_aa = random.choice([aa for aa in amino_acids if aa != WT_CDR[pos]])
        mutant[pos] = new_aa
    
    return ''.join(mutant)

# %%
def generate_mutants(WT_CDR, num_mutants, amino_acids, mutation_type="double"):
    """
    Generate a specified number of double or triple mutants for a given wild-type CDR sequence.

    Parameters:
    WT_CDR (str): Wild-type CDR sequence.
    num_mutants (int): Number of mutants to generate.
    amino_acids (list): List of amino acids.
    mutation_type (str): Type of mutation ('double' or 'triple').

    Returns:
    list: List of mutant sequences.
    """
    mutants = set()
    while len(mutants) < num_mutants:
        mutants.add(create_double_triple_mutants(WT_CDR, amino_acids, mutation_type))
    
    mutants = list(mutants)
    print(f"Number of {mutation_type} mutant sequences generated: {len(mutants)}")
    return mutants

# %%
def get_mutant_df(data_list):
    """
    Create a DataFrame from a list of mutant sequences with their respective mutation types.

    Parameters:
    data_list (list): List of tuples containing mutant sequences and their mutation type.

    Returns:
    DataFrame: Combined DataFrame of mutants.
    """
    mutant_dfs = []
    for data, m_type in data_list:
        df = pd.DataFrame(data, columns=['sequence'])
        df["dist"] = m_type
        mutant_dfs.append(df)
        
    combined_mutants = pd.concat(mutant_dfs, ignore_index=True)
    return combined_mutants

# %%

def calculate_mutants(n):
    """
    Calculate the total number of possible double and triple mutants for a given sequence length.

    Parameters:
    n (int): Length of the amino acid sequence.

    Returns:
    int: Number of possible double mutants.
    int: Number of possible triple mutants.
    """
    # 19 possible mutations at each position (excluding the wild-type amino acid)
    double_mutants = 0 if n < 2 else math.comb(n, 2) * 19 * 19
    triple_mutants = 0 if n < 3 else math.comb(n, 3) * 19 * 19 * 19

    return double_mutants, triple_mutants

# %%
def get_absolut_input(CDR_11_mer, perc_double_mut, perc_triple_mut, filename, output_dir):
    """
    Generate mutants and save them in specified file formats for Absolut! input.

    Parameters:
    CDR_11_mer (str): Wild-type CDR sequence.
    perc_double_mut (float): Percentage of double mutants to generate.
    perc_triple_mut (float): Percentage of triple mutants to generate.
    filename (str): Base name for output files.
    output_dir (str): Directory to save the output files.
    """
    WT_CDR = CDR_11_mer

    single_mutants, number_single_mutants = create_single_mutants(WT_CDR, amino_acids)

    n = len(WT_CDR)
    double_mutants_num, triple_mutants_num = calculate_mutants(n)

    print(f"Total number of double mutants for sequence length {n}: {double_mutants_num}")
    print(f"Total number of triple mutants for sequence length {n}: {triple_mutants_num}")

    num_double_mutants = int(double_mutants_num * perc_double_mut)
    double_mutants = generate_mutants(WT_CDR, num_double_mutants, amino_acids, 'double')

    num_triple_mutants = int(triple_mutants_num * perc_triple_mut)
    triple_mutants = generate_mutants(WT_CDR, num_triple_mutants, amino_acids, 'triple')

    data_list = [(single_mutants, 'single'), (double_mutants, 'double'), (triple_mutants, 'triple')]

    all_mutants = get_mutant_df(data_list)
    all_mutants = all_mutants.rename(columns={'sequence': 'x'})

    mapping = {'single': 1, 'double': 2, 'triple': 3}
    all_mutants['dist'] = all_mutants['dist'].map(mapping)

    print(all_mutants)

    output_path = os.path.join(output_dir, filename)
    all_mutants['x'].to_csv(f'{output_path}.txt', header=None, sep="\t")
    all_mutants.to_csv(f'{output_path}_with_dist.txt', sep="\t", index=False)

    # Check if the file was created and is not empty
    if os.path.exists(f'{output_path}.txt') and os.path.getsize(f'{output_path}.txt') > 0:
        print(f"Mutant file generated successfully: {f'{output_path}.txt'}")
    else:
        print(f"Failed to generate the mutant file: {f'{output_path}.txt'}")

# %%
def main():
    parser = argparse.ArgumentParser(description="Generate CDR mutants for Absolut! input.")
    parser.add_argument('--sequence', type=str, help="Wild-type CDR sequence.")
    parser.add_argument('--perc_double_mutants', type=float, help="Percentage of double mutants to generate.")
    parser.add_argument('--perc_triple_mutants', type=float, help="Percentage of triple mutants to generate.")
    parser.add_argument('--filename', type=str, help="Base name for output files.")
    parser.add_argument('--output_dir', type=str, default='.', help="Directory to save output files.")

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    get_absolut_input(args.sequence, args.perc_double_mutants, args.perc_triple_mutants, args.filename, args.output_dir)



# %%
if __name__ == '__main__':
    main()
