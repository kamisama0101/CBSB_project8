import random
import numpy as np
import itertools
from affinity import affinity_to_probability
from SELEX_simulation import SELEX_simulation

if __name__ == "__main__":
    seed = 139948
    # Set the seed
    random.seed(seed)

    binding_energy_matrix = np.array([
        [-4.722516, -5.729347, 0.000000, -6.251779],
        [-7.447426, -5.981440, 0.000000, -16.853690],
        [0.000000, -6.946246, -15.701235, -8.529272],
        [-7.746046, -15.548042, -12.535315, 0.000000],
        [-7.989755, -7.201358, -24.708969, 0.000000],
        [0.000000, -9.611195, -8.497223, -5.336888],
        [-0.505663, -19.926999, 0.000000, -4.445374],
        [-1.836787, -0.228140, 0.000000, -0.945140],
        [-1.841359, -1.612913, 0.000000, -1.417988],
        [-1.431632, -1.539663, 0.000000, -0.235633]
    ])

    # the length of kmer
    l = 10
    # the length of oligonucleotides (sequences)
    k = 15
    # the amount of oligonucleotides pool
    S = 100000000
    # the amount of targets
    T = 1000000
    # the round of SELEX
    r = 3

    seq_and_energy, probabilities = affinity_to_probability(binding_energy_matrix, l)

    # Write the results into a file
    with open('./SELEX_Data/propability.txt', 'w') as f:
        for (energy, sequence), probability in zip(seq_and_energy, probabilities):
            f.write(f'Sequence: {sequence}, Total binding energy: {energy}, Binding probability: {probability}\n')

        # Read the results file and create a dictionary mapping 10bp sequences to their probabilities
    kmer_probabilities = {}
    with open('./SELEX_Data/propability.txt', 'r') as f:
        for line in f:
            parts = line.strip().split(', ')
            kmer = parts[0][10:]
            probability = float(parts[2][20:])
            binding_eneergy = float(parts[1][21:])
            kmer_probabilities[kmer] = probability

    # Normalize the probabilities so they sum to 1
    total_probability = sum(kmer_probabilities.values())
    for kmer in kmer_probabilities:
        kmer_probabilities[kmer] /= total_probability


    sequence_type_counts = SELEX_simulation(k, S, T, r, kmer_probabilities)
    most_frequent_sequence = max(sequence_type_counts, key=sequence_type_counts.get)
    num_frequent_sequence = sequence_type_counts[most_frequent_sequence]

    # write the results in sequences.txt
    with open('./SELEX_Data/sequences.txt', 'w') as file:
        file.write(f'The most frequent sequence is: {most_frequent_sequence}\t {num_frequent_sequence}\n')
        for sequence, count in sequence_type_counts.items():
            file.write(f'{sequence}\t{count}\n')
        file.close()
