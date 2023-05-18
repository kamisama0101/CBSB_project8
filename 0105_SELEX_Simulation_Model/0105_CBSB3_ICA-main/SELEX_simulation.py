import numpy as np
import random
import multiprocessing as mp
import time

# Function to generate a random k bp sequence
def generate_kbp_sequence(k):
    return ''.join(random.choice('ATCG') for _ in range(k))

def complement(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[base] for base in sequence)

def max_kmer_probability(sequence, kmer_probabilities):
    return max(kmer_probabilities.get(sequence[i:i + 10], 0) for i in range(len(sequence) - 9))

# Function to process a single sequence
def process_sequence(_, k, kmer_probabilities):
    sequence = generate_kbp_sequence(k)
    antisense_sequence = complement(sequence)
    max_probability = max(
        max_kmer_probability(sequence, kmer_probabilities),
        max_kmer_probability(sequence[::-1], kmer_probabilities),
        max_kmer_probability(antisense_sequence, kmer_probabilities),
        max_kmer_probability(antisense_sequence[::-1], kmer_probabilities)
    )
    return (sequence, max_probability)

def SELEX_simulation(k, S, T, r, kmer_probabilities):
    # Set up a multiprocessing pool and generate and process S sequences
    with mp.Pool() as pool:
        sequence_probabilities = pool.starmap(process_sequence, [(i, k, kmer_probabilities) for i in range(S)])

    # Normalize the sequence probabilities so they sum to 1
    total_probability = sum(probability for sequence, probability in sequence_probabilities)
    sequence_probabilities = [(sequence, probability / total_probability) for sequence, probability in sequence_probabilities]

    # Repeat the selection and replication process for r rounds
    for round in range(r):
        # Select T sequences according to their probabilities
        sequences = [sequence for sequence, probability in sequence_probabilities]
        probabilities = [probability for sequence, probability in sequence_probabilities]
        selected_sequences = np.random.choice(sequences, size=T, p=probabilities)

        # Simulate the PCR process to replicate the selected sequences back up to S
        pcr_probabilities = [1 / len(selected_sequences)] * len(selected_sequences)
        sequence_counts = np.random.multinomial(S, pcr_probabilities)

        # Update the sequence probabilities based on the maximum probability kmer within each sequence
        sequence_probabilities = [(sequence, probability) for sequence, probability in sequence_probabilities if sequence in selected_sequences]

        # Normalize the updated sequence probabilities so they sum to 1
        total_probability = sum(probability for sequence, probability in sequence_probabilities)
        sequence_probabilities = [(sequence, probability / total_probability) for sequence, probability in sequence_probabilities]

    final_sequences_pool = {sequence: count for sequence, count in zip(sequences, sequence_counts)}
    # Randomly select 1 million sequences from the final pool
    selected_sequences = random.choices(list(final_sequences_pool.keys()), k=1000000)

    # Count the number of each sequence type
    sequence_type_counts = {sequence: selected_sequences.count(sequence) for sequence in set(selected_sequences)}
    return sequence_type_counts


if __name__ == "__main__":
    time1 = time.time()

    # # Generate a seed
    # seed = random.randint(0, 10 ** 9)

    seed = 139948
    # Set the seed
    random.seed(seed)

    print(f"Random seed: {seed}")

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

    # Run the simulation with parameters k=16, S=1000000000, T=8000000, r=3 and find the most frequent sequence in results
    sequence_type_counts = SELEX_simulation(15, 10000000, 100000, 3, kmer_probabilities)
    most_frequent_sequence = max(sequence_type_counts, key=sequence_type_counts.get)
    num_frequent_sequence = sequence_type_counts[most_frequent_sequence]

    # write the results in sequences.txt
    with open('./SELEX_Data/sequences.txt', 'w') as file:
        file.write(f'The most frequent sequence is: {most_frequent_sequence}\t {num_frequent_sequence}\n')
        for sequence, count in sequence_type_counts.items():
            file.write(f'{sequence}\t{count}\n')
        file.close()

    time2 = time.time()
    print("Program time:", time2-time1)



