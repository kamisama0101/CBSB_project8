from matplotlib import pyplot as plt

from SELEX_simulation import max_kmer_probability, complement

def max_kmer(sequences, kmer_probabilities):
    for i in range(len(sequences)-9):
        max_probability = max(kmer_probabilities.get(sequence[i:i + 10], 0))
        return max_probability

def get_key (dict, value):
    return list(dict.keys())[list(dict.values()).index(value)]



if __name__ == "__main__":
    # bar plot of kmer probability before and after SELEX simulation
    # calculate the kmer probabilities before SELEX simulation
    pre_simulation_probs = {}
    with open('./SELEX_Data/propability.txt', 'r') as f:
        for line in f:
            parts = line.strip().split(', ')
            kmer = parts[0][10:]
            probability = float(parts[2][20:])
            binding_eneergy = float(parts[1][21:])
            pre_simulation_probs[kmer] = probability
    # normalization
    total_pre_probability = sum(pre_simulation_probs.values())
    for kmer in pre_simulation_probs:
        pre_simulation_probs[kmer] /= total_pre_probability

    # print(list(kmer_binding_energy.items())[:5])

    # calculate the kmer probabilities after SELEX simulation
    count_sequence = {}
    sequences = []
    with open("./SELEX_Data/sequences.txt",'r') as file:
        next(file)
        for line in file:
            parts = line.strip().split('\t')
            sequence = parts[0][:]
            numbers = int(parts[1][:])
            count_sequence[sequence] = numbers
        file.close()

    # print(sequences[:5])

    post_simulation_probs ={}
    for sequence in count_sequence.keys():
        antisense_sequence = complement(sequence)
        max_probability = max(
            max_kmer_probability(sequence, pre_simulation_probs),
            max_kmer_probability(sequence[::-1], pre_simulation_probs),
            max_kmer_probability(antisense_sequence, pre_simulation_probs),
            max_kmer_probability(antisense_sequence[::-1], pre_simulation_probs)
        )
        kmer = get_key(pre_simulation_probs, max_probability)
        if kmer not in post_simulation_probs:
            post_simulation_probs[kmer] = count_sequence[sequence]
        else:
            post_simulation_probs[kmer] += count_sequence[sequence]

    total_post_probability = sum(post_simulation_probs.values())
    for kmer in post_simulation_probs:
        post_simulation_probs[kmer] /= total_post_probability

    # Sort k-mers by pre-simulation probabilities
    kmers = sorted(pre_simulation_probs.keys(), key=lambda k: pre_simulation_probs[k], reverse=True)

    # Limit to top N k-mers
    N = 50
    kmers = kmers[:N]

    pre_probs = [pre_simulation_probs[kmer] for kmer in kmers]
    post_probs = [post_simulation_probs.get(kmer, 0) for kmer in kmers]  # Use 0 if k-mer not in post_simulation_probs

    x = range(N)

    plt.figure(figsize=(12, 6))

    plt.bar(x, pre_probs, width=0.4, label='Pre-simulation')
    plt.bar(x, post_probs, width=0.4, label='Post-simulation', color='red', align='edge')

    plt.xlabel('k-mers')
    plt.ylabel('Probability')
    plt.title('k-mer probabilities before and after SELEX simulation')
    plt.xticks(x, kmers, rotation=90)  # Rotate labels for clarity
    plt.legend()

    plt.tight_layout()
    plt.show()




