import numpy as np
import itertools
import heapq

# # Initialize an empty heap
# heap = []
#
# # Generate all possible 10bp sequences and compute total binding energy
# for sequence in itertools.product(range(4), repeat=10):
#     total_energy = sum(binding_energy_matrix[i, base] for i, base in enumerate(sequence))
#     # Use negative total_energy because heapq is a min-heap
#     heapq.heappush(heap, (-total_energy, "".join(base_dict[base] for base in sequence)))
#
# # Now the heap contains all sequences sorted by their total binding energy in descending order
# # The top 10 sequences and their total binding energies can be retrieved as follows:
# for _ in range(10):
#     energy, sequence = heapq.heappop(heap)
#     print(f"Sequence: {sequence}, Total binding energy: {-energy}")

def affinity_to_probability(affinity_matrix, l):
    # Initialize a list to store energy and sequence
    seq_and_energy = []

    # Create a dictionary of base pairs
    base_dict = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}

    # Generate all possible 10bp sequences and compute total binding energy
    for sequence in itertools.product(range(4), repeat=l):
        total_energy = sum(affinity_matrix[i, base] for i, base in enumerate(sequence))
        seq_and_energy.append((total_energy, "".join(base_dict[base] for base in sequence)))

    # Sort the list in descending order of energy
    seq_and_energy.sort(reverse=True)

    # Normalize the energy to convert it into probability
    energies = np.array([energy for energy, _ in seq_and_energy])
    probabilities = (energies - np.min(energies)) / (np.max(energies) - np.min(energies))
    return seq_and_energy, probabilities

if __name__ == "__main__":
    # Provided binding energy matrix
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

    l = 10

    seq_and_energy, probabilities = affinity_to_probability(binding_energy_matrix, l)

    # Write the results into a file
    with open('./SELEX_Data/propability.txt', 'w') as f:
        for (energy, sequence), probability in zip(seq_and_energy, probabilities):
            f.write(f'Sequence: {sequence}, Total binding energy: {energy}, Binding probability: {probability}\n')
