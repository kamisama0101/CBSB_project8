import logomaker
import pandas as pd
import matplotlib.pyplot as plt

def nucleotide_frequencies(sequences):
    # Initialize a dataframe to store the counts
    df = pd.DataFrame(columns=['A', 'T', 'C', 'G'], index=range(1, len(sequences[0]) + 1)).fillna(0)

    # Loop through each sequence
    for sequence in sequences:
        # Loop through each base in the sequence
        for i, base in enumerate(sequence, start=1):
            df.at[i, base] += 1

    # Divide by the total number of sequences to get the frequencies
    df = df / len(sequences)

    return df

if __name__ == "__main__":
    # read the final sequences
    sequences = []
    with open("./SELEX_Data/sequences.txt", 'r') as file:
        next(file)
        for line in file:
            parts = line.strip().split('\t')
            sequence = parts[0][:]
            numbers = int(parts[1][:])
            for i in range(numbers):
                sequences.append(sequence)
        file.close()

    print(len(sequences))
    df = nucleotide_frequencies(sequences)
    print(df)
    df.to_csv("./SELEX_Data/base_probabilities.csv")
    # if we use this df to generate logo directly, it will report an error: logomaker.src.error_handling.LogomakerError: column 9 is '10' and has length 2; must have length 1.
    # Renaming columns to have single-character names
    # new_columns = {}
    # for col in df.columns:
    #     new_columns[col] = f"{col}A"
    #
    # df = df.rename(columns=new_columns)

    logomaker.Logo(df, font_name='Arial')

    plt.show()
