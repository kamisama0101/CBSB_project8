# 0105_CBSB3_ICA
This repository is a sub-project of the CBSB3 course ICA, offered by ZJE-UOE institute. The aim is to reproduce and optimize the mathematical model of SELEX technique developed by Juli et al, 2012 via python tools. 

**If you want to come across the source paper, please click the superlink:['Juli et al, 2012'](https://arxiv.org/abs/1205.1819)**
>**workflow for replication**
>
> 1. Extract the data from your SELEX experiment and preprocess it, including aligning the oligonucleotides and applying any necessary transformation.
> 2. Define the Gibbs free energy equation and the conditional probability distribution based on the SELEX data.
> 3. Use Monte Carlo simulation to estimate the denominator in conditional probability distribution.
> 4. Use the downhill simplex method to optimize the likelihood function of the model, which is the negative of the log-likelihood.
> 5. Extract the Gibbs free energy matrix, which represents the affinity of DNA binding sequences to the transcription factor.
> 6. Test the SELEX model with simulated and real SELEX data to validate its performance.

## Change of Plan
Given that we could not get experimental SELEX data. We change our plan from replication of Juli et al.'s model to simulation the whole SELEX workflow, simplistically.

 >**workflow for simulation**
>
> 1. The total binding energy of all possible 10bp kmer is calculated using the base-position Gibbs free energy matrix, and it is normalized as the selection probability of kmer.
> 2. Randomly generate a billion sequences of 15bp length as the initial oligonucleotides pool.
> 3. The kmer corresponding to each sequence is found, and the extraction probability of kmer is regarded as the extraction probability of this sequence
> 4. One million sequences were selected according to the extraction probability to simulate the process of combining targets with sequences.
> 5. The 1 million selected sequences were amplified back to 1 billion in a simulated PCR replication. 
> 6. The extraction and PCR replication processes were repeated three times to simulate the three rounds of enrichment in the actual experiment.
> 7. A million sequences were randomly selected (not based on extraction probability) to simulate the final sequencing process in the SELEX experiment.
> 8. visualize the final sequences results via probabilities bar plot and sequence logo. 

### Program manual
install dependency
```sh
python -m pip install -r requirements.txt 
```

run the main program of the models (Noted that it will take up a lot of memory and time)
```sh
python main.py 
```

If you want to visualize the sequencing data, you can run the two visualization programs
```sh
# bar plot
python visualization_bar.py
# sequences logo (PWM)
python visualization_logo.py
```


