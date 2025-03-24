# PP-QADMM
PP-QADMM: A Privacy-Preserving and Quantized ADMM Framework for Ensuring Privacy and Energy Efficiency in Federated Learning

# All main files run the algorithms and save results in mat files. Then corrsponding plot files are used to plot the generated results 

main.m: runs a comparison between PPQADMM, and considered baselines on i.i.d data distribution using California dataset (pre-processed, normalized, and clipped).

main_noniid.m: runs a comparison between PPQADMM, and considered baselines on a non i.i.d data distribution.

main_variousSigmas.m: Generate results for different choices for $\sigma_n^0$.

main_syntheticdata.m : runs a comparison between PPQADMM, and considered baselines on a synthetically generated dataset (10,0000 samples, 50 features each)

sensitivity_analysis.m : Can be edited to generate results for different number of workers, values for $\rho$, and $\epsilon$.



