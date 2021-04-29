# network-flow-model
Contains code used in the analysis of the project on modeling flow distribution in a bipartite flow network

## graphing+numerical_sol.py
Contains code for
1. Graphing log10(data/model) for specific networks (toggle line 199)
2. Generating MaxEnt model flow rates through numerial solution of Eq. 6ab (toggle line 200); stored in file "log10_data_model.txt"
3. Generating random networks with the same topology and aggregate flow rates as the Glen Canyon network (toggle line 201); values of Diatom-Simuliidae pair are stored in file "glen_diat_simul_random.txt"

## log10_data_model.txt
Stores all 210 values of log10(data/model) (data = empirical flow rate; model = MaxEnt model flow rate).

## glen_diat_simul_random.txt
Stores 1000 values of the random model flow rate for the Diatom-Simuliidae pair of the Glen Canyon network.

## histogram.py
Code for generating a histogram from values of log10(data/model); also generates some statistics (e.g. RMS, Q1, Q3).
