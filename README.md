# network-flow-model
Contains code used in the analysis of the project on modeling flow distribution in a bipartite flow network

### graphing+numerical_sol.py
Contains code for
1. Graphing log10(data/model) for specific networks (toggle line 199)
2. Generating MaxEnt model flow rates through numerial solution of Eq. 6ab (toggle line 200); stored in file "log10_data_model.txt"
3. Generating random networks with the same topology and aggregate flow rates as the Glen Canyon network (toggle line 201); values of Diatom-Simuliidae pair are stored in file "glen_diat_simul_random.txt"

### energy_flow_data
Contains the data files that "graphing+numerical_sol.py" imports for graphing and comparison with model predictions.

### log10_data_model.txt
Stores all 210 values of log10(data/model) (data = empirical flow rate; model = MaxEnt model flow rate).

### glen_diat_simul_random.txt
Stores 1000 values of the random model flow rate for the Diatom-Simuliidae pair of the Glen Canyon network.

### histogram.py
Generates a histogram from values of log10(data/model); also generates some statistics (e.g. RMS, Q1, Q3).

### random_glen.py
Generates some statistics (mean, stdev, min, max) of random/model; hypothesis testing for how much empirical flow rate of Diatom-Simuliidae of the Glen Canyon network deviates from model prediction.

### example_networks.py
Generates model flow rates using Eq. 7b for example networks; the graphs are Fig. 3a-d of Main Text.

## References

Nicknames for the food webs (for naming the corresponding energy flow data files) are given in square brakets.

1.	\[Glen Canyon\] W. F. Cross, C. V. Baxter, K. C. Donner, E. J. Rosi-Marshall, T. A. Kennedy, R. O. Hall, Jr., H. A. Wellard Kelly, R. S. Rogers, Ecosystem ecology meets adaptive management: food web response to a controlled flood on the Colorado River, Glen Canyon. *Ecol. Appl.* **21**, 2016–2033 (2011).
2.	\[Forest streams\] P. M. Bumpers, A. D. Rosemond, J. C. Maerz, J. P. Benstead, Experimental nutrient enrichment of forest streams increases energy flow to predators along greener food-web pathways. *Freshw. Biol.* **62**, 1794–1805 (2017).
3.	\[Coral reef\] E. C. de la Morinière, B. J. A. Pollux, I. Nagelkerken, M. A. Hemminga, A. H. L. Huiskes, G. van der Velde, Ontogenetic dietary changes of coral reef fishes in the mangrove-seagrass-reef continuum: stable isotopes and gut-content analysis. *Mar. Ecol. Prog. Ser.* **246**, 279–289 (2003).
4.	\[Crustacean\] D. Rudnick, V. Resh, Stable isotopes, mesocosms and gut content analysis demonstrate trophic differences in two invasive decapod crustacea. *Freshw. Biol.* **50**, 1323–1336 (2005).
5.	\[Earthworm\] M. Judas, Gut content analysis of earthworms (Lumbricidae) in a beechwood. *Soil Biol. Biochem.* **24**, 1413–1417 (1992).
6.	\[Chesapeake\] D. Baird, R. E. Ulanowicz, The seasonal dynamics of the Chesapeake Bay ecosystem. *Ecol. Monogr.* **59**, 329–364 (1989).
7.	\[Caribbean\] L. Goldwasser, J. Roughgarden, Construction and analysis of a large Caribbean food web. *Ecology* **74**, 1216–1233 (1993).
8.	\[Seagrass\] R. R. Christian, J. J. Luczkovich, Organizing and understanding a winter’s seagrass foodweb network through effective trophic levels. *Ecol. Modell.* **117**, 99–124 (1999).
