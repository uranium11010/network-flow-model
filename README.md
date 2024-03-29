# Empirically testing a MaxEnt network flow model for a bipartite food web

This is one of two repositories associated with the paper
"Nestedness Promotes Stability in Maximum-Entropy Bipartite Food Webs" by Li and Harte, 2023.
The other repository is located [https://github.com/uranium11010/nestedness-stability](here).

If you use our code in your work, please cite our paper:
```
@article{li2023nestedness,
  title={Nestedness Promotes Stability in Maximum-Entropy Bipartite Food Webs},
  author={Li, Zhening and Harte, John},
  year={2023}
}
```

## Dependencies

All code was implemented in Python 3.

To install the required dependencies, run the following command in your command line:
```
pip install numpy scipy matplotlib xlrd==1.2.0 pandas xlsxwriter
```

## Code files

Each `.py` file is a standalone file that can be run on its own. Please see the detailed descriptions below.

### graphing+numerical_sol.py
Contains code for
1. Graphing log10(data/model) for specific networks (toggle line 199)
2. Generating MaxEnt model flow rates through numerial solution of Eq. 6ab (toggle line 200); stored in file `log10_data_model.txt`; `histogram.py` draws histogram from these values (Fig. 4 of Main Text). See section _Empirical Testing_ of Main Text.
3. Generating random networks with the same topology and aggregate flow rates as the Glen Canyon network (toggle line 201); values of Diatom-Simuliidae pair are stored in file `glen_diat_simul_random.txt`. See section _Discussion: Discovering Interaction-Specific Mechanisms_ of Main Text.

### energy_flow_data
Contains the data files that `graphing+numerical_sol.py` imports for graphing and comparison with model predictions.

### log10_data_model.txt
Stores all 210 values of log10(data/model) (data = empirical flow rate; model = MaxEnt model flow rate) generated from `graphing+numerical_sol.py`; used by `histogram.py` to generate a histogram (Fig. 4 of Main Text). See section _Empirical Testing_ of Main Text.

### histogram.py
Generates a histogram from values of log10(data/model) (Fig. 4 of Main Text); also generates some statistics (e.g. RMS, Q1, Q3). See section _Empirical Testing_ of Main Text.

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
