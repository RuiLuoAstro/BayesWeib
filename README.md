# BayesWeib
A simple Bayesian package to calculate the repeating burst rate based on the Weibull distribution

## References

If you would like to use this code to study the non-Possionian process of repeating FRBs, please cite the paper [Luo et al. 2020, Nature, 586, 693](https://ui.adsabs.harvard.edu/abs/2020Natur.586..693L/abstract).

## Mathematical Basis

Please refer to the statisical derivation for non-Possionian process in the paper [Oppermann et al. 2018, MNRAS, 475, 5109](https://ui.adsabs.harvard.edu/abs/2018MNRAS.475.5109O/abstract)

## Dependencies

Python (3.6.x), PyMultiNest (see https://github.com/JohannesBuchner/PyMultiNest for more details)

## Usage

Run by ``` ./run_nest.sh ``` 

**Notes: Sampling is very fast, no need to run it in parallel on any computer cluster. The posterior outputs will be saved on the directory ```./nest_out/```**


## Plot making

Available in ```pltpost.ipynb```
