# Filtering data

In this stage, we clean the data by filtering out fixing erroneous data, filtering out data with high-viscosity components, and running an initial filter for properties of interest.

## Erroneous data

In version 1.2 of the NIST ThermoML dataset (labelled v2020-09-30, last revised 2021-09-09), we discovered several cases of erroneous data. Full details are available in the `clean-data/` directory. They include data filtered out for several reasons:

- the mole fractions of binary mixtures were mixed up
- a high pressure set of experiments was entered as atmospheric pressure
- esters were incorrectly specified

## High viscosity components

We filtered out any data with a component over 0.3 Pascal seconds. This experimental data was obtained from ThermoML using the scripts in the `viscosities/scripts/` directory. The code in `viscosity_stub.py` was necessary to parse ThermoML for viscosity values and resulted in the files in `viscosities/`.

## Initial intermediate filter

The data was then filtered using all filters as laid out in `filter-data-initial.py`. 