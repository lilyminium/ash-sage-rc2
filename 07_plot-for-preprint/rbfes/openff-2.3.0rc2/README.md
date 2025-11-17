# OpenFF 2.3.0rc2 benchmarks

## Overview

Relative binding free energy benchmarks of the openff-2.3.0rc2 release using the JACS set from (https://github.com/OpenFreeEnergy/IndustryBenchmarks2024).

### Contents

#### Benchmark variants

The following benchmark variants were run:
  * `openff-2.3.0rc2_elf10`: simulations using OpenFF 2.3.0rc2 and AM1BCCELF10 charges using the OpenEye Toolkit.
  * `openff-2.3.0rc2_nagl`: simulations using OpenFF 2.3.0rc2 and AshGC v0.1rc3 charges

For additional validation, an extra 3 replica were run for both variants across all sets. Results have been collated in separate
files which are labelled with the tag `6replica` in their name.

#### RBFE data

The `data` folder contains all the ddG data for each transformation in each set.


##### Per system data files

These are stored as TSV files within subdirectories for each system in the JACS set.
Reference CSV data (labelled as `{system}_ref.csv`), taken from the Ross et al. 2023 dataset (https://github.com/schrodinger/public_binding_free_energy_benchmark)
can also be found.

##### Aggregate dG and ddG files

CSV files with all the dG and ddG values for all datasets can be found in `data`.

### Simulation details

* All simulations were carried out using the RelativeHybridTopology protocol in OpenFE v1.4.
* Mappings were created using Kartograf, these are the same as those used in the OpenFE industry benchmarks.
* Transformation edges were created using Lomap, these are the same as those used in the OpenFE industry benchmarks.
* Simulations used openfe default settings, except the following:
  * Dodecahedron solvent box with 1 nm padding in complex and 1.5 nm padding in solvent
  * 2.5 ps HREX exchange rate
  * 0.9 nm short range interaction cutoff
  * 3 replica were originally run, and then a further 3 were added to the `6replica` files.

