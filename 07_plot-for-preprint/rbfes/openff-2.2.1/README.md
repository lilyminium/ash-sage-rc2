# OpenFF 2.2.1 benchmarks

## Overview

Relative binding free energy benchmarks of the openff-2.2.1 release using the JACS set from (https://github.com/OpenFreeEnergy/IndustryBenchmarks2024).

### Contents

#### Benchmark variants

The following benchmark variants were run:
  * `openff-2.2.1_elf10`: simulations using OpenFF 2.2.1 and AM1BCCELF10 charges using the OpenEye Toolkit.
  * `openff-2.2.1_nagl`: simulations using OpenFF 2.2.1 and AshGC v0.1rc3 charges
  * `openff-2.2.1_am1bcc`: simulations using OpenFF 2.2.1 and AM1BCC charges calculated using ambertools. Please note that these charges were not assigned to the input conformer, but rather a conformer that was generated at charge assignment time via RDKit through the OpenFF toolkit.

#### RBFE data

The `data` folder contains all the ddG data for each transformation in each set. These are stored as TSV files within subdirectories for each system in the JACS set.

##### Per system data files

These are stored as TSV files within subdirectories for each system in the JACS set.
Reference CSV data (labelled as `{system}_ref.csv`), taken from the Ross et al. 2023 dataset (https://github.com/schrodinger/public_binding_free_energy_benchmark)
can alsos be found.

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

