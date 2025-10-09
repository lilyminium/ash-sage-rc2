# Valence fitting

The subdirectories here handle steps involved with a ForceBalance valence fit, i.e. re-fitting bonds, angles, propers, impropers to QCArchive data (optimizations and torsiondrives).

The expected steps to be run first are:
- downloading QM data in ../01_download-data
- re-fitting vdW terms in ../02_fit-vdw


The steps here are:
- Generating new parameter splits, or modifying chemical partitioning by modifying SMARTS, in `01_generate-forcefield`
- Generating training and validation data set for broad coverage in `02_curate-data`
- Generating initial valence parameter values for bonds and angles in `03_generate-msm`
- Actually fitting the force field in `04_fit-forcefield`.
