# Compare AshGC training data to OpenFE FE benchmarks

The scripts and files in this directory compare the common core of ligands in each RBFE benchmarking system, following the standard OpenFE file structure. This expects protein-ligand targets with the following file structure:

```
.
└── <group name> (e.g. charge_annihilation_set)
        │
        └── <system name> (e.g. tyk2)
                │
                └── ligands.sdf (containing all system ligands)

```

The ligands in `ligands.sdf` are loaded and a common core is determined using RDKit's maximum common substructure algorithm. This core is then compared against the entire AshGC training set using OpenFF's substructure matching (powered either by OpenEye or RDKit) to look for matches. This gives somewhat more permissive results than looking for exact matches to a ligand.

The training set SMILES is uploaded here as a `.tar.gz`; untar to use with the script. An example of use is shown in `run-compare-training-overlap.sh` and an example output is in `compare-training-overlap.log`. 

On the OpenFE 2024 Industry benchmarks set, the below following systems and groups did not match the training set. Output is saved in `output/`. Core images are saved in `images/`.

- bayer_macrocycles:
    - brd4

- charge_annihilation_set:
    - cdk2
    - egfr
    - ephx2
    - irak4_s2
    - irak4_s3
    - itk
    - jnk1
    - ptp1b
    - tyk2

- jacs_set:
    - cdk2
    - jnk1
    - ptp1b
    - thrombin
    - tyk2

- janssen_bace:
    - bace_ciordia_prospective
    - bace_p3_arg368_in
    - keranen_p2

- mcs_docking_set:
    - renin

- merck:
    - cdk8
    - cmet
    - pfkfb3
    - syk- 

- miscellaneous_set:
    - cdk8
    - galectin
    - hiv1_protease

- scaffold_hopping_set:
    - factor_xa

- water_set:
    - brd4
    - chk1
    - hsp90_kung
    - taf12
    - thrombin
  