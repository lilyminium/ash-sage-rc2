# vdW refit

This directory contains scripts for re-fitting a force field using pre-equilibrated data.

- `generate-forcefield.py` generates a force field without the ToolkitAM1BCC section, replacing it with NAGL charges. It saves it to `forcefield/force-field.offxml`
- `optimize.in` was copied from the Sage 2.0 fit. It specifies ForceBalance fitting options, which have not changed.

Note that to re-execute this fit, you need to copy over the `stored_data` directory generated in the `equilibration` directory.