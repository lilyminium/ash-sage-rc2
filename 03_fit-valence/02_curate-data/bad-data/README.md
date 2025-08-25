# Filtering out data to ignore

QM data is excluded on the basis of the following reasons:

- AM1-BCC charges cannot be applied (`filter-out-smiles-by-charge.py`, `failed-charge-cmiles.dat`)
- NAGL charges cannot be applied (forbidden bond patterns in `smarts-to-exclude.dat`)
- seesaw SX4 geometries are difficult to get right with our current force field form (`filter-seesaw-sx4[-tds].py`, `failed-seesaw-sx4[-tds].dat`, also some previously determined in Sage 2.2.0 and copied to `existing-filtered-ids.dat`)
- bad 7-membered ring geometries in `existing-filtered-ids.dat` from Sage 2.2.0
