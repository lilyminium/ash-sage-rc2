#!/bin/bash

# rsync -arv -P *.py run.sh hpc3:/dfs9/dmobley-lab/lilyw7/ash-sage-rc2/04_benchmark/phys-prop/
# rsync -arv -P mappings hpc3:/dfs9/dmobley-lab/lilyw7/ash-sage-rc2/04_benchmark/phys-prop/
rsync -arv -P --exclude="*.pkl" hpc3:/dfs9/dmobley-lab/lilyw7/ash-sage-rc2/04_benchmark/phys-prop/training .
rsync -arv -P --exclude="*.pkl" hpc3:/dfs9/dmobley-lab/lilyw7/ash-sage-rc2/04_benchmark/phys-prop/validation .