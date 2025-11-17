#!/bin/bash

# python monitor.py --scope_key scoped-keys/mnsol_fb-fit-v1-single-mean-k100_key.dat --restart
# python monitor.py --scope_key scoped-keys/fsolv_fb-fit-v1-single-mean-k100_key.dat --restart
# python monitor.py --scope_key scoped-keys/mnsol_fb-fit-v3-single-mean-k100_key.dat --restart
# python monitor.py --scope_key scoped-keys/fsolv_fb-fit-v3-single-mean-k100_key.dat --restart
# python monitor.py --scope_key scoped-keys/mnsol_fb-fit-v0-single-mean-k20_key.dat --restart
# python monitor.py --scope_key scoped-keys/fsolv_fb-fit-v0-single-mean-k20_key.dat --restart
# python monitor.py --scope_key scoped-keys/mnsol_fb-fit-v3-single-mean-k20_key.dat --restart
# python monitor.py --scope_key scoped-keys/fsolv_fb-fit-v3-single-mean-k20_key.dat --restart
python monitor.py --scope_key scoped-keys/mnsol_fit-iter-1_key.dat --restart

python monitor.py --scope_key scoped-keys/mnsol_vdw-refit_key.dat --restart
# python monitor.py --scope_key scoped-keys/fsolv_vdw-refit_key.dat --restart
