# Intro
This directory contains congested information that comes from parsed CFOUR 
output files. 

## Database
Various files in this directory, together, form the database.

### Recreate the database
Steps below allow to build a database that works as an input to the `input_mk`
script.

1. EOM states and TDM
a) Go to the directory that stores outputs of the EOM energy calculations and
run this:
```bash
cfour_parser -j output.c4 > output.json
processors/print_roots.py -j output.json | jq > eom.json
mv eom.json path/db/mol_at_geo/
```

b) Go to the TDM directory and generate TDMs:
```bash
# for each state
xcfour.py -j output.c4 > output.json
processors/print_tdm.py -j output.json | jq > tdm_state.json
```
Move all collected TDMs into the database and merge them into a single file:
```bash
mkdir path/db/mol_at_geo/tdms
# for each state
mv tdm_state.json path/db/mol_at_geo/tdms
cd path/db/mol_at_geo/tdms
cat tdm_*.json | jq -s '.' > mol_tdms.json
```

c) Inject the tdm values into the EOM states database:

Use `prepare.py` script, whose job is to add transition properties into the
database of EOM states: the addition works as this: 
- the TDM state matches `irrep #`, `eom model`, and `['energy']['total']['au']`
  of states from the database.
- The TMD values of `['energy']['transition']` are added to the matched state.
- The TDM state's transition properties are added to the matched state.

```bash
cd mol_at_geo
cp eom.json mol_db.json
prepare.py -j mol_db.json -f tdms/mol_tdms.json > tmp.json
mv tmp.json mol_db.json
```

2. Normal modes
```bash
cfour_parser -j output.c4 > output.json
processors/print_normal_coordinates.py -j output.json | jq > normal_coordinates.json
```

3. Lambdas 
```bash
cfour_parser -j out.c4 | jq > output.json
processors/print_gradient.py -j output.json | jq > lambda_state+state.json
mkdir db/lambdas
mv lambda_*.json db/lambdas
cd db/lambdas 
cat lambda_*.json | jq -s '.' > all_lambdas.json
```

4. Linear kappas

Literally the same as with lambdas. Use the `kappas_linear` name for the
directory.

5. Quadratic kappas
(automated -- questionable reliability)
Run the `generate_kappas_quadratic.sh` script 

(manual)
On each of the `xquadmodelout` files run 
```bash
kappa_quadratic.py -j <xquadmodelout> <state_name> | jq > kappa2nd_${a}.json
cp kappa2nd_${a}.json db/kappa_quadratic
```
After all states are ready, do
```bash
cd db/kappa_quadratic
cat kappa2nd_*.json | jq -s . > all_kappa2nd.json
```

## How to use the data
Data collected in this directory is used by other programs, most likely,
`input_mk.py` which has the job of creating an input file for the xsim program.
