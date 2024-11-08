# xsim mk
This directory contains scripts that enable an automated generation of an
input file for the `xsim` program of Dr. Stanton. 

The `input_mk.py` script expects as an input a parsed and lightly processed
output from CFOUR's calculations. The collection of these "parsed and lightly
processed ouputs" is referred to as a database and it is stored in the `db`
directory. The `db` directory contains among others the information on how to
re-create the database.

## How to use
1. Produce the database using instructions inside the `db` directory.
2. Go to the `config` directory and prepare a config file for your molecule.
3. Run `input_mk.py`.