# xsim mk config
A config is stored as a `.toml` file. The main purpose of the config file is to
point to the files that store the "parsed and lightly processed" data from
CFOUR. Here is an example
```toml
# Default location of data needed by xsim.
states = "/path/to/mol_db.json"
# state_nicknames = "/path/to/nicknames.json"  # optional
# better_energies =   # optional
normal_modes = "/path/to/mol_normal_modes.json"
kappa_linear = "/path/to/all_kappa.json"
# kappa_quadratic = 
lambdas = "/path/to/all_lambdas.json"

active_states_and_modes = "/path/to/mol_active_modes.toml"
```

## files from the config
The `db` directory stores information on how to produce most the files
listed in the example above. The remaining files are listed below 

### `active_states_and_modes`
This file is used to describe the states and modes active in the KDC
simulation.

The `active_states_and_modes` file uses tables to define the KDC model
```toml
[first]
state_idxs = [1, 3]  
mode_idxs = [3, 11, 15, 20, 24]
```
The header is irrelevant, it serves almost as a comment. Mode numbering starts
at 1. The `active_states_and_modes` file can contain many tables
```toml
[first]
state_idxs = [1, 3]  
mode_idxs = [3, 11, 15, 20, 24]

[second]
state_idxs = [2]  
mode_idxs = [1, 7]
```
The active states in such model will be a union of the states defined in each
table. For the example above it will be `[1, 2, 3]`. Each table can miss any of
the `state_idxs` or `mode_idxs`.
```toml
[first]
state_idxs = [1, 3]  
mode_idxs = [3, 11, 15, 20, 24]

[extra_state]
state_idxs = [2]  

[extra_coupling]
mode_idxs = [1]
```

### `state_nicknames` (optional)
By default the `input_mk.py` script generates a figure showing the states and
couplings. If for any reason you would like the state names to be displayed as
something different create the `state_nicknames` file like so 
```json
[
  {
    "irrep": {
      "name": "A1",
      "energy #": 1
    },
    "nickname": "X"
  },
  {
    "irrep": {
      "name": "B1",
      "energy #": 1
    },
    "nickname": "A"
  }
]
```
As you can see above for each state all that matters is that `irrep` dictionary
stores the `name` and `energy #` that match the state for which name will be
changed.

### `mode_names_Mulliken` (optional)
Ab initio programs often use non-Mulliken orientation for the molecule's
geometry, but using the Mulliken convention is good. This file allows for a
manual assignment of mode names. Relavant only for plotting. Example file:
```json
[
  {
    "frequency, cm-1": 1000,
    "Mulliken": {
      "name": "A1",
      "number": 1
    }
  },
  {
    "frequency, cm-1": 700,
    "Mulliken": {
      "name": "A1",
      "number": 2
    }
  },
  {
    "frequency, cm-1": 800,
    "Mulliken": {
      "name": "B2",
      "number": 3
    }
  }
]
```
Modes are matched by their frequency.
