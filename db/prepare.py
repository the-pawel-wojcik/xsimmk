#!/usr/bin/env python3

import argparse
import json
import math as m
import sys


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('-f', '--tdm',
                        help='JSON file with TDM and the like data.')
    parser.add_argument('-j', '--json', default=False,
                        action='store_true', help='Print final db in JSON.')
    parser.add_argument('-v', '--verbose', default=False, action="store_true")
    args = parser.parse_args()
    return args


def energies_match(one: float, two: float, tolerance: float = 1e-5):
    """
    Returns True if the energies are the same.
    The energies are the same if their difference is smaller than tolerance.
    If energies are in a.u. the tolerance of 1e-5 a.u. is about 0.3 meV.
    """
    return m.fabs(one-two) < tolerance


def states_are_different(one, two):
    """
    Returns False if irrep #, model and energy (up to 1e-5 au) match.

    -- paramerters
    `one` and `two` are states. A minimal state looks like this:

    {
     "irrep": {"#": 8},
     "model": "EOMEE-CCSD",
     "energy": {
         "total": {
             "au": -263.72252668779663
             }
         }
     }


    """
    one_model = one['model']
    two_model = two['model']
    if two_model != one_model:
        return True

    one_irrep_no = one['irrep']['#']
    two_irrep_no = two['irrep']['#']
    if two_irrep_no != one_irrep_no:
        return True

    one_energy = one['energy']['total']['au']
    two_energy = two['energy']['total']['au']
    if not energies_match(one_energy, two_energy):
        return True

    return False


def add_tdm_to_db(db, tdm, verbose=False):
    """
    Example of the db file:
    [{
        "model": "EOMEE-CCSD",
        "irrep": {
          "#": 5,
          "# of roots": 2,  # new part
          "name": "B1u",  # new part
          "energy #": 1
        },
        "converged root": { "singles": [...], "doubles": [...], }, # new part
        "energy": {
          "excitation": {
            "au": 0.161843216187947,
            "eV": 4.403977806379308
          },
          "total": {
            "au": -263.75449282684275
          }
        },
        "ids": {
          "#": 8,
          "energy #": 0
        }
    },...]
    Example of the tdm file:
    [{
      "irrep": { "#": 1 },
      "energy": {
        "total": { "au": -263.5574035814958 },
        "transition": {
          "eV": 9.767,
          "nm": 126.9414,
          "cm-1": 78776.525,
        }
      },
      "model": "EOMEE-CCSD",
      "Right TDM": { "x": -0, "y": 0, "z": 0 },
      "Left TDM": { "x": -0, "y": 0, "z": 0 },
      "Dipole strength": { "x": 0, "y": 0, "z": 0 },
      "Oscillator strength": { "x": 0, "y": 0, "z": 0 },
      "Norm of oscillator strength": 0
    },]
    """

    for state in db:
        if states_are_different(state, tdm):
            continue

        st_number = state['ids']['#']
        if verbose is True:
            print(f"Found TDMs for a database state # {st_number}. ", end='')
            print(f"It's an {state['model']} state from irrep #"
                  f"{state['irrep']['#']} {state['irrep']['energy #']}"
                  f"{state['irrep']['name']}")

        transition_energy = tdm['energy']['transition']
        state['energy']['transition'] = transition_energy
        del state['energy']['excitation']

        keys = ['Right TDM', 'Left TDM', 'Dipole strength',
                'Oscillator strength', 'Norm of oscillator strength',]
        for key in keys:
            state[key] = tdm[key]

        return

    tdm_state_str = f"{tdm['model']} state in irrep #{tdm['irrep']['#']}"
    print(f"Warning! A TDM available for {tdm_state_str} "
          "does not match any state from the database.", file=sys.stderr)


def main():
    args = get_args()
    with open(args.database, 'r') as database:
        db = json.load(database)

    if args.tdm is not None:
        with open(args.tdm, 'r') as tdm_input:
            tdms = json.load(tdm_input)

        for tdm in tdms:
            add_tdm_to_db(db, tdm, args.verbose)

    if args.json is True:
        print(json.dumps(db))

    return


if __name__ == "__main__":
    main()
