#!/usr/bin/env python

import argparse
import math as m
import sys
from xsim.db.prepare import energies_match
from xsim.xsim_ids_processor import get_data_with_xsim_ids
from xsim.data_collector import get_state_nicknames
from vertical.energies_and_couplings.overview import visualize_the_couplings
from parsers.text import str_eom_state
from xsim.xsim_ids_processor import get_default_locations

eV2cm = 8065.54
ADIABATIC_ANALYSIS = """
Analyze
Convention
Diabatic
Optstate
2
Opttol
1
"""

CONICAL_INTERSECTION = """
Analyzer
Seamfind
"""


def get_args():
    parser = argparse.ArgumentParser(usage='Reads all data from config.toml')
    parser.add_argument('-a', '--adiabatic_analysis', default=False,
                        action='store_true', help='Add section that requests'
                        ' an analysis of the adiabatic states.')
    parser.add_argument('-p', '--print_vectors', default=False,
                        action='store_true', help='Add section that requests'
                        ' printing of solution vectors. BE CAREFUL: produces'
                        ' gigatinc files: LEVECS and LBASIS.')
    parser.add_argument('-v', '--visualize', default=True,
                        action='store_false', help='By default this program'
                        ' produces a figure visualizing the active states and'
                        ' modes of the model. Use this switch to disable this'
                        ' functionality.')
    args = parser.parse_args()
    return args


def prepare_xsim_input_1st_sec(states, modes, basis, lanczos):
    xsim_input = "Units\n"
    xsim_input += "1\n"
    xsim_input += "Dataset\n"
    xsim_input += "1\n"

    xsim_input += "States\n"
    xsim_input += f"{len(states)}\n"

    xsim_input += "Modes\n"
    xsim_input += f"{len(modes)}\n"

    xsim_input += "Basis Functions\n"
    for _ in modes:
        xsim_input += f"{basis} "
    xsim_input = xsim_input[:-1] + "\n"

    xsim_input += "Lanczos\n"
    xsim_input += f"{lanczos}\n"

    xsim_input += "Transition Moment\n"
    for state in states:
        if 'Norm of oscillator strength' not in state.keys():
            print("Error. Oscillator strength not available for state:")
            print(state)
            xsim_input += "xxx "
        else:
            # # The oscillator strength version
            # f = state['Oscillator strength']
            # xsim_input += f"{f:.3f} "
            # The dipole strength version
            tdm_vec = state['Dipole strength']
            tdm_square = sum([i**2 for i in tdm_vec.values()])
            tdm = m.sqrt(tdm_square)
            xsim_input += f"{tdm:.3f} "
    xsim_input = xsim_input[:-1] + "\n"

    xsim_input += "Vertical Energies\n"
    for state in states:
        if 'transition' not in state['energy'].keys():
            print("Error. transition energy not available for state:")
            print(state)
            xsim_input += f"{state['xsim #']} xxx\n"
        else:
            vertical = state['energy']['transition']['eV'] * eV2cm
            xsim_input += f"{state['xsim #']} {vertical:.2f}\n"

    return xsim_input


def prepare_xsim_input_2nd_sec(modes):
    """
    This section prints frequencies of the harmonic oscialltor basis.
    """
    xsim_input = "Harmonic Frequencies\n"
    for mode in modes:
        frequency = mode['frequency, cm-1']
        xsim_input += f"{mode['xsim #']} {frequency:.1f}\n"
    return xsim_input


def prepare_xsim_input_3rd_sec(data):
    """ This section corresponds to energy gradients.
    Includes both state gradients (linear kappas) and the interstate gradients
    (lambdas).
    """
    xsim_input = "Linear\n"

    linear_kappas = data['linear kappas']
    for kappa in linear_kappas:
        state_idx = kappa['EOM states'][0]['xsim #']

        if state_idx <= 0:
            continue

        for mode in kappa['gradient']:

            if mode['xsim #'] is None:
                print("Warning! Unable to process linear kappa for a mode with"
                      " an unkown xsim #; ɷ = "
                      f"{mode['frequency, cm-1']:.2f}", file=None)
                continue

            if mode['xsim #'] <= 0:
                continue

            amplitude = mode['gradient, cm-1']
            # don't bother clogging the input with zeros
            if energies_match(amplitude, 0.0, 0.1):
                continue

            xsim_input += f"{state_idx} {state_idx} {mode['xsim #']}"
            xsim_input += f" {amplitude:-7.1f}\n"

    lambdas = data['lambdas']
    for lmbd in lambdas:
        state_ids = [state['xsim #'] for state in lmbd['EOM states']]
        if len(state_ids) != 2:
            print("Error! There is a lambda with # EOM states != 2.",
                  file=sys.stderr)
            print(f"{lmbd['EOM states']=}", file=sys.stderr)
            sys.exit(1)

        # these states are not active in the simulation
        if state_ids[0] <= 0 or state_ids[1] <= 0:
            continue

        # these states ARE active
        for mode in lmbd['gradient']:

            if mode['xsim #'] is None:
                print("Warning! Unable to process lambda for a mode with"
                      " an unkown xsim #; ɷ = "
                      f"{mode['frequency, cm-1']:.2f}", file=sys.stderr)
                continue

            if mode['xsim #'] <= 0:
                continue

            amplitude = mode['gradient, cm-1']
            # don't bother clogging the input with zeros
            if energies_match(amplitude, 0.0, 0.1):
                continue

            xsim_input += f"{state_ids[0]} {state_ids[1]} {mode['xsim #']}"
            xsim_input += f" {amplitude:-7.1f}\n"

    return xsim_input


def prepare_xsim_input_4th_sec(quadratic_kappas, FULLY_SYMMETRIC):
    """
    This section prepares the second energy derivatives (quadratic kappas).
    The expansion coefficients are listed only along the fully symmetric modes.
    """
    xsim_input = "Quadratic\n"
    """
    state state mode mode coefficient
    """

    for k2 in quadratic_kappas:
        st_xsim = k2['state']['ids']['xsim #']
        pos2xsim = k2['pos2xsim#']
        pos2symmetry = k2['pos2symm']
        matrix = k2['kappa quadratic']
        if k2['state']['ids']['xsim #'] < 1:
            continue

        for row in range(len(matrix)):
            rowxsim = pos2xsim[row]
            row_symmetry = pos2symmetry[row]
            if rowxsim < 1:
                continue
            if row_symmetry != FULLY_SYMMETRIC:
                continue
            for col in range(row, len(matrix)):
                colxsim = pos2xsim[col]
                if colxsim < 1:
                    continue
                col_symmetry = pos2symmetry[col]
                if col_symmetry != FULLY_SYMMETRIC:
                    continue
                amplitude = matrix[row][col]
                # don't bother clogging the input with zeros
                if energies_match(amplitude, 0.0, 0.1):
                    continue
                xsim_input += f"{st_xsim} {st_xsim} {rowxsim} {colxsim}"
                xsim_input += f" {amplitude:-7.1f}"
                xsim_input += "\n"

    return xsim_input


def test_state_matched_better_state(state, ccsdt):
    """
    States match when "irrep"'s "name" and "energy #" match.
    """

    if state['irrep']['name'] != ccsdt['irrep']['name']:
        return False
    if state['irrep']['energy #'] != ccsdt['irrep']['energy #']:
        return False

    return True


def nickname_and_state_match(nickname: dict, state: dict) -> bool:
    """
    Helper for `inject_state_nicknames`.
    """
    if state['irrep']['name'] != nickname['irrep']['name']:
        return False

    if state['irrep']['energy #'] != nickname['irrep']['energy #']:
        return False

    return True


def inject_state_nicknames(states: dict, nicknames: dict):
    """
    Adds a "nickname" keyword to the "irrep" dictionary for each eom state.
    """

    # Print warnings for states without nicknames
    for state in states:
        state_has_a_nickname = False
        for nickname in nicknames:
            if nickname_and_state_match(nickname, state):
                state_has_a_nickname = True
                break
        if state_has_a_nickname is False:
            print(f"Warning: No nickname for {str_eom_state(state)}",
                  file=sys.stderr)

    # See if for every nickname there is exactly one matching target state
    for nickname in nicknames:
        nickname['matches'] = []
        for state in states:
            if not nickname_and_state_match(nickname, state):
                continue
            nickname['matches'] += [state]

    for nickname in nicknames:
        if len(nickname['matches']) == 0:
            print(f"Error: The nickname:\n\t{nickname}\n"
                  "\tdoes not correspond to any state.",
                  file=sys.stderr)
            continue
        if len(nickname['matches']) > 1:
            print("Warrning: The nickname:"
                  "\n\t"
                  f"{nickname}"
                  "\n\t"
                  "corresponds to more than one state.",
                  file=sys.stderr)
        for state in nickname['matches']:
            state['ids']['nickname'] = nickname['nickname']


def inject_better_energies(data):
    """
    Replaces the 'energy': 'transition': 'eV' value of the `active_states`
    with the correspondig "better energies" value found in `data`.

    The "better energies" are basically a vertical energies at the same
    geometry but at a better level of theory.

    Prints an error to stderr if the "better energies" value is missing, and
    leaves the original value unchanged.
    """

    if 'better energies' not in data:
        print("Error! Better energies not available.", file=sys.stderr)

    better_energies_states = data['better energies']
    states = data['eom states']
    for state in states:
        better_energy = None
        for better_state in better_energies_states:
            if not test_state_matched_better_state(state, better_state):
                continue
            better_energy = better_state['energy']['transition']['eV']
            break

        if better_energy is None:
            if state['ids']['xsim #'] >= 0:
                print('Missing "better" energy for an active state:\n\t'
                      + str_eom_state(state), file=sys.stderr)
            continue

        if 'transition' not in state['energy']:
            continue

        state['energy']['transition']['eV'] = better_energy


def prepare_xsim_input(data):
    basis = 15
    lanczos = 2000
    SECTION_SEP = "0 0 0 0 0 0 0 0 0 0 0 0\n"

    eom_states = data['eom states']
    active_states = [state for state in eom_states
                     if state['ids']['xsim #'] > 0]
    active_states.sort(key=lambda x: x['ids']['xsim #'])

    nmodes = data['nmodes']
    active_modes = [mode for mode in nmodes if mode['xsim #'] > 0]
    active_modes.sort(key=lambda x: x['xsim #'])

    xsim_input = prepare_xsim_input_1st_sec(
        active_states, active_modes, basis, lanczos)
    xsim_input += SECTION_SEP
    xsim_input += prepare_xsim_input_2nd_sec(active_modes)
    xsim_input += SECTION_SEP
    xsim_input += prepare_xsim_input_3rd_sec(data)
    xsim_input += SECTION_SEP

    FULLY_SYMMETRIC = 'Ag'
    quadratic_kappas = data['quadratic kappas']
    if quadratic_kappas is None:
        print("Warning: No quadratic kappas available", file=sys.stderr)
    else:
        print(f"Warning: Looking for {FULLY_SYMMETRIC} in "
              "quadratic_kappas kappas! Adjust for your point group.",
              file=sys.stderr)
        xsim_input += prepare_xsim_input_4th_sec(
            quadratic_kappas, FULLY_SYMMETRIC)
        xsim_input += SECTION_SEP

    return xsim_input


def main():
    args = get_args()
    default_locations = get_default_locations()
    data = get_data_with_xsim_ids()

    if 'better energies' in data:
        # TODO: move section that collectes "better energies" from
        # `get_data_xsim_ids` in here.
        print("Info: Using better energies", file=sys.stderr)
        inject_better_energies(data)

    # TODO: it would be nice to know if there are lambdas or kappas missing
    # for a state that is active
    xsim_input = prepare_xsim_input(data)

    if args.adiabatic_analysis is True:
        xsim_input = xsim_input[:-1]  # remove the last new line
        xsim_input += ADIABATIC_ANALYSIS

    if args.print_vectors is True:
        xsim_input += "Storevectors\n"
        xsim_input += "Plotroots\n"
        xsim_input += "<n_roots>\n"
        xsim_input += "lines from fort.20 for each root"

    print(xsim_input)

    if args.visualize is True:
        eom_states = data['eom states']
        lambdas = data['lambdas']
        save_svg = True
        if 'state_nicknames' in default_locations:
            print("Info: Using state nicknames.", file=sys.stderr)
            nicknames = get_state_nicknames(default_locations)
            inject_state_nicknames(eom_states, nicknames['state nicknames'])

        visualize_the_couplings(eom_states, lambdas, save=save_svg)


if __name__ == "__main__":
    main()
