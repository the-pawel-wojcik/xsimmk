#!/usr/bin/env python

import argparse
import math as m
import sys
from xsim.db.prepare import energies_match
from xsim.xsim_ids_processor import get_data_with_xsim_ids
from xsim.helpers.plot_overview import visualize_the_couplings

eV2cm = 8065.543937
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
    parser.add_argument(
        '-a',
        '--adiabatic_analysis',
        default=False,
        action='store_true',
        help='Add section that requests an analysis of the adiabatic'
        ' states.'
    )
    parser.add_argument(
        '--default_locations',
        type=str,
    )
    parser.add_argument(
        '-p',
        '--print_vectors',
        default=False,
        action='store_true',
        help='Add section that requests printing of solution vectors.'
             ' BE CAREFUL: produces gigantic files: LEVECS and LBASIS.',
    )
    parser.add_argument(
        '-v',
        '--visualize',
        default=True,
        action='store_false',
        help='By default this program produces a figure visualizing the'
             ' active states and modes of the model. Use this switch to'
             ' disable this functionality.',
    )
    parser.add_argument(
        '-f',
        '--show_frequencies',
        default=False,
        action='store_true',
        help='When -v, this flag allows to show mode frequencies together with'
        'couplings',
    )
    parser.add_argument(
        '--no_couplings',
        default=False,
        action='store_true',
        help='Prepare the input file but do NOT include in it the diabatic'
             ' couplings',
    )
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
        if 'Dipole strength' not in state.keys():
            print(
                "Error. Dipole strength not available for state:",
                file=sys.stderr
            )
            print(
                state,
                file=sys.stderr
            )
            xsim_input += "xxx "
        else:
            # # The oscillator strength version
            # f = state['Oscillator strength']
            # xsim_input += f"{f:.3f} "
            # The dipole strength version
            # # HINT: This is an old (incorrect) use of xsim
            # tdm_vec = state['Dipole strength']
            # tdm_square = sum([i**2 for i in tdm_vec.values()])
            # tdm = m.sqrt(tdm_square)
            # HINT: This is the new (correct) use of xsim
            #       use moments not squares
            dipole_strength = state['Dipole strength']
            dipole_strength_sum = sum([
                abs(i) for i in dipole_strength.values()
            ])
            moment = m.sqrt(dipole_strength_sum)
            xsim_input += f"{moment:.3f} "
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


def prepare_xsim_input_3rd_sec_gradients(linear_kappas) -> str:
    """ This section corresponds to energy gradients.
    Includes only state gradients (linear kappas)
    """
    xsim_input = "Linear\n"

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
    return xsim_input


def prepare_xsim_input_3rd_sec_couplings(lambdas) -> str:
    """ This section corresponds to energy gradients.
    Includes only interstate gradients (lambdas), i.e., linear diabatic
    coupling constants.
    """

    xsim_input = str()
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
    xsim_input = ""
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


def prepare_xsim_input_4th_sec_dibatic(quadratic_kappas_diabatic):
    """
    This section prepares the second energy derivatives (quadratic kappas)
    along the coupling modes, aka, along the diabatic part of potential that
    does not match the adiabatic part of the potential.
    """
    """
    state state mode mode coefficient
    """
    xsim_input = ""
    for k2 in quadratic_kappas_diabatic:
        state_ids = [state['ids']['xsim #'] for state in k2['EOM states']]

        # of no interest if any of the involved states is not part of the model
        if any([id < 1 for id in state_ids]):
            continue

        mode_ids = [mode['xsim #'] for mode in k2['normal modes']]
        if any([id < 1 for id in mode_ids]):
            continue

        amplitude = k2['kappa, cm-1']
        if energies_match(amplitude, 0.0, 0.1):
            continue

        xsim_input += " ".join([str(id) for id in state_ids]) 
        xsim_input += " "
        xsim_input += " ".join([str(id) for id in mode_ids])
        xsim_input += " "
        xsim_input += f"{amplitude:-7.1f}"
        xsim_input += "\n"

    return xsim_input


def prepare_xsim_input(data, no_couplings: bool = False) -> str:
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
    xsim_input += prepare_xsim_input_3rd_sec_gradients(data['linear kappas'])
    if no_couplings is False:
        xsim_input += prepare_xsim_input_3rd_sec_couplings(data['lambdas'])
    xsim_input += SECTION_SEP

    quadratic_kappas = data['quadratic kappas']
    quadratic_kappas_diabatic = data['quadratic kappas diabatic']
    if quadratic_kappas is not None or quadratic_kappas_diabatic is not None:
        xsim_input += "Quadratic\n"

    FULLY_SYMMETRIC = 'Ag'
    if quadratic_kappas is None:
        print("Info: No quadratic kappas available", file=sys.stderr)
    else:
        print(f"Info: Looking for {FULLY_SYMMETRIC} in "
              "quadratic_kappas kappas! Adjust for your point group.",
              file=sys.stderr)
        xsim_input += prepare_xsim_input_4th_sec(
            quadratic_kappas, FULLY_SYMMETRIC)

    if quadratic_kappas_diabatic is None:
        print("Info: No quadratic kappas diabatic available", file=sys.stderr)
    else:
        xsim_input += prepare_xsim_input_4th_sec_dibatic(
            quadratic_kappas_diabatic
        )

    if quadratic_kappas_diabatic is not None or quadratic_kappas is not None:
        xsim_input += SECTION_SEP

    xsim_input = xsim_input[:-1]  # remove the trailing newline
    return xsim_input


def main():
    args = get_args()
    data = get_data_with_xsim_ids(
        default_locations_filename=args.default_locations,
    )

    # TODO: it would be nice to know if there are lambdas or kappas missing
    # for a state that is active
    xsim_input = prepare_xsim_input(data, no_couplings=args.no_couplings)

    if args.adiabatic_analysis is True:
        xsim_input = xsim_input[:-1]  # remove the last new line
        xsim_input += ADIABATIC_ANALYSIS

    if args.print_vectors is True:
        xsim_input += '\n'
        xsim_input += "Storevectors\n"
        xsim_input += "Plotroots\n"
        xsim_input += "<n_roots>\n"
        xsim_input += "lines from fort.20 for each root"

    print(xsim_input)

    if args.visualize is True:
        eom_states = data['eom states']
        lambdas = data['lambdas']
        nmodes: list[dict] = data['nmodes']

        show_frequencies = args.show_frequencies
        save_svg = True
        visualize_the_couplings(
            eom_states=eom_states,
            normal_modes=nmodes,
            lambdas=lambdas,
            show_frequencies=show_frequencies,
            save=save_svg
        )


if __name__ == "__main__":
    main()
