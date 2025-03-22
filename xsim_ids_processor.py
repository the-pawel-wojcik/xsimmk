"""
Database readers.
"""

import os
import sys
import tomllib
from xsim.db.prepare import states_are_different, energies_match
from xsim.data_collector import get_eom_states, get_normal_modes, \
    get_linear_kappas, get_lambdas, get_active_states_and_modes, \
    get_kappas_quadratic, get_better_energies, get_Mulliken_mode_names, \
    get_state_nicknames, get_quadratic_kappas_diabatic
from parsers.text import str_eom_state
from typing import Any

# modes have the same frequency if they differ by less than FREQ_TOL cm-1
FREQ_TOL = 0.05


def get_default_locations(
        default_locations_filename: str | None = None
) -> dict:
    """
    Returns data from the `default_locations.toml` file.
    """
    if default_locations_filename is None:
        default_locations_path = "~/Code/chemistry/cfour/xsim/config/default_locations.toml"
    else:
        default_locations_path = default_locations_filename
    default_locations_path = os.path.expanduser(default_locations_path)
    with open(default_locations_path, 'rb') as default_locations_file:
        default_locations = tomllib.load(default_locations_file)

    # Make sure that shelltalk (i.e. '~') will work
    for key, value in default_locations.items():
        default_locations[key] = os.path.expanduser(value)

    return default_locations


def assign_xsim_idx_to_states(states, idxs):
    for state in states:
        # states that are not included in the simulation get -1
        if state['ids']['#'] not in idxs:
            state['xsim #'] = -1
            state['ids']['xsim #'] = -1
            continue

        # active states get numbers 1, 2, ...
        for xsim_idx, state_idx in enumerate(idxs):
            if state['ids']['#'] != state_idx:
                continue
            state['xsim #'] = xsim_idx + 1
            state['ids']['xsim #'] = xsim_idx + 1
            break


def add_Mulliken_names_to_modes(modes: list[dict], new_names: list[dict]):
    for new_name in new_names:
        new_freq = new_name['frequency, cm-1']
        found = False
        for mode in modes:
            old_freq = mode['frequency, cm-1']
            if not energies_match(new_freq, old_freq, tolerance=FREQ_TOL):
                continue
            if found is True:
                print(
                    "Warning: Detected multiple modes with frequency "
                    f"{new_freq:.0f}. The same Mulliken name will be assigned "
                    " to multiple modes. You do not want this.",
                    file=sys.stderr
                )
            found = True
            mode['Mulliken'] = new_name['Mulliken']

        if found is True:
            continue
        print(f"Warning: Mode with a Mulliken name with freq {new_freq:.0f} "
              "is absent in the normal modes list.", file=sys.stderr)

    for mode in modes:
        if 'Mulliken' in mode:
            continue
        print(
            f"Warning: Mode with freq {mode['frequency, cm-1']:.0f} "
            "is missing a Mulliken name.", file=sys.stderr
        )


def assign_xsim_idx_to_modes(modes, idxs):
    """
    Modes is a list of normal modes
    idxs is a list of modes that will be used in the simulation

    Add to every mode (which is a dictionary) a new key "xsim #".
    Modes which do not appear in the simulation get -1.
    """

    for mode_no, mode in enumerate(modes):
        if mode_no not in idxs:
            mode['xsim #'] = -1
            continue

        for xsim_idx, mode_idx in enumerate(idxs):
            if mode_idx != mode_no:
                continue
            mode['xsim #'] = xsim_idx + 1
            break


def match_eom_state(unknown, known_states):
    """
    See if unknown state matches any of the states listed in known_states.

    Retruns None if no match is found.

    Example of unknown:
        {"#": 0,
         "irrep": {"#": 8},
         "eom model": "EOMEE-CCSD",
         "energy": {
             "total": {
                 "au": -263.72252668779663
                 }
             }
         }

    """

    for state in known_states:
        if states_are_different(unknown, state):
            continue

        # states match their irrep and total energy -- I say:
        #   > these are the same states
        return state

    print("Error in match_eom_state:\n"
          "An EOM state does not match any known state.", file=sys.stderr)
    print(f"Unknown: {str_eom_state(unknown)}", file=sys.stderr)
    print("Known states:", file=sys.stderr)
    for state in known_states:
        print(str_eom_state(state, full=True), file=sys.stderr)

    return None


def assign_state_xsim_ids_to_gradient_eom_states(gradients, known_states):
    for gradient in gradients:
        # match eom state
        for test_state in gradient['EOM states']:
            matched_state = match_eom_state(test_state, known_states)

            if matched_state is None:
                print("Warning! "
                      "Unable to assign xsim # to a gradient EOM state.")
                print("Setting the state as inactive!")
                print(f"{gradient.keys()=}")
                print("\n\n")
                test_state['xsim #'] = -1
                test_state['ids'] = {'xsim #': -1}
                continue

            # old version
            test_state['xsim #'] = matched_state['xsim #']
            # new version
            test_state['ids'] = matched_state['ids']
            test_state['irrep'] = matched_state['irrep']

            if 'transition' in matched_state['energy'].keys():
                test_state['energy']['transition'] = \
                    matched_state['energy']['transition']


def assign_state_xsim_ids_to_k2_eom_states(kappas_quadratic, known_states):
    """
    Kappas have very baisc structure.
    This function adds to their 'state' dictionary the 'ids' section and
    populates it with the {'#': id} value. This is the EOM # that is stored
    in the EOM states database. The assignment is made based on the eom
    irrep name and irrep energy order.
    """
    for kappas2nd in kappas_quadratic:
        test_name = kappas2nd['state']['text name']
        try:
            test_irrep_no = int(test_name[0])
            test_name = test_name[1:]
        except ValueError:
            test_irrep_no = 1

        for known_state in known_states:
            known_name = known_state['irrep']['name']
            known_irrep_no = known_state['irrep']['energy #']
            if known_name != test_name or known_irrep_no != test_irrep_no:
                continue
            kappas2nd['state']['ids'] = known_state['ids']
            break

        if 'ids' not in kappas2nd['state']:
            print("Warning: unable to assing ids to the kappas state"
                  f" {kappas2nd['state']}",
                  file=sys.stderr)


def assign_state_ids_in_quadratic_kappas_diabatic(
    quad_kapps_diab,
    known_states
):
    """ This function adds ids to each item of the `EOM states` list.
    An example of the quad_kapps_diab json
    ```json
    {
        "EOM states": [
            {"energy": {"total": {"au": -263.75449188023947}}},
            {"energy": {"total": {"au": -263.75449188023947}}}
        ],
        "normal modes": [ 
            {"Mulliken": {"number": 20, "symmetry": "b3g"}},
            {"Mulliken": {"number": 20, "symmetry": "b3g"}}
        ],
        "kappa, cm-1": 1491.7
    },
    {
    ```
    """
    for kappa in quad_kapps_diab:
        for eom_state in kappa['EOM states']:
            new_energy = eom_state['energy']['total']['au']
            for known_state in known_states:
                known_energy = known_state['energy']['total']['au']
                if not energies_match(new_energy, known_energy, tolerance=1e-5):
                    continue
                eom_state['ids'] = known_state['ids']
                break

            if 'ids' not in eom_state:
                print("Warning: unable to assing ids the states "
                      f"(E tot = {new_energy:.5f} Ha) in quadratic diabatic"
                      " kappas.",
                      file=sys.stderr)
                eom_state['ids'] = None


def assign_mode_xsim_ids_to_gradient_modes(gradients, known_nmodes):
    """
    Modes are called the same if they have matching frequency.

    Adds mode's symmetry and Mulliken's name (if available).
    """

    for gradient in gradients:
        for gradient_mode in gradient['gradient']:
            gradient_mode_w = gradient_mode['frequency, cm-1']
            for db_mode in known_nmodes:
                db_mode_w = db_mode['frequency, cm-1']
                if not energies_match(gradient_mode_w, db_mode_w, FREQ_TOL):
                    continue
                # alright these modes are the same
                gradient_mode['xsim #'] = db_mode['xsim #']
                gradient_mode['symmetry'] = db_mode['symmetry']
                if 'Mulliken' in db_mode:
                    gradient_mode['Mulliken'] = db_mode['Mulliken']
                break

            if 'xsim #' not in gradient_mode.keys():
                print("Warning!\n"
                      "A mode from linear kappa (ɷ = "
                      f"{gradient_mode_w:.2f} cm-1) "
                      "doesn't match any known mode.", file=sys.stderr)
                gradient_mode['xsim #'] = None
    return



def modes_match_Mullikens(left, right):
    if left['Mulliken']['symmetry'] != right['Mulliken']['symmetry']:
        return False
    if left['Mulliken']['number'] != right['Mulliken']['number']:
        return False
    return True


def assign_mode_xsim_ids_in_quadratic_kappas_diabatic(
    quad_kapps_diab,
    known_nmodes
):
    """ This function adds ids to each item of the `normal modes` list.
    An example of the quad_kapps_diab json
    ```json
    {
        "EOM states": [
            {"energy": {"total": {"au": -263.75449188023947}}, "ids": ...},
            {"energy": {"total": {"au": -263.75449188023947}}, "ids": ...}
        ],
        "normal modes": [ 
            {"Mulliken": {"number": 20, "symmetry": "b3g"}},
            {"Mulliken": {"number": 20, "symmetry": "b3g"}}
        ],
        "kappa, cm-1": 1491.7
    },
    {
    ```
    """
    for kappa in quad_kapps_diab:
        for mode in kappa['normal modes']:
            for known_mode in known_nmodes:
                if not modes_match_Mullikens(mode, known_mode):
                    continue
                mode['xsim #'] = known_mode['xsim #']
                break
            if 'xsim #' not in mode:
                print(
                    "Warning! Unable to assing xsim # to a mode from"
                    " quadratic diabatic kappas.",
                    file=sys.stderr
                )
                mode['xsim #'] = None
    return
                

def assign_mode_xsim_ids_to_k2_freqs(kappas_quadratic, known_nmodes):
    """
    For each set of kappas add a dictionary item that translates between
    the matrix index and the xsim id.

    Work also for other data with the quadratic-kappas-like structure, i.e.,
    quadratic diabatic shifts.
    """

    # first six are translation/rotations so they get -1, all that follow
    # are the frequency ordered normal modes of the reference
    matrix_pos2xsim_id = {}
    matrix_pos2symmetry = {}
    matrix_pos2Mulliken = {}
    for i in range(6):
        matrix_pos2xsim_id[i] = -1
        matrix_pos2symmetry[i] = 'xxx'
        matrix_pos2Mulliken[i] = None
    normal_modes = sorted(known_nmodes, key=lambda x: x['frequency, cm-1'])
    for pos, db_mode in enumerate(normal_modes):
        matrix_pos2xsim_id[pos + 6] = db_mode['xsim #']
        matrix_pos2symmetry[pos + 6] = db_mode['symmetry']
        if 'Mulliken' in db_mode:
            matrix_pos2Mulliken[pos + 6] = db_mode['Mulliken']
        else:
            matrix_pos2Mulliken[pos + 6] = None

    for k2 in kappas_quadratic:
        k2['pos2xsim#'] = matrix_pos2xsim_id
        k2['pos2symm'] = matrix_pos2symmetry
        k2['pos2Mulliken'] = matrix_pos2Mulliken


def nickname_and_state_match(nickname: dict, state: dict) -> bool:
    """
    Helper for `inject_state_nicknames`.
    """
    if state['irrep']['name'] != nickname['irrep']['name']:
        return False

    if state['irrep']['energy #'] != nickname['irrep']['energy #']:
        return False

    return True


def assure_every_state_has_exactly_one_nickname(states: dict, nicknames: dict):
    # Print warnings for states with zero or more than one nicknames
    for state in states:
        state['n_nicknames'] = 0
        for nickname in nicknames:
            if nickname_and_state_match(nickname, state):
                state['n_nicknames'] += 1

        if state['n_nicknames'] == 0:
            print(f"Warning: No nickname for {str_eom_state(state)}",
                  file=sys.stderr)
        elif state['n_nicknames'] > 1:
            print(f"Error: Multiple nicknames for {str_eom_state(state)}",
                  file=sys.stderr)
            sys.exit(1)

        del state['n_nicknames']


def assure_every_nickname_finds_its_state(states: dict, nicknames: dict):
    """
    Print warnings if there are nicknames that correspond to zero or more than
    one states.
    """
    for nickname in nicknames:
        nickname['matches'] = []
        for state in states:
            if not nickname_and_state_match(nickname, state):
                continue
            nickname['matches'] += [state]

        if len(nickname['matches']) == 0:
            print(f"Warning: The nickname:\n\t{nickname}\n"
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

        del nickname['matches']


def inject_state_nicknames(states: dict, nicknames: dict,
                           match_all_nicknames: bool = True):
    """
    Adds a "nickname" keyword to the "ids" dictionary for each eom state.
    """
    assure_every_state_has_exactly_one_nickname(states, nicknames)
    if match_all_nicknames is True:
        assure_every_nickname_finds_its_state(states, nicknames)

    for nickname in nicknames:
        for state in states:
            if not nickname_and_state_match(nickname, state):
                continue
            state['ids']['nickname'] = nickname['nickname']


def inject_state_nicknames_to_lambdas(lambdas: dict, nicknames: dict):
    """
    Adds a "nickname" keyword to the "ids" dictionary for each eom state.
    """
    for lmbda in lambdas:
        states = lmbda["EOM states"]
        inject_state_nicknames(states, nicknames, match_all_nicknames=False)


def test_state_matched_better_state(state, ccsdt):
    """
    States match when "irrep"'s "name" and "energy #" match.
    """

    if state['irrep']['name'] != ccsdt['irrep']['name']:
        return False
    if state['irrep']['energy #'] != ccsdt['irrep']['energy #']:
        return False

    return True


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


def get_data_with_xsim_ids(
        default_locations_filename: str | None = None
) -> dict[str, Any]:
    """
    Returns a dictionary.

    Dabatbase data with additional information about its index in the xsim's
    simulation. Eg. for EOM states the new data is the 'xsim #'.
    These new ids will be set to -1 for modes or states that are incative in
    the simulation or unrecognized.
    """
    default_locations = get_default_locations(default_locations_filename)

    eom_states = get_eom_states(default_locations)
    state_idxs, mode_idxs = get_active_states_and_modes(default_locations)

    out_pack = {}
    assign_xsim_idx_to_states(eom_states, state_idxs)
    eom_states.sort(key=lambda x: x['ids']['xsim #'])
    out_pack['eom states'] = eom_states

    better_energies_dict = get_better_energies(default_locations)
    out_pack.update(better_energies_dict)

    nmodes = get_normal_modes(default_locations)
    assign_xsim_idx_to_modes(nmodes, mode_idxs)
    nmodes.sort(key=lambda x: x['xsim #'])
    out_pack['nmodes'] = nmodes

    mode_names_Mulliken = get_Mulliken_mode_names(default_locations)
    add_Mulliken_names_to_modes(nmodes, mode_names_Mulliken)

    # kappas, i.e., gradients
    linear_kappas = get_linear_kappas(default_locations)
    assign_state_xsim_ids_to_gradient_eom_states(linear_kappas, eom_states)
    assign_mode_xsim_ids_to_gradient_modes(linear_kappas, nmodes)
    linear_kappas.sort(key=lambda x: x['EOM states'][0]['ids']['xsim #'])
    out_pack['linear kappas'] = linear_kappas

    # Second order kappas
    kappas_quadratic = get_kappas_quadratic(default_locations)
    if kappas_quadratic is None:
        out_pack['quadratic kappas'] = None
    else:
        assign_state_xsim_ids_to_k2_eom_states(kappas_quadratic, eom_states)
        assign_mode_xsim_ids_to_k2_freqs(kappas_quadratic, nmodes)
        kappas_quadratic.sort(key=lambda x: x['state']['ids']['xsim #'])
        out_pack['quadratic kappas'] = kappas_quadratic

    # diabatic shifts second order
    quadratic_kappas_diabatic = get_quadratic_kappas_diabatic(
        default_locations
    )
    if quadratic_kappas_diabatic is None:
        out_pack['quadratic kappas diabatic'] = None
    else:
        assign_state_ids_in_quadratic_kappas_diabatic(
            quadratic_kappas_diabatic, 
            eom_states
        )
        assign_mode_xsim_ids_in_quadratic_kappas_diabatic(
            quadratic_kappas_diabatic,
            nmodes,
        )
        out_pack['quadratic kappas diabatic'] = quadratic_kappas_diabatic

    # lambdas are basically a transition gradient
    lambdas = get_lambdas(default_locations)
    assign_state_xsim_ids_to_gradient_eom_states(lambdas, eom_states)
    assign_mode_xsim_ids_to_gradient_modes(lambdas, nmodes)
    out_pack['lambdas'] = lambdas

    if 'state_nicknames' in default_locations:
        print("Info: Using state nicknames.", file=sys.stderr)
        nicknames = get_state_nicknames(default_locations)
        nicknames = nicknames['state nicknames']
        inject_state_nicknames(eom_states, nicknames)
        inject_state_nicknames_to_lambdas(lambdas, nicknames)

    if 'better energies' in out_pack:
        print("Info: Using better energies", file=sys.stderr)
        inject_better_energies(out_pack)

    return out_pack
