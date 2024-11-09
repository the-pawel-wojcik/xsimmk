import sys
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from xsim.xsim_ids_processor import get_data_with_xsim_ids
from cfour_parser.text import str_eom_state
import prettytable


def pprint_Mulliken_irrep(mode: dict) -> str:
    mulliken = mode['Mulliken']
    # freq = mode['frequency, cm-1']  # U-turn, I don't want to text this info
    irrep = mulliken['symmetry']
    number = mulliken['number']
    pretty = irrep[0] + r"$ _{" + irrep[1:] + r"}$"
    # output = f"{pretty:<3s}" + f"{number:3d}" + f" {freq:4.0f}"
    output = f"{pretty:<3s}" + f"{number:3d}"
    return output


def pretty_print_irrep(irrep):
    pretty = irrep[0] + r"$ _{" + irrep[1:] + r"}$"
    return pretty


def coupling_to_inactive_state(lmbda):
    if lmbda['EOM states'][0]['xsim #'] <= 0:
        return True

    if lmbda['EOM states'][1]['xsim #'] <= 0:
        return True

    return False


def show_text_lambdas_summary(lambdas):
    for lmbda in lambdas:
        if coupling_to_inactive_state(lmbda):
            continue

        print("Coupling between:")
        for state in lmbda['EOM states']:
            print(str_eom_state(state))

        active_lmbda = [
            mode for mode in lmbda['gradient'] if mode['xsim #'] > 0
        ]
        inactive_lmbda = [
            mode for mode in lmbda['gradient'] if mode['xsim #'] <= 0
        ]

        print("Couplings included in the model:")
        print_active_lambdas_table(active_lmbda)

        print("Couplings ommited in the model:")
        print_inactive_lambdas_table(inactive_lmbda)


def format_Mulliken(f, v):
    return f"{v['number']:2d}({v['symmetry']})"


def print_active_lambdas_table(lmbda):
    if len(lmbda) == 0:
        print("No active modes in the model.")
        return

    table = prettytable.PrettyTable()
    table.field_names = [key for key in lmbda[0].keys()]
    for row in lmbda:
        table.add_row([val for val in row.values()])
    table.set_style(prettytable.SINGLE_BORDER)

    for key in ['frequency, cm-1', 'gradient, cm-1']:
        table.custom_format[key] = lambda f, v: f"{v:.2f}"
        table.align[key] = 'r'

    table.custom_format['Mulliken'] = format_Mulliken

    table.custom_format['xsim #'] = lambda f, v: f"{v:-3d}"
    table.align[key] = 'r'
    print(table)


def print_inactive_lambdas_table(lmbda):
    if len(lmbda) == 0:
        print("No inactive modes in the model.")
        return

    table = prettytable.PrettyTable()
    table.field_names = [key for key in list(
        lmbda[0].keys()) if key != 'xsim #']
    for row in lmbda:
        table.add_row([val for val in list(row.values()) if val != -1])
    table.set_style(prettytable.SINGLE_BORDER)

    for key in ['frequency, cm-1', 'gradient, cm-1']:
        table.custom_format[key] = lambda f, v: f"{v:.2f}"
        table.align[key] = 'r'

    table.custom_format['Mulliken'] = format_Mulliken

    print(table)


def collect_sns_annotations(normal_modes, n_couplings: int):
    """
    Creates a matrix with annotations indicating which mode is active in the
    simulation.
    """

    # Annotate with active or not
    annotate = []
    for mode in normal_modes:
        if mode['xsim #'] is None:
            print("Warning! Unable to tell if mode ("
                  f"{mode['frequency, cm-1']:.2f}"
                  " cm-1) is active as it misses the xsim #.", file=sys.stderr)
            annotate += [["?"] * (n_couplings + 1)]
            continue

        active = mode['xsim #'] > 0
        if active:
            note = "x"
        else:
            note = ""
        annotate += [[note] * (n_couplings + 1)]

    return annotate


def coupling_to_B3u(lmbda) -> bool:
    """Customization for pyrazine"""

    if lmbda['EOM states'][0]['ids']['nickname'] == r'B$_{3u}$':
        return True

    if lmbda['EOM states'][1]['ids']['nickname'] == r'B$_{3u}$':
        return True

    return False


def collect_sns_matrices(lambdas: list, normal_modes: list):
    couplings_matrix = [[mode['frequency, cm-1']] for mode in normal_modes]
    labels = [r'freq']

    for lmbda in lambdas:
        if coupling_to_inactive_state(lmbda):
            continue
        # HACK: specific for Pyrazine shows only couplings to 1B3u
        if not coupling_to_B3u(lmbda):
            continue

        label = ""
        for state in lmbda['EOM states']:

            if 'nickname' in state['ids']:
                label += ' ' + state['ids']['nickname']

            elif 'energy #' in state['irrep'].keys()\
                    and 'name' in state['irrep'].keys():
                state_in_irrep_number = state['irrep']['energy #']
                label += ' '
                if state_in_irrep_number > 1:
                    label += str(state_in_irrep_number)
                label += pretty_print_irrep(state['irrep']['name'])

            elif 'name' in state['irrep'].keys():
                label += ' ' + pretty_print_irrep(state['irrep']['name'])

            else:
                label += ' irrep #' + state['irrep']['#']

        labels += [label]

        # Collect couplings
        matrix = [[mode['gradient, cm-1']] for mode in lmbda['gradient']]
        couplings_matrix = [current_row + new_row
                            for current_row, new_row
                            in zip(couplings_matrix, matrix)]

    n_couplings = len(couplings_matrix[0]) - 1
    annotations_matrix = collect_sns_annotations(normal_modes, n_couplings)

    return couplings_matrix, annotations_matrix, labels


def symmetry_match(left: dict, right: dict) -> bool:
    return left['Mulliken']['symmetry'] == right['Mulliken']['symmetry']


def prepare_yticklabels(
        normal_modes: list[dict],
        mulliken_to_idx: list[dict],
        use_Mulliken: bool = True,
) -> list[str]:

    symmetry_labels = False
    symmetry_labels = True
    if symmetry_labels is True:
        if use_Mulliken:
            pre_mode_symmetries = [
                {
                    'Mulliken': mode['Mulliken'],
                    'frequency, cm-1': mode['frequency, cm-1']
                }
                for mode in normal_modes
            ]

            pre_mode_symmetries = [
                pre_mode_symmetries[idx['idx']] for idx in mulliken_to_idx
            ]

            mode_it = iter(pre_mode_symmetries)
            last_mode = next(mode_it)
            mode_symmetries = [pprint_Mulliken_irrep(last_mode)]
            while True:
                try:
                    next_mode = next(mode_it)
                except StopIteration:
                    break

                if symmetry_match(last_mode, next_mode):
                    lbl = str(next_mode['Mulliken']['number'])
                    # lbl += f" {next_mode['frequency, cm-1']:4.0f}"
                    mode_symmetries.append(lbl)
                else:
                    mode_symmetries.append(pprint_Mulliken_irrep(next_mode))

                last_mode = next_mode

        else:
            mode_symmetries = [
                pretty_print_irrep(mode['symmetry'])
                for mode in normal_modes
            ]

    return mode_symmetries, None


def sort_lambdas(lambdas: list[dict]):
    for lmbda in lambdas:
        lmbda['EOM states'].sort(key=lambda x: x['energy']['transition']['eV'])

    lambdas.sort(key=lambda x: sum(
        eom['energy']['transition']['eV'] for eom in x['EOM states'])
    )


def show_sns_lambdas_summary(ax, normal_modes, lambdas, **heatmap_kwargs):
    use_Mulliken = all(
        all(
            'Mulliken' for grd in lmbda['gradient']
        ) for lmbda in lambdas
    )

    sort_lambdas(lambdas)
    normal_modes.sort(key=lambda x: x['frequency, cm-1'])

    couplings, annotations, labels = collect_sns_matrices(
        lambdas=lambdas,
        normal_modes=normal_modes,
    )

    max_row = max(couplings, key=lambda x: abs(max(x, key=lambda y: abs(y))))
    max_value = abs(max(max_row, key=lambda x: abs(x)))

    # defaults are None
    parameters = {
        "vmin": -max_value,
        "vmax": max_value,
        "center": 0.0,
        # "cmap": "YlGnBu",
        # "cmap": "coolwarm",  # divergin
        "cmap": "vlag",  # divergin
        "fmt": '',  # string formatting options for displaying annot
        "linewidth": 0.5,  # width of the lines that divide cells.
        "square": True,
        "cbar": True,  # draw color bar
        "ax": ax,  # Axes where the figure will be drawn
        "cbar_ax": None,  # Axis for the color bar, i.e., the legend
        "xticklabels": labels,
    }

    mulliken_to_idx = list()
    if use_Mulliken is True:
        for idx, mode in enumerate(normal_modes):
            mulliken_to_idx.append({
                'Mulliken #': mode['Mulliken']['number'],
                'idx': idx,
            })

        mulliken_to_idx.sort(key=lambda x: x['Mulliken #'])
        couplings = [couplings[idx['idx']] for idx in mulliken_to_idx]
        annotations = [annotations[idx['idx']] for idx in mulliken_to_idx]

    parameters.update({
        "annot": annotations,  # bool or matrix
    })

    mode_symmetries = prepare_yticklabels(
        normal_modes, mulliken_to_idx, use_Mulliken
    )

    parameters['yticklabels'] = mode_symmetries

    parameters.update(heatmap_kwargs)
    ax: mpl.axes.Axes = sns.heatmap(couplings, **parameters)

    return ax


def main():
    data = get_data_with_xsim_ids()
    lambdas = data['lambdas']
    nmodes = data['nmodes']

    show_text_lambdas_summary(lambdas)

    fig, ax = plt.subplots(layout='constrained')
    show_sns_lambdas_summary(ax, normal_modes=nmodes, lambdas=lambdas)
    plt.show()


if __name__ == "__main__":
    main()
