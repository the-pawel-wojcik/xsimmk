import argparse
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from xsim.xsim_ids_processor import get_data_with_xsim_ids
from cfour_parser.text import str_eom_state
import prettytable
from prettytable import TableStyle


def get_args():
    parser = argparse.ArgumentParser()
    which_gradient = parser.add_mutually_exclusive_group()
    which_gradient.add_argument(
        '--kappas',
        default=False,
        action='store_true',
        help='Regular energy gradients.',
    )
    which_gradient.add_argument(
        '--lambdas',
        default=False,
        action='store_true',
        help='"Inter-state" gradients, a.k.a., linear diabatic couplings.',
    )
    parser.add_argument(
        '--save_figure',
        default=False,
        action='store_true',
        help='TODO: not implemented',
    )
    parser.add_argument(
        '--show_frequencies',
        default=False,
        action='store_true',
        help='TODO: not implemented',
    )
    parser.add_argument(
        '--no_pictures',
        default=False,
        action='store_true',
        help='Don\'t show any images. So far the only implemented option',
    )
    parser.add_argument(
        '--quiet',
        help='Don\'t print summaries.',
        default=False,
        action='store_true',
    )
    parser.add_argument(
        '--take_abs',
        help='Broken! Display only the absolute values of the couplings',
        default=False,
        action='store_true'
    )
    args = parser.parse_args()
    return args


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


def inactive_state_involved(gradient):
    for state in gradient['EOM states']:
        if state['xsim #'] <= 0:
            return True

    return False


def show_text_gradient_summary(gradients, name: str):
    sort_gradients(gradients)
    for gradient in gradients:
        if inactive_state_involved(gradient):
            continue

        print("Gradient calculated for the following state(s):")
        for state in gradient['EOM states']:
            print(str_eom_state(state))

        active_grad = [
            mode for mode in gradient['gradient'] if mode['xsim #'] > 0
        ]
        inactive_grad = [
            mode for mode in gradient['gradient'] if mode['xsim #'] <= 0
        ]

        print(f"{name[0].upper() + name[1:]}s included in the model:")
        print_active_gradient_table(active_grad)

        print(f"{name[0].upper() + name[1:]}s ommited in the model:")
        print_inactive_gradient_table(inactive_grad)


def format_Mulliken(_, v):
    return f"{v['number']:2d}({v['symmetry']})"


def print_active_gradient_table(gradient: list):
    if len(gradient) == 0:
        print("No active modes in the model.")
        return

    table = prettytable.PrettyTable()
    table.field_names = [key for key in gradient[0].keys()]
    for row in sorted(gradient, key=lambda x: x['Mulliken']['number']):
        table.add_row([val for val in row.values()])
    table.set_style(TableStyle.SINGLE_BORDER)

    for key in ['frequency, cm-1', 'gradient, cm-1']:
        table.custom_format[key] = lambda _, v: f"{v:.2f}"
        table.align[key] = 'r'

    table.custom_format['Mulliken'] = format_Mulliken

    table.custom_format['xsim #'] = lambda _, v: f"{v:-3d}"
    table.align['xsim #'] = 'r'
    print(table)


def print_inactive_gradient_table(gradient):
    """ Differs slightly from `print_active_gradient_table` as this one does
    not print the `xsim #` column.
    """
    if len(gradient) == 0:
        print("No inactive modes in the model.")
        return

    table = prettytable.PrettyTable()
    table.field_names = [
        key for key in list(gradient[0].keys()) if key != 'xsim #'
    ]
    for row in sorted(gradient, key=lambda x: x['Mulliken']['number']):
        table.add_row([val for val in list(row.values()) if val != -1])
    table.set_style(TableStyle.SINGLE_BORDER)

    for key in ['frequency, cm-1', 'gradient, cm-1']:
        table.custom_format[key] = lambda _, v: f"{v:.0f}"
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


def collect_sns_matrices(
    gradients: list,
    normal_modes: list,
    take_abs: bool = False
):
    gradients_matrix = [[mode['frequency, cm-1']] for mode in normal_modes]
    labels = [r'freq']

    for grad in gradients:
        if inactive_state_involved(grad):
            continue
        # # HACK: specific for Pyrazine shows only couplings to 1B3u
        # if not coupling_to_B3u(lmbda):
        #     continue

        label = ""
        for state in grad['EOM states']:

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

        # Collect componets
        matrix = [[mode['gradient, cm-1']] for mode in grad['gradient']]
        gradients_matrix = [current_row + new_row
                            for current_row, new_row
                            in zip(gradients_matrix, matrix)]

    n_couplings = len(gradients_matrix[0]) - 1
    annotations_matrix = collect_sns_annotations(normal_modes, n_couplings)

    if take_abs is True:
        gradients_matrix = [
            # [row[0]] + [val * 10 for val in row[1:]]  #  HACK: for phenoxide
            [row[0]] + [val for val in row[1:]]
            for row in gradients_matrix
        ]

    return gradients_matrix, annotations_matrix, labels


def symmetry_match(left: dict, right: dict) -> bool:
    return left['Mulliken']['symmetry'] == right['Mulliken']['symmetry']


def prepare_yticklabels(
        normal_modes: list[dict],
        mulliken_to_idx: list[dict],
        use_Mulliken: bool = True,
) -> list[str] | None:

    symmetry_labels = all(mode['symmetry'] != '???' for mode in normal_modes)
    symmetry_labels = False  # Hack
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
                    lbl = ""
                    # lbl = str(next_mode['Mulliken']['number'])
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
    else:
        return None

    return mode_symmetries


def sort_gradients(gradients: list[dict]):
    for grad in gradients:
        grad['EOM states'].sort(key=lambda x: x['energy']['transition']['eV'])

    gradients.sort(key=lambda x: sum(
        eom['energy']['transition']['eV'] for eom in x['EOM states']
    ))


def show_sns_gradient_summary(
        ax: Axes,
        normal_modes,
        gradients,
        take_abs: bool = True,
        show_frequencies: bool = True,
        **heatmap_kwargs
):
    use_Mulliken = all(
        all(
            'Mulliken' in grd for grd in gradient['gradient']
        ) for gradient in gradients
    )

    # Show gradients in the increasing energy order
    sort_gradients(gradients)

    normal_modes.sort(key=lambda x: x['frequency, cm-1'])

    couplings, annotations, labels = collect_sns_matrices(
        gradients=gradients,
        normal_modes=normal_modes,
        take_abs=take_abs,
    )

    if show_frequencies is False:
        couplings = [row[1:] for row in couplings]
        annotations = [row[1:] for row in annotations]
        labels = labels[1:]

    max_row = max(couplings, key=lambda x: abs(max(x, key=lambda y: abs(y))))
    max_value = abs(max(max_row, key=lambda x: abs(x)))

    # defaults are None
    parameters = {
        "center": 0.0,
        "vmin": -max_value,
        # "vmin": 0.0,
        "vmax": max_value,
        # "cmap": "hot_r",
        # "cmap": "copper_r",
        "cmap": "vlag",  # diverging (top choice for diverging)
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

    if mode_symmetries is not None:
        parameters['yticklabels'] = mode_symmetries

    parameters.update(heatmap_kwargs)

    # Hack: show 1-based mode names
    yticks_one_based = list()
    for idx in range(1, len(couplings)+1):
        if idx % 2 == 0:
            yticks_one_based.append('')
        else:
            yticks_one_based.append(str(idx))
    parameters['yticklabels'] = yticks_one_based

    sns.heatmap(couplings, **parameters)

    # Hack: show 1-based mode names
    ax.tick_params(axis='y', labelrotation=0)

    return ax


def main():
    args = get_args()
    data = get_data_with_xsim_ids()
    if args.kappas is True:
        gradients = data['linear kappas']
        name = 'kappas'
    elif args.lambdas is True:
        gradients = data['lambdas']
        name = 'lambdas'
    else:
        print(
            "Error. Pick gradient, i.e., --kappas or --lambdas",
            file=sys.stderr
        )
        sys.exit(1)

    nmodes = data['nmodes']

    quiet = args.quiet
    if not quiet:
        show_text_gradient_summary(gradients, name)

    no_pictures = args.no_pictures
    if not no_pictures:
        fig, ax = plt.subplots(layout='constrained')
        show_frequencies = args.show_frequencies
        take_abs = args.take_abs
        show_sns_gradient_summary(
            ax,
            normal_modes=nmodes,
            gradients=gradients,
            take_abs=take_abs,
            show_frequencies=show_frequencies,
        )

        save = args.save_figure
        if save is False:
            plt.show()
        else:
            fname='gradients.pdf'
            fig.savefig(fname=fname)
            print(f"Info: Figure saved as {fname}", file=sys.stderr)


if __name__ == "__main__":
    main()
