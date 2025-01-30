import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from cfour_parser.text import str_eom_state
from xsim.db.prepare import energies_match
from xsim.xsim_ids_processor import get_data_with_xsim_ids
from xsim.helpers.plot_diabatic_couplings import show_sns_lambdas_summary
from templates.save import save_to_file
import sys
from adjustText import adjust_text


def max_energy(states):
    """
    the max value that I consider is the smallest of the max energies in each
    irrep
    """
    max_energies = dict()
    for state in states:
        irrep = state['irrep']['#']
        energy = state['energy']['transition']['eV']

        if irrep not in max_energies:
            max_energies[irrep] = energy
            continue

        if energy > max_energies[irrep]:
            max_energies[irrep] = energy

    max_energy = min(max_energies.values())
    return max_energy


def state_plot(ax: mpl.axes.Axes, state: dict) -> mpl.text.Text:
    if 'transition' not in state['energy'].keys():
        print("Warning: State (" + str_eom_state(state) + ") is missing "
              "transition energy.", file=sys.stderr)
        return

    if 'nickname' in state['ids']:
        text = state['ids']['nickname']
    elif 'name' in state['irrep'].keys():
        name = state['irrep']['name']
        text = name[0] + r"$ _{" + name[1:] + r"}$"
    else:
        text = '???'

    textkwargs = {
        'verticalalignment': 'center_baseline',
        's': text
    }

    energy = state['energy']['transition']['eV']

    # HACK: Pyrazine
    if energy > 8.0:
        return
    ys = [energy] * 2

    kwargs = {
        'linewidth': 3
    }
    if 'xsim #' in state.keys() and state['xsim #'] > 0:
        kwargs['color'] = 'tab:blue'
    else:
        kwargs['color'] = 'tab:gray'

    tdm = state['Norm of oscillator strength']

    if energies_match(tdm, 0.0, 0.0001):
        # Dark state
        xs = [0.55, 0.95]
        kwargs['linestyle'] = ':'
        # kwargs['color'] = 'grey'
        text_obj = ax.text(xs[1], ys[0], ha='left', **textkwargs)
    else:
        # Bright state
        xs = [0.05, 0.45]
        kwargs['linestyle'] = '-'
        # kwargs['color'] = 'grey'
        text_obj = ax.text(xs[0] - 0.05, ys[0], ha='right', **textkwargs)

    ax.plot(xs, ys, **kwargs)
    return text_obj


def mk_filename(lambdas):
    coupled_states = {}
    for lmbda in lambdas:
        for lambda_state in lmbda['EOM states']:
            if lambda_state['xsim #'] <= 0:
                continue
            lambda_irrep_name = lambda_state['irrep']['name']
            if lambda_irrep_name in coupled_states.keys():
                continue
            state_energy_eV = lambda_state['energy']['transition']['eV']
            coupled_states[lambda_irrep_name] = state_energy_eV

    coupled_states = [
        {'name': key, 'energy': val} for key, val in coupled_states.items()
    ]
    coupled_states.sort(key=lambda x: x['energy'])

    file_name = ''
    for state in coupled_states:
        lambda_irrep = state['name']
        file_name += lambda_irrep + '+'

    file_name = file_name[:-1]
    return file_name


def decongest_state_lables(
        state_labels: list[mpl.text.Text],
        ax: mpl.axes.Axes,
):
    only_move = {
        "text": "y",
        "static": "y",
        "explode": "y",
        "pull": "y",
    }
    adjust_text(
        texts=state_labels, ax=ax,
        avoid_self=False,
        only_move=only_move,
    )


def visualize_the_couplings(
        eom_states,
        normal_modes,
        lambdas,
        show_frequencies: bool = True,
        save: bool = False
):
    height = 3.75
    width = 16/9 * height
    fig = plt.figure(figsize=(width, height), layout='constrained')
    gs0: mpl.gridspec.GridSpec = fig.add_gridspec(
        nrows=1,
        ncols=3,
        width_ratios=[1, 1, 0.25],
    )
    gs1 = gs0[0].subgridspec(1, 1)
    gs2 = gs0[1].subgridspec(1, 1)
    gs3 = gs0[2].subgridspec(1, 1)
    ax = fig.add_subplot(gs1[0])
    ax_nmmodes = fig.add_subplot(gs2[0])
    ax_cbar = fig.add_subplot(gs3[0])

    states = [state for state in eom_states if 'transition' in state['energy']]

    state_labels = list()
    for state in states:
        state_label = state_plot(ax, state)
        if state_label is not None:
            state_labels.append(state_label)

    heatmap_kwargs = {
        'cbar': True,
        'cbar_ax': ax_cbar,
    }
    show_sns_lambdas_summary(
        ax=ax_nmmodes,
        lambdas=lambdas,
        normal_modes=normal_modes,
        take_abs=True,
        show_frequencies=show_frequencies,
        **heatmap_kwargs
    )

    ax_cbar.set_xlabel(r'cm$^{-1}$')
    ax.set_ylabel("Vertical energy/eV")
    ax.minorticks_on()
    ax.set_xlim((-0.5, 1.5))
    ax.set_xticks([])
    ax.spines[['right', 'top', 'bottom']].set_visible(False)

    decongest_state_lables(state_labels, ax)

    if save is False:
        plt.show()
    else:
        filename = mk_filename(lambdas)
        save_to_file(fig, outname=filename, filetype='svg')
        print(f"Info: Figure saved as {filename}.svg", file=sys.stderr)


def main():
    data = get_data_with_xsim_ids()
    eom_states = data['eom states']
    lambdas = data['lambdas']
    nmodes = data['nmodes']
    # visualize_the_couplings(eom_states, nmodes, lambdas, save=True)
    visualize_the_couplings(eom_states, nmodes, lambdas, save=False)


if __name__ == "__main__":
    main()
