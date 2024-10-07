#!/usr/bin/env python3

import argparse
import json
import re
import sys


def show_kappas(matrix):
    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, ax = plt.subplots(figsize=(2, 4))

    parameters = {
        # "vmin": 0.0,
        # "vmax": 3.0,
        "center": 0.0,
        # "cmap": "YlGnBu",
        "cmap": "vlag",  # diverging
        "linewidth": 0.5,  # width of the lines that divide cells.
        "square": True,
        "cbar": True,  # draw color bar
        "ax": ax,  # Axes where the figure will be drawn
    }

    ax = sns.heatmap(matrix, **parameters)
    plt.title("Quadratic kappas")
    plt.show()


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('xquadmodeloutput', help='xquadmodel output.'
                        'Parses quadratic kappas.')
    parser.add_argument('state_name',
                        help='A string that identifies the state.')
    parser.add_argument('-j', '--json', default=False, action='store_true')
    parser.add_argument('-v', '--verbose', default=False, action='store_true')
    args = parser.parse_args()
    return args


def parse_quadmodel_freq(xqmout):
    """
    Parses the output of xquadmodel executable (part of CFOUR).
    It looks for the 'Final state harmonic frequencies' and returns them
    as a list.
    """

    freq_header = "Final state harmonic frequencies (negative is imaginary)"
    freq_terminator = "force constant matrix in reference reduced normal"
    frequencies = []
    while True:
        line = xqmout.readline()

        # EOF
        if line == "":
            break

        if line.strip() == freq_header:
            while True:
                line = xqmout.readline().strip()
                if line == freq_terminator:
                    break
                line = line.split()
                for freq in line:
                    frequencies += [float(freq)]
            break

    return frequencies


def parse_quadmodel_fcm(xqmout, freqs):
    """
    Parses the output of xquadmodel executable (part of CFOUR).
    It looks for the 'Final state harmonic frequencies' and returns them
    as a list.
    """

    kappas_header = "force constant matrix in reference reduced normal"
    kappas_header += " coordinates (cm-1)"
    kappas = [[0.0 for i in range(len(freqs))] for j in range(len(freqs))]

    column_pattern = re.compile(r'(?:\s*COLUMN\s+(\d+))')
    while True:
        line = xqmout.readline()

        # EOF
        if line == "":
            print("Warning! "
                  "End of file reached before parsing of second order kappas.")
            break

        if line.strip() == kappas_header:
            while True:
                xqmout.readline()  # skip one
                line = xqmout.readline()
                match = column_pattern.findall(line)
                if match is None:
                    print("Error in parsing second order kappas!",
                          file=sys.stderr)
                    return kappas
                # Note that xquadmodel prints rows and column numbers starting
                # with 1 while python uses numbering that starts with 0
                cols = [int(col) for col in match]
                for row in range(len(freqs)):
                    line = xqmout.readline().split()
                    if int(line[1]) - 1 != row:
                        print("Error in parsing second order kappas!"
                              " Rows don't match.", file=sys.stderr)
                    for col_no, col in enumerate(cols):
                        kappas[row][col-1] = float(line[2 + col_no])

                if cols[-1] == len(freqs):
                    break
            break

    return kappas


def main():
    args = get_args()
    with open(args.xquadmodeloutput, 'r') as xqmout:
        quadmodel_freq = parse_quadmodel_freq(xqmout)

    with open(args.xquadmodeloutput, 'r') as xqmout:
        quadmodel_fcm = parse_quadmodel_fcm(xqmout, quadmodel_freq)

    if args.verbose is True:
        show_kappas(quadmodel_fcm)

    if args.json is True:
        outpack = {
            'kappa quadratic': quadmodel_fcm,
            'frequencies': quadmodel_freq,
            "state": {'text name': args.state_name, },
        }
        print(json.dumps(outpack))


if __name__ == "__main__":
    main()
