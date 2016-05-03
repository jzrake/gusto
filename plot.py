import argparse
import matplotlib.pyplot as plt
import numpy as np
import h5py
import gusto_dataset



def list_from_csv(csv):
    if ',' not in csv: return [csv]
    else: return csv.split(',')


def plot_1d(ax, filename, args):
    dset = gusto_dataset.GustoDataset(filename)
    x = dset.get_cell_variable(args.x_data, row=args.row)

    for key in list_from_csv(args.y_data):
        y = dset.get_cell_variable(key, row=args.row)
        ax.plot(x, y, args.linestyle, label=key)

    ax.set_xscale('log' if args.logx else 'linear')
    ax.set_yscale('log' if args.logy else 'linear')
    ax.set_xlabel(args.x_data)
    ax.legend()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('command')
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-x', '--x-data', default='X')
    parser.add_argument('-y', '--y-data', default='dg')
    parser.add_argument('--cmap', default='jet')
    parser.add_argument('--row', default=0, type=int)
    parser.add_argument('--logx', action='store_true', default=False)
    parser.add_argument('--logy', action='store_true', default=False)
    parser.add_argument('-s', '--linestyle', default='-o')

    args = parser.parse_args()

    plots = {'1d': plot_1d}
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for n, f in enumerate(args.filenames):
        plots[args.command](ax, f, args)

    plt.show()
