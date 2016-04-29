import argparse
import matplotlib.pyplot as plt
import numpy as np
import h5py
import gusto_dataset


global_vars = dict(num_files=0, file_index=0)



def plot_1d(filename, args):
    dset = gusto_dataset.GustoDataset(filename)
    x = dset.get_cell_variable('x3', row=args.row)
    y = dset.get_cell_variable(args.data, row=args.row)
    c = float(global_vars['file_index']) / global_vars['num_files']
    plt.plot(x, y, '-o', c=[c]*3)
    plt.xlabel('r')
    plt.ylabel(args.data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('command')
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-d', '--data', default='dg')
    parser.add_argument('--cmap', default='jet')
    parser.add_argument('--row', default=0, type=int)
    parser.add_argument('--log', action='store_true', default=False)

    args = parser.parse_args()

    plots = {'1d': plot_1d}

    global_vars['num_files'] = len(args.filenames)
    for n, f in enumerate(args.filenames):
        global_vars['file_index'] = n
        plots[args.command](f, args)
    plt.show()
