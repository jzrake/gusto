import argparse
import matplotlib.pyplot as plt
import numpy as np
import h5py
import gusto_dataset


global_vars = dict(num_files=0, file_index=0)


def plot_faces(filename):
    dset = gusto_dataset.GustoDataset(filename)
    segments = dset.get_face_segments()
    for seg in segments:
        plt.plot(seg[:,0], seg[:,2], '-o', c='k')
    plt.axis('equal')


def plot_1d(filename):
    dset = gusto_dataset.GustoDataset(filename)
    x = dset.get_cell_variable('x1')
    y = dset.get_cell_variable('u1')
    c = float(global_vars['file_index']) / global_vars['num_files']
    plt.plot(x, y, '-o', c=[c]*3)


def triangle_variable_plot(filename):
    dset = gusto_dataset.GustoDataset(filename)
    x = dset.get_cell_variable('x1')
    z = dset.get_cell_variable('x3')
    f = dset.get_cell_variable('dg')
    plt.tripcolor(x, z, f)
    plt.axis('equal')


def triangle_mesh_plot(filename):
    dset = gusto_dataset.GustoDataset(filename)
    x = dset.get_cell_variable('x1')
    z = dset.get_cell_variable('x3')
    plt.triplot(x, z)
    plt.axis('equal')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('command')
    parser.add_argument('filenames', nargs='+')

    args = parser.parse_args()

    plots = {'1d': plot_1d,
             'trimesh': triangle_mesh_plot,
             'triplot': triangle_variable_plot,
             'faces': plot_faces}

    global_vars['num_files'] = len(args.filenames)
    for n, f in enumerate(args.filenames):
        global_vars['file_index'] = n
        plots[args.command](f)
    plt.show()
