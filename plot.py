import argparse
import matplotlib.pyplot as plt
import numpy as np
import h5py
import gusto_dataset


global_vars = dict(num_files=0, file_index=0)


def plot_faces(filename, args):
    dset = gusto_dataset.GustoDataset(filename)
    segments = dset.get_face_segments()
    for seg in segments:
        plt.plot(seg[:,0], seg[:,2], '-o', c='k')
    plt.axis('equal')
    plt.margins(0.1)


def plot_1d(filename, args):
    dset = gusto_dataset.GustoDataset(filename)
    x = dset.get_cell_variable('x1')
    y = dset.get_cell_variable('u1')
    p = dset.get_cell_variable('pg')
    d = dset.get_cell_variable('dg')
    s = np.log(p / d**(4./3))
    #y = s
    c = float(global_vars['file_index']) / global_vars['num_files']
    plt.plot(x, y, '-', c=[c]*3)


def triangle_variable_plot(filename, args):
    dset = gusto_dataset.GustoDataset(filename)
    x = dset.get_cell_variable('x1')
    z = dset.get_cell_variable('x3')
    f = dset.get_cell_variable('dg')
    plt.tripcolor(x, z, f)
    plt.axis('equal')
    plt.colorbar()


def triangle_mesh_plot(filename, args):
    dset = gusto_dataset.GustoDataset(filename)
    x = dset.get_cell_variable('x1')
    z = dset.get_cell_variable('x3')
    plt.triplot(x, z)
    plt.axis('equal')


def mesh_plot(filename, args):
    from matplotlib.collections import PolyCollection
    dset = gusto_dataset.GustoDataset(filename)
    vert = dset.get_cell_polygons()
    data = dset.get_cell_variable(args.data)
    cells = PolyCollection(vert, array=data, cmap=args.cmap, linewidths=0.0)
    fig = plt.figure()
    ax0 = fig.add_subplot(1, 1, 1)
    ax0.add_collection(cells)
    ax0.autoscale_view()
    ax0.set_aspect('equal')
    fig.colorbar(cells)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('command')
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-d', '--data', default='dg')
    parser.add_argument('--cmap', default='jet')

    args = parser.parse_args()

    plots = {'1d': plot_1d,
             'trimesh': triangle_mesh_plot,
             'triplot': triangle_variable_plot,
             'mesh': mesh_plot,
             'faces': plot_faces}

    global_vars['num_files'] = len(args.filenames)
    for n, f in enumerate(args.filenames):
        global_vars['file_index'] = n
        plots[args.command](f, args)
    plt.show()
