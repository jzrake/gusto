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
    plt.margins(0.1)


def plot_1d(filename):
    dset = gusto_dataset.GustoDataset(filename)
    x = dset.get_cell_variable('x1')
    y = dset.get_cell_variable('u1')
    p = dset.get_cell_variable('pg')
    d = dset.get_cell_variable('dg')
    s = np.log(p / d**(4./3))
    #y = s
    c = float(global_vars['file_index']) / global_vars['num_files']
    plt.plot(x, y, '-', c=[c]*3)


def triangle_variable_plot(filename):
    dset = gusto_dataset.GustoDataset(filename)
    x = dset.get_cell_variable('x1')
    z = dset.get_cell_variable('x3')
    f = dset.get_cell_variable('b2')
    plt.tripcolor(x, z, f)
    plt.axis('equal')
    plt.colorbar()


def triangle_mesh_plot(filename):
    dset = gusto_dataset.GustoDataset(filename)
    x = dset.get_cell_variable('x1')
    z = dset.get_cell_variable('x3')
    plt.triplot(x, z)
    plt.axis('equal')


def mesh_plot(filename):
    from matplotlib.collections import PolyCollection
    from matplotlib import cm
    cmap = cm.get_cmap('jet')
    scal = cm.ScalarMappable(cmap=cmap)
    dset = gusto_dataset.GustoDataset(filename)
    vert = dset.get_cell_polygons()
    data = dset.get_cell_variable('dg')
    colors =  scal.to_rgba(data)
    cells = PolyCollection(vert, linewidths=0.0, facecolors=colors)
    fig = plt.figure()
    ax0 = fig.add_subplot(1, 1, 1)
    ax0.add_collection(cells)
    ax0.autoscale_view()
    ax0.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('command')
    parser.add_argument('filenames', nargs='+')

    args = parser.parse_args()

    plots = {'1d': plot_1d,
             'trimesh': triangle_mesh_plot,
             'triplot': triangle_variable_plot,
             'mesh': mesh_plot,
             'faces': plot_faces}

    global_vars['num_files'] = len(args.filenames)
    for n, f in enumerate(args.filenames):
        global_vars['file_index'] = n
        plots[args.command](f)
    plt.show()
