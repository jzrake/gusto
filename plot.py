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
    x = dset.get_cell_variable('x3', row=args.row)
    y = dset.get_cell_variable(args.data, row=args.row)
    c = float(global_vars['file_index']) / global_vars['num_files']
    plt.plot(x, y, '-o', c=[c]*3)
    plt.xlabel('z')
    plt.ylabel(args.data)


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


def triangle_vert_plot(filename, args):
    dset = gusto_dataset.GustoDataset(filename)
    x = dset.get_vert_variable('x1')
    z = dset.get_vert_variable('x3')
    f = dset.get_vert_variable(args.data)
    plt.tripcolor(x, z, f)
    plt.axis('equal')


def mesh_plot(filename, args):
    from matplotlib.collections import PolyCollection
    dset = gusto_dataset.GustoDataset(filename)
    vert = dset.get_cell_polygons()
    data = dset.get_cell_variable(args.data, log=args.log)

    # u0 = dset.get_cell_variable('u0')
    # b0 = dset.get_cell_variable('b0')
    # b1 = dset.get_cell_variable('b1')
    # b2 = dset.get_cell_variable('b2')
    # b3 = dset.get_cell_variable('b3')
    # dg = dset.get_cell_variable('dg')
    # bb = b1*b1 + b2*b2 + b3*b3 - b0*b0
    # data = np.log10(bb / (dg * u0))

    # u1 = dset.get_cell_variable('u1')
    # u3 = dset.get_cell_variable('u3')
    # up = (u1**2 + u3**2)**0.5
    # data = up
    # data = np.log10(data)

    cells = PolyCollection(vert, array=data, cmap=args.cmap,
                           linewidths=0.0, antialiased=True)
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
    parser.add_argument('--row', default=0, type=int)
    parser.add_argument('--log', action='store_true', default=False)

    args = parser.parse_args()

    plots = {'1d': plot_1d,
             'vert': triangle_vert_plot,
             'trimesh': triangle_mesh_plot,
             'triplot': triangle_variable_plot,
             'mesh': mesh_plot,
             'faces': plot_faces}

    global_vars['num_files'] = len(args.filenames)
    for n, f in enumerate(args.filenames):
        global_vars['file_index'] = n
        plots[args.command](f, args)
    plt.show()
