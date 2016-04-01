import argparse
import matplotlib.pyplot as plt
import numpy as np
import h5py
import gusto_dataset



def plot_faces(filename):
    dset = gusto_dataset.GustoDataset(filename)
    segments = dset.get_face_segments()

    for seg in segments:
        plt.plot(seg[:,0], seg[:,2], '-o', c='k')

    plt.axis('equal')



def plot_1d(filename):
    dset = gusto_dataset.GustoDataset(filename)
    x = dset.get_cell_variable('x3')[0]
    y = dset.get_cell_variable('dg')[0]
    plt.plot(x, y, '-o')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    args = parser.parse_args()

    #plot_faces(args.filenames[0])

    for filename in args.filenames:
        plot_1d(filename)

    plt.show()
