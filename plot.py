import argparse
import matplotlib.pyplot as plt
import numpy as np
import h5py
import gusto_dataset



def plot_ascii_1d():
    A = np.loadtxt('gusto.dat')
    x3 = A[:,0]
    dg = A[:,9]

    plt.plot(x3, dg, '-o')
    plt.show()



def plot_h5_1d(filenames):
    zloc = [ ]
    for f in filenames:
        h5f = h5py.File(f, 'r')
        z = h5f['rows']['row_000000']['x3'][1:-1]
        y = h5f['rows']['row_000000']['dg'][1:-1]
        zloc.append( z[len(z) / 2] )
        h5f.close()
        #plt.plot(z, y, '-o')
    plt.plot(zloc, '-o')
    plt.show()



def plot_h5():
    h5f = h5py.File('chkpt.0000.h5', 'r')
    x = [ ]
    y = [ ]
    z = [ ]

    for row_name in h5f['rows']:
        x1 = h5f['rows'][row_name]['x1'][:]
        x3 = h5f['rows'][row_name]['x3'][:]
        dg = h5f['rows'][row_name]['dg'][:]
        x.append(x1)
        y.append(x3)
        z.append(dg)
    h5f.close()

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    plt.pcolormesh(x, y, z)
    plt.axis('equal')
    plt.show()



def plot_faces(filename):
    dset = gusto_dataset.GustoDataset(filename)
    segments = dset.get_face_segments()

    for seg in segments:
        plt.plot(seg[:,0], seg[:,2], '-o', c='k')

    plt.axis('equal')
    plt.show()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    args = parser.parse_args()

    #plot_h5_1d(args.filenames)
    plot_faces(args.filenames[0])
