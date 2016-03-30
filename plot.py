import matplotlib.pyplot as plt
import numpy as np
import h5py

def plot_ascii_1d():
    A = np.loadtxt('gusto.dat')
    x3 = A[:,0]
    dg = A[:,9]

    plt.plot(x3, dg, '-o')
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

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    # for i in range(x.shape[0] / 2):
    #     plt.plot(x[i,x.shape[1]/2:], y[i,x.shape[1]/2:], c='k')

    # for j in range(x.shape[1] / 2):
    #     plt.plot(x[:-x.shape[0]/2,-j], y[:-x.shape[0]/2:,-j], c='k')

    plt.pcolormesh(x, y, z)
    plt.axis('equal')
    plt.show()


plot_ascii_1d()
