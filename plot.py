import matplotlib.pyplot as plt
import numpy as np
import h5py

def plot_ascii_1d():
    A = np.loadtxt('gusto.dat')
    x3 = A[:,0]
    u3 = A[:,4]
    dg = A[:,9]

    plt.plot(x3, u3)
    plt.plot(x3, dg)
    plt.show()



def plot_h5():
    h5f = h5py.File('chkpt.0000.h5', 'r')
    for row_name in h5f['rows']:
        x1 = h5f['rows'][row_name]['x1'][:]
        x3 = h5f['rows'][row_name]['x3'][:]
        plt.plot(x1, x3, '-', c='k')
    plt.show()


plot_ascii_1d()
