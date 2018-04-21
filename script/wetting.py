from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np


class Profile:
    def __init__(self, file):
        self.rootgrp = Dataset(file, "r", format="NETCDF4")

    def get_wetting_frac(self, ncut=0):
        return np.array(self.rootgrp.variables["wetting_fraction"])[ncut:]

    def show_wetting_frac(self):
        t = self.rootgrp.variables["time"]
        wetting_frac = self.get_wetting_frac()
        # print(wetting_frac[0])
        plt.plot(t, wetting_frac[:, 0])
        plt.plot(t, wetting_frac[:, 1])
        plt.show()
        plt.close()


def plot_wetting_frac_vs_alpha():
    # tumbling rate
    alpha = np.array([0.005, 0.01, 0.02, 0.03])
    # order para for the wetting to partial-wetting trasition
    Zm = np.zeros_like(alpha)
    ncut = 300
    for i in range(alpha.size):
        filename = r"D:/code/Swimmer/BD2D_MPI/data/profile_%g_0.04_10.nc" % (
            alpha[i])
        prof = Profile(filename)
        Zm[i] = 1 - np.mean(prof.get_wetting_frac(ncut))
        print(alpha[i], Zm[i])
    plt.loglog(alpha, Zm, "o")
    plt.show()
    plt.close()


if __name__ == "__main__":
    # path = r"D:/code/Swimmer/BD2D_MPI/data/0.005_0.04_10/"
    # x = Profile(path + "profile.nc")
    # x.show_wetting_frac()
    plot_wetting_frac_vs_alpha()
