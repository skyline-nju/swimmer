from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os


def read_dat(file, ncut=0):
    with open(file, "r") as f:
        lines = f.readlines()[ncut:]
        frac = np.zeros((len(lines), 2))
        for i, line in enumerate(lines):
            s = line.replace("\n", "").split("\t")
            frac[i, 0] = float(s[1])
            frac[i, 1] = float(s[2])
    return frac


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


def plot_wetting_frac_vs_alpha(k, ncut=1000):
    # tumbling rate
    alpha = np.array([
        0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.013, 0.016, 0.02, 0.025,
        0.03
    ])
    # order para for the wetting to partial-wetting trasition
    Zm = np.zeros_like(alpha)
    for i in range(alpha.size):
        filename = r"profile_%g_%d_0.04.nc" % (alpha[i], k)
        if os.path.exists(filename):
            prof = Profile(filename)
            frac_arr = prof.get_wetting_frac(ncut)
        else:
            frac_arr = read_dat(filename.replace("nc", "dat"), ncut)
        Zm[i] = 1 - np.mean(frac_arr)
        print(alpha[i], Zm[i])
    plt.loglog(alpha, Zm, "o")
    plt.show()
    plt.close()


if __name__ == "__main__":
    # path = r"D:/code/Swimmer/BD2D_MPI/data/0.005_0.04_10/"
    # x = Profile(path + "profile.nc")
    # x.show_wetting_frac()
    os.chdir(r"E:\data\roughening\wetting\1500_1500")
    # pro = Profile("profile_0.005_15_0.04.nc")
    # pro.show_wetting_frac()
    plot_wetting_frac_vs_alpha(5, 1500)
