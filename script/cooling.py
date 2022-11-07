# -*- coding: utf-8 -*-
"""
Created on Mon May  3 10:50:09 2021

@author: anna
"""

import itertools
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, colors

import mydir

import plot_config

import gnedincooling as gc

gc.frtinitcf(0, os.path.join(mydir.script_dir, "input_data", "cf_table.I2.dat"))


save_to_file = False
save_plot = True
load_from_file = False
plotting = True


if __name__ == "__main__":
    quantity: str = input("which quantity (n,plw,ph1,pg1,pc6)?")

    n = float(input("n value [default to unity] : ") or "1.")
    plw = float(input("Plw value [default to zero] : ") or "0.")
    ph1 = float(input("Ph1 value [default to zero] : ") or "0.")
    pg1 = float(input("Pg1 value [default to zero] : ") or "0.")
    pc6 = float(input("Pc6 value [default to zero] : ") or "0.")

    # r1 = plw
    # r2 = (ph1/plw)**0.353 * (pg1/plw)**0.923 * (pc6/plw)**0.263
    # r3 = (ph1/plw)**-0.103 * (pg1/plw)**-0.375 * (pc6/plw)**0.976

    # print("#############################")
    # print("r1 =", r1)
    # print("r2 =", r2)
    # print("r3 =", r3)
    # print("#############################")

    values = np.asarray([n, plw, ph1, pg1, pc6])

temperatures = np.logspace(3, 8, 1000)

densities = np.logspace(-5, 3, 9)

plws = np.logspace(-15, -5, 11)
ph1s = np.logspace(-15, -5, 11)
pg1s = np.logspace(-15, -5, 11)

pc6s = np.logspace(-25, -15, 11)


# plws = np.insert(plws, 0, 0.)
# ph1s = np.insert(ph1s, 0, 0.)
# pg1s = np.insert(pg1s, 0, 0.)


def cooling(T, n, Plw, Ph1, Pg1, Pc6=0.):
    """
        Cooling function (calling gnedincooling from Fortran functions)
        
        Parameters
        ==========
        T: float
            temperature
        n: float
            density
        Plw: float
            Lyman-Wernner PD rate
        Ph1: float
            H photoionization rate
        Pg1: float
            He photoionization rate
        """

    Zeta = 1.

    return gc.frtgetcf_cool(T, n, Zeta, Plw, Ph1, Pg1, Pc6)


def heating(T, n, Plw, Ph1, Pg1, Pc6=0.):
    """
        Heating function (calling gnedincooling from Fortran functions)
        
        Parameters
        ==========
        T: float
            temperature
        n: float
            density
        Plw: float
            Lyman-Wernner PD rate
        Ph1: float
            H photoionization rate
        Pg1: float
            He photoionization rate    
        """

    Zeta = 1.

    return gc.frtgetcf_heat(T, n, Zeta, Plw, Ph1, Pg1, Pc6)


def cooling_vectorized(T_vec, n, Plw, Ph1, Pg1, Pc6):
    cooling_values = [cooling(T, n, Plw, Ph1, Pg1, Pc6) for T in T_vec]

    print("INFO: loading cooling function with values: n={:.1e}, plw={:.1e}, ph1={:.1e}, pg1={:.1e}, pc6={:.1e}". \
          format(n, Plw, Ph1, Pg1, Pc6))

    return T_vec, np.asarray(cooling_values)


def heating_vectorized(T_vec, n, Plw, Ph1, Pg1, Pc6):
    heating_values = [heating(T, n, Plw, Ph1, Pg1, Pc6) for T in T_vec]

    print("INFO: loading heating function with values: n={:.1e}, plw={:.1e}, ph1={:.1e}, pg1={:.1e}, pc6={:.1e}". \
          format(n, Plw, Ph1, Pg1, Pc6))

    return T_vec, np.asarray(heating_values)


def to_file(name_base, var_vec, T_vec, n, Plw, Ph1, Pg1, Pc6):
    folder = 'data_cooling'

    if not os.path.exists(os.path.join(mydir.data_dir, folder)):
        os.mkdir(os.path.join(mydir.data_dir, folder))

    name = "{}_n{:.1e}_plw{:.1e}_ph{:.1e}_pg{:.1e}_pc{:.1e}.dat" \
        .format(name_base, n, Plw, Ph1, Pg1, Pc6)

    path = os.path.join(mydir.data_dir, folder, name)

    np.savetxt(path, (T_vec, var_vec))

    print("INFO: saving data to {}".format(path))
    return


def from_file(name_base, n, Plw, Ph1, Pg1, Pc6):
    folder = 'data_cooling'

    name = "{}_n{:.1e}_plw{:.1e}_ph{:.1e}_pg{:.1e}_pc{:.1e}.dat" \
        .format(name_base, n, Plw, Ph1, Pg1, Pc6)

    path = os.path.join(mydir.data_dir, folder, name)

    T_vec, var_vec = np.loadtxt(path)

    print("INFO: loading data from {}".format(path))

    return T_vec, var_vec


def plot_quantity(quantity, quantity_array, values):
    if quantity == "n":
        n, plw, ph1, pg1, pc6 = values
        title = "Plw={:.1e}, Ph1={:.1e}, Pg1={:.1e}, Pc6={:.1e}".format(plw, ph1, pg1, pc6)
        label_cmap = r'density [cm^-3]'
        name = "plot_lamda_vs_T_vs_{}_plw{:.1e}_ph{:.1e}_pg{:.1e}_pc{:.1e}.png" \
            .format(quantity, plw, ph1, pg1, pc6)

    elif quantity == "plw":
        n, plw, ph1, pg1, pc6 = values
        title = "n={:.1e}, Ph1={:.1e}, Pg1={:.1e}, Pc6={:.1e}".format(n, ph1, pg1, pc6)
        label_cmap = r'plw [s^-1]'
        name = "plot_lamda_vs_T_vs_{}_n{:.1e}_ph{:.1e}_pg{:.1e}_pc{:.1e}.png" \
            .format(quantity, n, ph1, pg1, pc6)

    elif quantity == "ph1":
        n, plw, ph1, pg1, pc6 = values
        title = "n={:.1e}, Plw={:.1e}, Pg1={:.1e}, Pc6={:.1e}".format(n, plw, pg1, pc6)
        label_cmap = r'ph1 [s^-1]'
        name = "plot_lamda_vs_T_vs_{}_n{:.1e}_plw{:.1e}_pg{:.1e}_pc{:.1e}.png" \
            .format(quantity, n, plw, pg1, pc6)

    elif quantity == "pg1":
        n, plw, ph1, pg1, pc6 = values
        title = "n={:.1e}, Plw={:.1e}, Ph1={:.1e}, Pc6={:.1e}".format(n, plw, ph1, pc6)
        label_cmap = r'pg1 [s^-1]'
        name = "plot_lamda_vs_T_vs_{}_n{:.1e}_plw{:.1e}_ph{:.1e}_pc{:.1e}.png" \
            .format(quantity, n, plw, ph1, pc6)

    elif quantity == "pc6":
        n, plw, ph1, pg1, pc6 = values
        title = "n={:.1e}, Plw={:.1e}, Ph1={:.1e}, Pg1={:.1e}".format(n, plw, ph1, pg1)
        label_cmap = r'pc6 [s^-1]'
        name = "plot_lamda_vs_T_vs_{}_n{:.1e}_plw{:.1e}_ph{:.1e}_pg{:.1e}.png" \
            .format(quantity, n, plw, ph1, pg1)

    else:
        raise ValueError("No correct quantity given")

    fig, ax = plt.subplots(figsize=(8.27, 5.))

    cmap_rend_col = cm.get_cmap('viridis')

    v_min = min(quantity_array)
    v_max = max(quantity_array)

    norm = colors.LogNorm(vmin=v_min / 10 ** 0.5, vmax=v_max * 10 ** 0.5)

    if load_from_file == True:
        T_vec, cooling_vec = from_file("cooling", n, plw, ph1, pg1, pc6)
        T_vec, heating_vec = from_file("heating", n, plw, ph1, pg1, pc6)
    else:
        T_vec, cooling_vec = cooling_vectorized(temperatures, n, plw, ph1, pg1, pc6)
        T_vec, heating_vec = heating_vectorized(temperatures, n, plw, ph1, pg1, pc6)

    ax.plot(T_vec, cooling_vec, color="gray")
    ax.plot(T_vec, heating_vec, color="gray", linestyle=':')

    for q in quantity_array:

        if quantity == "n":
            n = q
        elif quantity == "plw":
            plw = q
        elif quantity == "ph1":
            ph1 = q
        elif quantity == "pg1":
            pg1 = q
        elif quantity == "pc6":
            pc6 = q

        if load_from_file == True:
            T_vec, cooling_vec = from_file("cooling", n, plw, ph1, pg1, pc6)
            T_vec, heating_vec = from_file("heating", n, plw, ph1, pg1, pc6)
        else:

            T_vec, cooling_vec = cooling_vectorized(temperatures, n, plw, ph1, pg1, pc6)
            T_vec, heating_vec = heating_vectorized(temperatures, n, plw, ph1, pg1, pc6)

        ax.plot(T_vec, cooling_vec,
                color=cmap_rend_col((np.log10(q) - np.log10(v_min)) / (np.log10(v_max) - np.log10(v_min))))
        ax.plot(T_vec, heating_vec,
                color=cmap_rend_col((np.log10(q) - np.log10(v_min)) / (np.log10(v_max) - np.log10(v_min))), \
                linestyle=':')

    ax.set_xlabel("T [K]")
    ax.set_ylabel(r"$|\Lambda$(T)| $[\mathrm{erg}\,\mathrm{cm}^3\,\mathrm{s}^{-1}]$")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim((1e-27, 1e-21))

    plt.subplots_adjust(left=0.15,  # the left side of the subplots of the figure
                        right=0.98,  # the right side of the subplots of the figure
                        bottom=0.15,  # the bottom of the subplots of the figure
                        top=0.9,  # the top of the subplots of the figure
                        wspace=0.1,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace=0.1)  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height

    cmap = cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])

    cb = fig.colorbar(cmap, orientation='vertical')

    ax.set_title(title)
    cb.set_label(label_cmap, rotation=90.)

    fig.tight_layout()

    if save_plot:
        folder = 'plot_cooling'

        if not os.path.exists(os.path.join(mydir.plot_dir, folder)):
            os.mkdir(os.path.join(mydir.plot_dir, folder))

        path = os.path.join(mydir.plot_dir, folder, name)

        plt.savefig(path)

        print("INFO: saving plot to {}".format(path))


if save_to_file == True:

    import gnedincooling as gc

    gc.frtinitcf(0, os.path.join(mydir.script_dir, "input_data", "cf_table.I2.dat"))

    print("###############################################")
    print("STARTING ITERATIONS")
    print("###############################################")

    i = 0

    for n, plw, ph1, pg1, pc6 in itertools.product(densities, plws, ph1s, pg1s, pc6s):
        temperatures, cooling_values = cooling_vectorized(temperatures, n, plw, ph1, pg1, pc6)
        temperatures, heating_values = heating_vectorized(temperatures, n, plw, ph1, pg1, pc6)

        to_file("cooling", cooling_values, temperatures, n, plw, ph1, pg1, pc6)
        to_file("heating", heating_values, temperatures, n, plw, ph1, pg1, pc6)

        i += 1
        print("ITERATION {} OF {}".format(i, densities.size * plws.size * ph1s.size * pg1s.size))

if plotting:

    if load_from_file == False:
        import gnedincooling as gc

        gc.frtinitcf(0, os.path.join(mydir.script_dir, "input_data", "cf_table.I2.dat"))

    if quantity == "n":
        quantity_array = densities

    elif quantity == "plw":
        quantity_array = plws

    elif quantity == "ph1":
        quantity_array = ph1s

    elif quantity == "pg1":
        quantity_array = pg1s

    elif quantity == "pc6":
        quantity_array = pc6s

    else:
        raise ValueError("No correct quantity given")
    plot_quantity(quantity, quantity_array, values=values)

plt.show()
