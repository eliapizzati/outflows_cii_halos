# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 11:38:20 2021

@author: anna
"""




import numpy as np
import matplotlib.pyplot as plt
import natconst as nc


from my_utils import get_circular_velocity_profile_NFW_and_disk, get_circular_velocity_profile_NFW, \
                        get_vc_from_virial_mass, get_virial_radius, get_mass_profile_disk, mstar_behroozi
from cmb_suppression import eta_func, T_spin_func

import matplotlib



matplotlib.rcParams.update({
        "font.size": 16.0,
        "font.family": 'sans-serif',
#        "font.sans-serif": ['Helvetica'],
        "axes.titlesize": 14.,
        "axes.labelsize": 14.,
        "xtick.labelsize": 14.,
        "ytick.labelsize": 14.,
        "xtick.major.size": 4.0,
        "ytick.major.size": 4.0,
        "xtick.minor.size": 2.0,
        "ytick.minor.size": 2.0,
        "legend.fontsize": 13.0,
        "legend.frameon": False
 #       "figure.dpi": 200,
 #       "savefig.dpi": 200,
#        "text.usetex": True
})

plot1 = True  # gravity
plot2 = False      # exp profile

if plot1:
    """
    # NFW vc + comparison
    """

    params = dict([("DM_model", "NFW"),
                   ("beta", 5.0),
                   ("SFR", 50.),
                   ("f_esc_ion", 0.),
                   ("f_esc_FUV", 0.),
                   ("v_c", 250.),
                   ("redshift", 5.),
                   ("Zeta", 1.0),
                   ("alfa", 1.0),
                   ("R_in", 0.3)])

    fig, [ax_nfw, ax_v] = plt.subplots(1,2,figsize=(1.5*8.27,1.4*3.2))


    from sol_modules import get_profiles

    M_vir = np.linspace(1., 9, 4)



    cmap_rend_col = matplotlib.cm.get_cmap('viridis')

    norm = matplotlib.colors.Normalize(vmin=M_vir.min(), vmax=M_vir.max())

    ax_v.set_ylabel(r"$v(r)$ [km/s]")
    ax_v.set_xlabel("log (r [kpc])")

    ax_v.set_xlim((np.log10(0.3), np.log10(30)))
    ax_v.set_xlim((np.log10(0.3), np.log10(30)))

    for i in range(len(M_vir)):
        print("M_vir=", M_vir[i])

        v_c = get_vc_from_virial_mass(M_vir[i] * 1e11, params["redshift"])/1e5
        print(v_c)

        params.update(DM_model="NFW")

        params.update(v_c=v_c)

        profiles = get_profiles(params, resol=1000)

        r = profiles.r
        v = profiles.v

        ax_v.plot(np.log10(r / (1000 * nc.pc)), v / 10 ** 5,
                  color=cmap_rend_col((M_vir[i] - M_vir.min()) / (M_vir.max() - M_vir.min())))

        params.update(DM_model="NFW+disk")

        params.update(v_c=v_c)

        profiles = get_profiles(params, resol=1000)

        r = profiles.r
        v = profiles.v

        ax_v.plot(np.log10(r / (1000 * nc.pc)), v / 10 ** 5,
                  color=cmap_rend_col((M_vir[i] - M_vir.min()) / (M_vir.max() - M_vir.min())), linestyle="-.")



        params.update(DM_model=None)

        profiles = get_profiles(params, resol=1000)

        params.update(v_c=v_c)


        r = profiles.r
        v = profiles.v

        alpha = 0.3

        ax_v.plot(np.log10(r / (1000 * nc.pc)), v / 10 ** 5,
                 color=cmap_rend_col((M_vir[i] - M_vir.min()) / (M_vir.max() - M_vir.min())), alpha=alpha, linestyle="--")







    cax_cc = plt.axes([0.915, 0.15, 0.014, 0.8])


    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])

    cb = fig.colorbar(cmap, ticks=M_vir, cax=cax_cc, orientation='vertical')
    cb.set_label(r'$M_{vir}$ [10$^{11}$ M$_\odot$]', rotation=90., labelpad=5)
    cb.ax.tick_params(length=4)

    cb.set_ticks(np.arange(1.0, 10., 1.0))
    cb.set_ticklabels(np.arange(1.0, 10., 1.0))



    # NFW

    radius = np.linspace(0.3, 15, 1000)

    for i in range(len(M_vir)):
        v_c = get_circular_velocity_profile_NFW(radius * nc.pc * 1e3, M_vir[i] * 1e11, z=5)
        v_c_disk = get_circular_velocity_profile_NFW_and_disk(radius * nc.pc * 1e3, M_vir[i] * 1e11, z=5)

        v_c_global = get_vc_from_virial_mass(M_vir[i] * 1e11, z=5)
        #print(get_virial_radius(M_vir[i] * 1e11, z=5)/nc.pc/1e3)

        ax_nfw.plot(np.log10(radius), v_c / 1e5, color=cmap_rend_col((M_vir[i] - M_vir.min()) / (M_vir.max() - M_vir.min())))
        ax_nfw.plot(np.log10(radius), v_c_disk / 1e5, color=cmap_rend_col((M_vir[i] - M_vir.min()) / (M_vir.max() - M_vir.min())),\
                    linestyle ="-.")

        alpha = 0.3
        ax_nfw.plot(np.log10(radius), [v_c_global / 1e5] * len(radius),
                    color=cmap_rend_col((M_vir[i] - M_vir.min()) / (M_vir.max() - M_vir.min())), \
                    linestyle="--", alpha = 0.3)

    ax_nfw.set_ylabel(r"$v_c(r)$ [km/s]")
    ax_nfw.set_xlabel("log( r [kpc] )")
    ax_nfw.tick_params(length=4, which="both")

    cax_nfw = plt.axes([0.075, 0.15, 0.014, 0.8])

    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])

    cb = fig.colorbar(cmap, ticks=M_vir, cax=cax_nfw, orientation='vertical')
    cb.set_label(r'$M_{vir}$ [10$^{11}$ M$_\odot$]', rotation=90., labelpad=5)
    cb.ax.tick_params(length=4)

    cb.set_ticks(np.arange(1.0, 10., 1.0))
    cb.set_ticklabels(np.arange(1.0, 10., 1.0))

    cax_nfw.yaxis.set_ticks_position('left')
    cax_nfw.yaxis.set_label_position('left')



    plt.subplots_adjust(left=0.165,  # the left side of the subplots of the figure
                        right=0.9,  # the right side of the subplots of the figure
                        bottom=0.15,  # the bottom of the subplots of the figure
                        top=0.95,  # the top of the subplots of the figure
                        wspace=0.2,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace=0.1)  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height

    plt.show()




if plot2:
    """
    #check exp
    """

    M_vir = 7e11
    M_star = mstar_behroozi(M_vir)

    fig, ax = plt.subplots()

    ax.set_xlabel(r"r [kpc]")
    ax.set_ylabel(r"log( n [g/cm^3] )")

    n_points = 1000

    r = np.linspace(0.1, 10 * nc.pc * 1e3, n_points)

    R = 0.3 * nc.pc * 1e3

    x = r/R

    ax.plot(r / 1e3 / nc.pc, np.log10((M_star*nc.ms)/(8*np.pi*R**3)*np.exp(-x)), color = "gray")

    ax.axvline(R/1e3/nc.pc, color = "gray", linestyle = "--")

    ax_vc = ax.twinx()

    ax_vc.set_ylabel(r"v_c [km/s]")

    ax_vc.plot(r / 1e3 / nc.pc, np.sqrt(nc.gg * get_mass_profile_disk(r, M_star=M_star) * nc.ms / r) / 1e5,\
               label="disk")
    ax_vc.plot(r / 1e3 / nc.pc, get_circular_velocity_profile_NFW_and_disk(r, M_vir=M_vir, z=5)/ 1e5,
               label="disk + NFW")
    ax_vc.plot(r / 1e3 / nc.pc, get_circular_velocity_profile_NFW(r, M_vir=M_vir, z=5)/ 1e5,
               label="NFW")
    ax_vc.axhline(get_vc_from_virial_mass(M_vir=M_vir, z=5) / 1e5,
               label="iso sphere", linestyle = "--", color="C6")

    ax_vc.legend()

    plt.subplots_adjust(left=0.15,  # the left side of the subplots of the figure
                        right=0.85,  # the right side of the subplots of the figure
                        bottom=0.15,  # the bottom of the subplots of the figure
                        top=0.95,  # the top of the subplots of the figure
                        wspace=0.43,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace=0.1)  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height

    plt.show()

