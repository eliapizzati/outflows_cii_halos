# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 11:38:20 2021

@author: anna
"""




import numpy as np
import matplotlib.pyplot as plt
import natconst as nc


from my_utils import get_circular_velocity_profile_NFW_and_disk, get_vc_from_virial_mass, get_virial_radius,\
                     get_circular_velocity_profile_NFW, to_habing_flux
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


plot1 = False  # gravity
plot2 = False      # cmb
plot3 = False # gravity horizontla
plot4 = True # cmb horizontal


params = dict([("DM_model", "NFW+disk"),
               ("beta", 1.0),
               ("SFR", 50.),
               ("f_esc_ion", 0.),
               ("f_esc_FUV", 0.),
               ("v_c", 200.),
               ("redshift", 5.),
               ("Zeta", 1.0),
               ("alfa", 1.0),
               ("R_in", 0.3)])

#betas = np.arange(0.2, 8, 3.2)
betas = np.asarray([0.2, 0.4,  1., 1.5, 2.4,4.0, 6.6])

if plot1:
    """
    # NFW vc + comparison
    """



    fig, [ax_nfw, ax_v] = plt.subplots(1,2,figsize=(1.5*8.27,1.4*3.2))


    from sol_modules import get_profiles

    cmap_rend_col = matplotlib.cm.get_cmap('viridis_r')

    norm = matplotlib.colors.Normalize(vmin=0., vmax=7.)

    ax_v.set_ylabel("v [1000 km/s] ")
    ax_v.set_xlabel("log (r [kpc])")

    ax_v.set_xlim((np.log10(0.3), np.log10(30)))
    ax_v.set_xlim((np.log10(0.3), np.log10(30)))

    for i in range(len(betas)):
        print("beta=", betas[i])
        params.update(DM_model="NFW+disk")

        params.update(beta=betas[i])

        profiles = get_profiles(params, resol=1000)

        r = profiles.r
        v = profiles.v

        ax_v.plot(np.log10(r / (1000 * nc.pc)), v / 10 ** 8,
                  color=cmap_rend_col((betas[i] - betas.min()) / (betas.max() - betas.min())))


        params.update(DM_model=None)

        profiles = get_profiles(params, resol=1000)

        r = profiles.r
        v = profiles.v

        alpha = 1.0


        ax_v.plot(np.log10(r / (1000 * nc.pc)), v / 10 ** 8,
                 color=cmap_rend_col((betas[i] - betas.min()) / (betas.max() - betas.min())), alpha=alpha, linestyle="--")



        print("beta=", betas[i])




    cax_cc = plt.axes([0.915, 0.15, 0.014, 0.8])


    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])

    cb = fig.colorbar(cmap, orientation='vertical', cax=cax_cc)
    cb.set_label(r'$\eta$', rotation=90., labelpad=3)
    cb.set_ticks(np.arange(1.0, 10., 1.0))
    cb.set_ticklabels(np.arange(1.0, 10., 1.0))


    # NFW

    M_vir = np.linspace(1, 10, 11)

    cmap_rend_col = matplotlib.cm.get_cmap('inferno_r')

    norm = matplotlib.colors.Normalize(vmin=M_vir.min(), vmax=M_vir.max())

    radius = np.linspace(1, 55, 10000)


    for i in range(len(M_vir)):
        v_c = get_circular_velocity_profile_NFW_and_disk(radius * nc.pc * 1e3, M_vir[i] * 1e11, z=5)
        v_c_global = get_vc_from_virial_mass(M_vir[i] * 1e11, z=5)

        rvir = get_virial_radius(M_vir[i] * 1e11, z=5)/nc.pc/1e3

        ax_nfw.plot(radius[radius<rvir], v_c[radius<rvir] / 1e5, color=cmap_rend_col((M_vir[i] - M_vir.min()) / (M_vir.max() - M_vir.min())))
        ax_nfw.plot(radius, [v_c_global / 1e5] * len(radius),
                    color=cmap_rend_col((M_vir[i] - M_vir.min()) / (M_vir.max() - M_vir.min())), \
                    linestyle="--")

    ax_nfw.set_ylabel(r"$v_c(r)$ [km/s]")
    ax_nfw.set_xlabel("r [kpc]")
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
    #cmb emission
    """
    fig, [ax_texc, ax_sigma] = plt.subplots(1,2,figsize=(1.5*8.27,1.4*3.2))



    ax_texc.set_xlabel(r"log( $n_e$ [cm$^{-3}$] )")
    ax_texc.set_ylabel(r"log( I$_{UV}$ [erg s$^{-1}$cm$^{-2}$Hz$^{-1}$sr$^{-1}$] )")


    n_points = 1000

    n_es = np.linspace(-5, 3, n_points)
    I_UVs = np.linspace(-11, -19, n_points)

    nn, II = np.meshgrid(n_es, I_UVs)

    z = 5

    T_spin = T_spin_func(n=10**nn / 0.5, T=1e3, I_UV=10**II, x_e=0.5, z=z)
    z=5.
    T_CMB_z = nc.CMB_temperature * (1 + z)

    im = ax_texc.imshow(np.log10(T_spin / T_CMB_z), cmap=matplotlib.cm.inferno, aspect="auto", \
                        norm=matplotlib.colors.PowerNorm(gamma=0.3), interpolation="nearest",\
                        extent=[-5-8/2/n_points,3+8/2/n_points,-19-8/2/n_points,-11+8/2/n_points])

    # Create colorbar

    cax = plt.axes([0.06, 0.15, 0.014, 0.8])

    ticks = [0.01, 0.1,0.4,1.0,2.0]
    cbar = ax_texc.figure.colorbar(im, cax=cax, ticks = ticks)
    cbar.ax.set_ylabel(r"log( T$_{exc}$ / T$_{\rm CMB}$(z) )", rotation=90., labelpad=0)
    cbar.ax.set_yticklabels(ticks)

    cax.yaxis.set_ticks_position('left')
    cax.yaxis.set_label_position('left')

    # ax_texc.set_title("T_spin dependence on n_e, I_UV; (T = 10^3 K, z = 5.)")

    ax2 = ax_texc.twinx()

    ax_texc.set_ylim(-19,-11)
    ax2.set_ylim(np.log10(to_habing_flux(1e-19)), np.log10(to_habing_flux(1e-11)))

    habing_g0 = 4 * np.pi * 10**I_UVs * (13.6 - 6.) * nc.ev / nc.hh / nc.cc / 5.29e-14

    ax2.set_yticks([1.,2.,3.,4.,5.,6.,7.,8.])
    ax2.set_ylabel(r"log G$_0 $")

    #ax_texc.text(0.17, 0.87, s=r'z=5', fontsize=16, bbox={'facecolor': 'lightgray', 'alpha': 0.8, \
    #                                                      'boxstyle': "round", 'edgecolor': 'none', 'pad': 0.5}, \
    #             horizontalalignment='center', verticalalignment='center', \
    #             transform=ax_texc.transAxes)


    from sol_modules import get_profiles
    from post_sol_modules import get_ionization_states, get_surface_density

    cmap_rend_col = matplotlib.cm.get_cmap('viridis_r')

    norm = matplotlib.colors.Normalize(vmin=0., vmax=7.)



    ax_sigma.set_xlabel("b [kpc]")
    ax_sigma.set_ylabel(r"log ($\Sigma_{\rm CII}$ [erg cm$^{-2}$ s$^{-1}$])", labelpad=-1)
    ax_sigma.set_yscale("log")
    ax_sigma.set_xlim((0.3, 15.5))
    ax_sigma.set_ylim(1e-7, 5e-2)

    for i in range(len(betas)):
        print("beta=", betas[i])

        params.update(beta=betas[i])

        profiles = get_profiles(params, resol=1000)
        ionization_state = get_ionization_states(profiles, params)
        sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500,
                                        add_CMB_suppression=True)
        sigma_CII_no_cmb = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500,
                                               add_CMB_suppression=False)

        ax_sigma.plot(sigma_CII.h / (1e3 * nc.pc), sigma_CII.var,
                      color=cmap_rend_col((betas[i] - betas.min()) / (betas.max() - betas.min())))
        ax_sigma.plot(sigma_CII_no_cmb.h / (1e3 * nc.pc), sigma_CII_no_cmb.var,
                      color=cmap_rend_col((betas[i] - betas.min()) / (betas.max() - betas.min())), linestyle="--")


    cax_cc = plt.axes([0.925, 0.15, 0.014, 0.8])


    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])

    cb = fig.colorbar(cmap, orientation='vertical', cax=cax_cc,ticks=betas)
    cb.set_label(r'$\eta$', rotation=90., labelpad=3)
    cb.set_ticks(np.arange(1.0, 10., 1.0))
    cb.set_ticklabels(np.arange(1.0, 10., 1.0))






    plt.subplots_adjust(left=0.15,  # the left side of the subplots of the figure
                        right=0.91,  # the right side of the subplots of the figure
                        bottom=0.15,  # the bottom of the subplots of the figure
                        top=0.95,  # the top of the subplots of the figure
                        wspace=0.43,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace=0.1)  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height

    plt.show()




if plot3:
    """
    # NFW vc + comparison
    """



    fig, [ax_nfw, ax_v] = plt.subplots(1,2,figsize=(1.5*8.27,1.4*3.7))


    from sol_modules import get_profiles

    cmap_rend_col = matplotlib.cm.get_cmap('viridis_r')

    norm = matplotlib.colors.Normalize(vmin=0., vmax=7.)

    ax_v.set_ylabel("v [1000 km/s] ")
    ax_v.set_xlabel("log (r [kpc])")

    ax_v.set_xlim((np.log10(0.3), np.log10(30)))
    ax_v.set_xlim((np.log10(0.3), np.log10(30)))

    for i in range(len(betas[::])):
        print("beta=", betas[i])
        params.update(DM_model="NFW+disk")

        params.update(beta=betas[i])

        profiles = get_profiles(params, resol=1000)

        r = profiles.r
        v = profiles.v

        ax_v.plot(np.log10(r / (1000 * nc.pc)), v / 10 ** 8,
                  color=cmap_rend_col((betas[i] - betas.min()) / (betas.max() - betas.min())))


        params.update(DM_model=None)

        profiles = get_profiles(params, resol=1000)

        r = profiles.r
        v = profiles.v

        alpha = 1.0


        ax_v.plot(np.log10(r / (1000 * nc.pc)), v / 10 ** 8,
                 color=cmap_rend_col((betas[i] - betas.min()) / (betas.max() - betas.min())), alpha=alpha, linestyle="--")



        print("beta=", betas[i])




    cax_cc = plt.axes([0.582, 0.93, 0.383, 0.03])


    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])

    cb = fig.colorbar(cmap, orientation='horizontal', cax=cax_cc)
    cb.set_label(r'$\eta$', rotation=0., labelpad=2)
    cb.set_ticks(np.arange(1.0, 7., 1.0))
    cb.set_ticklabels(np.arange(1.0, 7., 1.0))



    # NFW

    M_vir = np.linspace(1, 10, 11)

    cmap_rend_col = matplotlib.cm.get_cmap('inferno_r')

    norm = matplotlib.colors.Normalize(vmin=M_vir.min(), vmax=M_vir.max())

    radius = np.linspace(1, 55, 10000)


    for i in range(len(M_vir)):
        v_c = get_circular_velocity_profile_NFW_and_disk(radius * nc.pc * 1e3, M_vir[i] * 1e11, z=5)
        v_c_global = get_vc_from_virial_mass(M_vir[i] * 1e11, z=5)

        rvir = get_virial_radius(M_vir[i] * 1e11, z=5)/nc.pc/1e3

        ax_nfw.plot(radius[radius<rvir], v_c[radius<rvir] / 1e5, color=cmap_rend_col((M_vir[i] - M_vir.min()) / (M_vir.max() - M_vir.min())))
        ax_nfw.plot(radius, [v_c_global / 1e5] * len(radius),
                    color=cmap_rend_col((M_vir[i] - M_vir.min()) / (M_vir.max() - M_vir.min())), \
                    linestyle="--")

    ax_nfw.set_ylabel(r"$v_c(r)$ [km/s]")
    ax_nfw.set_xlabel("r [kpc]")
    ax_nfw.tick_params(length=4, which="both")

    cax_nfw = plt.axes([0.08, 0.93, 0.383, 0.03])

    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])

    cb = fig.colorbar(cmap, ticks=M_vir, cax=cax_nfw, orientation='horizontal')
    cb.set_label(r'$M_{vir}$ [10$^{11}$ M$_\odot$]', rotation=0., labelpad=2)
    cb.ax.tick_params(length=4)

    cb.set_ticks(np.arange(2.0, 10., 1.0))
    cb.set_ticklabels(np.arange(2.0, 10., 1.0))



    plt.subplots_adjust(left=0.08,  # the left side of the subplots of the figure
                        right=0.97,  # the right side of the subplots of the figure
                        bottom=0.15,  # the bottom of the subplots of the figure
                        top=0.8,  # the top of the subplots of the figure
                        wspace=0.3,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace=0.1)  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height

    plt.show()


if plot4:
    """
    #cmb emission
    """
    fig, [ax_texc, ax_sigma] = plt.subplots(1,2,figsize=(1.5*8.27,1.4*3.7))



    ax_texc.set_xlabel(r"log( $n_e$ [cm$^{-3}$] )")
    ax_texc.set_ylabel(r"log( I$_{UV}$ [erg s$^{-1}$cm$^{-2}$Hz$^{-1}$sr$^{-1}$] )", labelpad=-1)


    n_points = 1000

    n_es = np.linspace(-5, 3, n_points)
    I_UVs = np.linspace(-11, -19, n_points)

    nn, II = np.meshgrid(n_es, I_UVs)

    z = 5

    T_spin = T_spin_func(n=10**nn / 0.5, T=1e3, I_UV=10**II, x_e=0.5, z=z)
    z=5.
    T_CMB_z = nc.CMB_temperature * (1 + z)

    im = ax_texc.imshow(np.log10(T_spin / T_CMB_z), cmap=matplotlib.cm.inferno, aspect="auto", \
                        norm=matplotlib.colors.PowerNorm(gamma=0.3), interpolation="nearest",\
                        extent=[-5-8/2/n_points,3+8/2/n_points,-19-8/2/n_points,-11+8/2/n_points])

    # Create colorbar

    cax = plt.axes([0.08, 0.93, 0.37, 0.03])

    ticks = [0.01, 0.1,0.4,1.0,2.0]
    cbar = ax_texc.figure.colorbar(im, cax=cax, ticks = ticks, orientation='horizontal')
    cbar.ax.set_xlabel(r"log( T$_{exc}$ / T$_{\rm CMB}$(z) )", rotation=0., labelpad=2)
    cbar.ax.set_xticklabels(ticks)


    # ax_texc.set_title("T_spin dependence on n_e, I_UV; (T = 10^3 K, z = 5.)")

    ax2 = ax_texc.twinx()

    ax_texc.set_ylim(-19,-11)
    ax2.set_ylim(np.log10(to_habing_flux(1e-19)), np.log10(to_habing_flux(1e-11)))

    habing_g0 = 4 * np.pi * 10**I_UVs * (13.6 - 6.) * nc.ev / nc.hh / nc.cc / 5.29e-14

    ax2.set_yticks([1.,2.,3.,4.,5.,6.,7.,8.])
    ax2.set_ylabel(r"log G$_0 $")

    #ax_texc.text(0.17, 0.87, s=r'z=5', fontsize=16, bbox={'facecolor': 'lightgray', 'alpha': 0.8, \
    #                                                      'boxstyle': "round", 'edgecolor': 'none', 'pad': 0.5}, \
    #             horizontalalignment='center', verticalalignment='center', \
    #             transform=ax_texc.transAxes)


    from sol_modules import get_profiles
    from post_sol_modules import get_ionization_states, get_surface_density

    cmap_rend_col = matplotlib.cm.get_cmap('viridis_r')

    norm = matplotlib.colors.Normalize(vmin=0., vmax=7.)



    ax_sigma.set_xlabel("b [kpc]")
    ax_sigma.set_ylabel(r"log ($\Sigma_{\rm CII}$ [erg cm$^{-2}$ s$^{-1}$])", labelpad=-1)
    ax_sigma.set_yscale("log")
    ax_sigma.set_xlim((0.3, 15.5))
    ax_sigma.set_ylim(1e-7, 5e-2)

    for i in range(len(betas)):
        print("beta=", betas[i])

        params.update(beta=betas[i])

        profiles = get_profiles(params, resol=1000)
        ionization_state = get_ionization_states(profiles, params)
        sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500,
                                        add_CMB_suppression=True)
        sigma_CII_no_cmb = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500,
                                               add_CMB_suppression=False)

        ax_sigma.plot(sigma_CII.h / (1e3 * nc.pc), sigma_CII.var,
                      color=cmap_rend_col((betas[i] - betas.min()) / (betas.max() - betas.min())))
        ax_sigma.plot(sigma_CII_no_cmb.h / (1e3 * nc.pc), sigma_CII_no_cmb.var,
                      color=cmap_rend_col((betas[i] - betas.min()) / (betas.max() - betas.min())), linestyle="--")


    cax_cc = plt.axes([0.597, 0.93, 0.37, 0.03])


    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])

    cb = fig.colorbar(cmap, orientation='horizontal', cax=cax_cc,ticks=betas)
    cb.set_label(r'$\eta$', rotation=0., labelpad=2)
    cb.set_ticks(np.arange(1.0, 7., 1.0))
    cb.set_ticklabels(np.arange(1.0, 7., 1.0))





    plt.subplots_adjust(left=0.08,  # the left side of the subplots of the figure
                        right=0.97,  # the right side of the subplots of the figure
                        bottom=0.15,  # the bottom of the subplots of the figure
                        top=0.8,  # the top of the subplots of the figure
                        wspace=0.4,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace=0.1)  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height


    plt.show()
