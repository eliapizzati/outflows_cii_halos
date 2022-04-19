import numpy as np
import matplotlib.pyplot as plt

import matplotlib
import natconst as nc


from sol_modules import get_profiles
from post_sol_modules import get_ionization_states

params = dict([("DM_model", "NFW+disk"),
               ("beta", 5.0),
               ("SFR", 50.),
               ("f_esc_ion", 0.),
               ("f_esc_FUV", 0.),
               ("v_c", 200.),
               ("redshift", 5.),
               ("Zeta", 1.0),
               ("alfa", 1.0),
               ("R_in", 0.3)])


profiles = get_profiles(params, resol=1000)


#r = profiles.r
v = profiles.v
n = profiles.n
T = abs(profiles.T)

delta_r = profiles.r[1] - profiles.r[0]

ionization_state = get_ionization_states(profiles, params)

x_e = ionization_state.x_e
x_CII = ionization_state.x_CII

# computing the emissivity

# n_CII = nc.A_C * Zeta * x_CII * n
Zeta = 1

epsilon = 7.9e-20 * n ** 2 * (nc.A_C * Zeta) * nc.A_H * x_e * x_CII * np.exp(-92. / T) / 92. ** 0.5
print(epsilon)
cos_thetas = np.linspace(-1.,1., 1000)
delta_cos_theta = cos_thetas[1] - cos_thetas[0]

velocities_bins = np.linspace(-500,500, 51)

bins_lo = velocities_bins[:-1]
bins_hi = velocities_bins[1:]

velocities  = 0.5 * (bins_hi + bins_lo)


flux_v = np.zeros_like(velocities)

#print(profiles.r[100])

for v_index, v in enumerate(velocities):
    for r_index, r in enumerate(profiles.r):

    #print(velocities_bins[v_index],  velocities_bins[v_index+1], v)
        for cos_theta in cos_thetas:
            if profiles.v[r_index] * cos_theta/1e5 > velocities_bins[v_index]   and profiles.v[r_index] * cos_theta/1e5 < velocities_bins[v_index+1]:
                flux_v[v_index] += 2*np.pi*delta_cos_theta*r**2*delta_r*epsilon[r_index]


plt.plot(velocities, flux_v)

print(flux_v)
plt.show()