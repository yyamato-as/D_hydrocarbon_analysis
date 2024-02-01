import analysis_utils as au
from diskdictionary import disk_dict
from linedictionary_v2 import line_dict
from specdata import SpectroscopicData, PartitionFunction
from qdisk.classes import FitsImage
import numpy as np
from MCMC_utils import Parameter, ParameterSet, log_prior, EmceeHammer
import spectral_model as specmodel
import multiprocessing
import os
import pickle
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

##################################################### Fit set-up #################################################
# basic info
source = "MWC_480"
species = ["C2H", "C2D"]
lines = ["C2H_1-0_fs1", "C2H_3-2_fs1", "C2H_3-2_fs2", "C2D_4-3_fs1", "C2D_4-3_fs2"]
molecular_weight = {"C2H": 25, "C2D": 26}

# files of radially-resolved spectra
specfilename = f"./data/spectrum/{source}_C2H_C2D_radial_spectra_shift.pkl"
specfilename_wcont = f"./data/spectrum/{source}_C2H_C2D_radial_spectra.pkl"

# directory to store the fit result
fitdir = f"./data/fit/{source}/C2H_C2D_LTE/"
if not os.path.exists(fitdir):
    os.makedirs(fitdir)

# if run or read the fit result files
run = True

# MCMC walker parameters
nprocess = 16  # number of MPI process
nwalker = 50
nburnin = 1000
nstep = 1500

# miscelleous parameters for fit
deltav = 0.01  # dv of the resampled axis in km/s; need to be enough small to resolve the expected thermal line width (~0.05 km/s at 10 K)
native_dchan = 122e3  # in Hz; common for all lines
line_vrange = 6  # velocity range to contain the line (+/- w.r.t. vsys); used to select line-free channels

# property of disk
disk = au.DiskProperty(
    PA=disk_dict[source]["PA_gofish"],
    incl=disk_dict[source]["incl"],
    Mstar=disk_dict[source]["M_star"],
    vsys=disk_dict[source]["v_sys"],
    distance=disk_dict[source]["distance"],
)

line_dict = line_dict[source]


#################################################### Useful functions ####################################################
# function to calculate thermal line width
def thermal_line_width(T, mu):
    return np.sqrt(au.k_B * T / (mu * au.m_p)) * 1e-5  # in km/s


############################################## Observtaional data loading ################################################
# shifted, line-only radial spectra
with open(specfilename, "rb") as f:
    data = pickle.load(f)

radii = data["radii"] * disk.distance
spectra = data["spectra"]

# line + continuum spectra without shift
with open(specfilename_wcont, "rb") as f:
    data = pickle.load(f)

assert len(data["radii"]) == len(radii), "Input data are not consistent"
spectra_wcont = data["spectra_wcont"]

# frequency ranges which contain line emission
line_freq_ranges = {}
for line in lines:
    line_freq_ranges[line] = []
    nu0 = float(line_dict[line]["restfreq"].replace("GHz", "")) * 1e9
    freqs = [float(nu.replace("GHz", "")) * 1e9 for nu in line_dict[line]["freqs"]]
    for nu in freqs:
        nu += au.dv2dnu(disk.vsys, nu0=nu0)
        dnu = au.dv2dnu_abs(dv=line_vrange, nu0=nu0)
        line_freq_ranges[line].append((nu - dnu, nu + dnu))


########################################## Spectroscopic data preparation ################################################
freq_ranges = []
nu_list = []
for line in lines:
    spec = spectra[line][0]
    nu_rest = au.unshift_freq(spec.nu, vsys=disk.vsys)
    freq_ranges.append((nu_rest.min(), nu_rest.max()))
    nu_list.append(spec.nu.copy())

transitions = {}
for s in species:
    filename = f"./data/spectroscopy/{s}_hfs.dat"
    specdata = SpectroscopicData()
    specdata.read_LAMDA_file(filename=filename, QNs=["N", "J", "F"])

    # select data
    restfreqs = [trans.restfreq for trans in specdata.radiative.transitions]
    indices = []
    for nurange in freq_ranges:
        numin, numax = nurange
        idx = [i for i, nu in enumerate(restfreqs) if (nu >= numin) and (nu <= numax)]
        indices.extend(idx)
    trans = [
        trans for i, trans in enumerate(specdata.radiative.transitions) if i in indices
    ]
    transitions[s] = trans

specdata = {}
for s in species:
    nu0 = np.array([trans.restfreq for trans in transitions[s]])
    Aul = np.array([trans.EinsteinA for trans in transitions[s]])
    Eu = np.array([trans.Eup for trans in transitions[s]])
    gu = np.array([trans.gup for trans in transitions[s]])

    T, Q = np.loadtxt(f"./data/spectroscopy/{s}_partition_function.dat", unpack=True)
    Q = PartitionFunction(specie=s, T=T, Q=Q)
    specdata[s] = au.SpectroscopicData(specie=s, nu0=nu0, Aul=Aul, Eu=Eu, gu=gu, Q=Q)


def calc_tau(nu, paramset):
    tau_l = 0.0
    for s in species:
        data = specdata[s]

        sigma_nu = au.dv2dnu_abs(
            thermal_line_width(paramset.Tex.value, molecular_weight[s]),
            data.nu0,
        )
        nu_c = data.nu0 + au.dv2dnu(disk.vsys, data.nu0)

        # optical depth
        logN = getattr(paramset, "logN_" + s)
        tau_l += specmodel.calc_line_optical_depth(
            nu=nu,
            nu0=nu_c,
            sigma=sigma_nu,
            Tex=paramset.Tex.value,
            logN=logN.value,
            Aul=data.Aul,
            gu=data.gu,
            Eu=data.Eu,
            Q=data.Q,
        )
    return tau_l


def calc_full_model(paramset):
    model_list = []
    model_cont_list = []
    for line, nu in zip(lines, nu_list):
        tau_l = 0
        dnu = au.dv2dnu_abs(deltav, nu.mean())
        nu_resampled = np.arange(nu.min() - dnu, nu.max() + dnu, dnu)
        tau_l = calc_tau(nu_resampled, paramset)

        # continuum model; linear function of frequency
        I_cont = 10 ** getattr(
            paramset, "log_I_cont_Band" + line_dict[line]["Band"]
        ).value * (nu_resampled / nu_resampled.min())

        # intensity
        I = (
            specmodel.Bnu(nu_resampled, paramset.Tex.value)
            - specmodel.Bnu_CMB(nu_resampled)
            - I_cont
        ) * (1 - np.exp(-tau_l))

        # correct for the beam broadening
        I = au.convolve_Gaussian(
            x=nu_resampled,
            y=I,
            sigma=au.dv2dnu_abs(
                paramset.sigma.value,
                nu_resampled.mean(),
            ),
        )

        # correct for the spectral response function
        I = au.convolve_boxcar(x=nu_resampled, y=I, w=native_dchan)
        I_cont = au.convolve_boxcar(x=nu_resampled, y=I_cont, w=native_dchan)

        # linear interpolation mimicking cvel
        I = interp1d(nu_resampled, I, kind="linear")(nu)
        I_cont = interp1d(nu_resampled, I_cont, kind="linear")(nu)

        # append to list
        model_list.append(I)
        model_cont_list.append(I_cont)
    return model_list, model_cont_list


if __name__ == "__main__":
    ############################################ Parameter Set-up ##############################################################
    param_list = [
        Parameter(
            name="Tex", value=20, bound=(5, 100), free=True, label="$T_\mathrm{ex}$ (K)"
        ),
        Parameter(
            name="log_I_cont_Band3",
            value=8.0,
            bound=(6.0, 12.0),
            free=True,
            label="$\sigma_v$ (km s$^{-1}$)",
        ),
        Parameter(
            name="log_I_cont_Band6",
            value=9.0,
            bound=(6.0, 12.0),
            free=True,
            label="$\sigma_v$ (km s$^{-1}$)",
        ),
        Parameter(
            name="log_I_cont_Band7",
            value=9.0,
            bound=(6.0, 12.0),
            free=True,
            label="$\sigma_v$ (km s$^{-1}$)",
        ),
        Parameter(
            name="sigma",
            value=0.5,
            bound=(0.1, 2),
            free=True,
            label="$\sigma_v$ (km s$^{-1}$)",
        ),
    ]

    param_list += [
        Parameter(
            name="logN_" + specie,
            value=12.0,
            bound=(10.0, 17.0),
            free=True,
            label=f"log$_{{10}}\,N$(" + au.specie_label[specie] + ") (cm$^{{-2}}$) \n",
        )
        for specie in species
    ]

    paramset = ParameterSet(param_list)

    # save the parameter set
    if run:
        filename = fitdir + "parameters.pkl"
        with open(filename, "wb") as f:
            pickle.dump(paramset, f)

    ########################################################## Fit #############################################################
    for i, r in enumerate(radii):
        I_list = []
        dI_list = []
        I_cont_list = []
        dI_cont_list = []
        for line in lines:
            # shifted line spectrum
            spec = spectra[line][i]
            I_list.append(spec.I.copy() / spec.Omega_beam)
            dI_list.append(spec.dI.copy() / spec.Omega_beam)

            # unshifted continuum spectra
            spec = spectra_wcont[line][i]
            I = spec.I.copy()
            dI = spec.dI.copy()

            # replace any possible zero to nan to avoid the overflow at the likelihood calculation
            I[I == 0] = np.nan
            dI[dI == 0] = np.nan
            for nurange in line_freq_ranges[line]:
                I[(spec.nu >= np.min(nurange)) & (spec.nu <= np.max(nurange))] = np.nan
                dI[(spec.nu >= np.min(nurange)) & (spec.nu <= np.max(nurange))] = np.nan

            I_cont_list.append(I / spec.Omega_beam)
            dI_cont_list.append(dI / spec.Omega_beam)

        # MCMC functions

        def log_likelihood(param):
            for name, val in zip(paramset.free_param_name, param):
                getattr(paramset, name).set_value(val)

            model, model_cont = calc_full_model(paramset)

            line_chi2 = np.nansum(
                (
                    (np.concatenate(I_list) - np.concatenate(model))
                    / np.concatenate(dI_list)
                )
                ** 2
            )

            cont_chi2 = np.nansum(
                (
                    (np.concatenate(I_cont_list) - np.concatenate(model_cont))
                    / np.concatenate(dI_cont_list)
                )
                ** 2
            )

            ll = -0.5 * (line_chi2 + cont_chi2)

            return ll

        def log_probability(param):
            lp = log_prior(param, paramset.bound)
            if not np.isfinite(lp):
                return -np.inf
            ll = log_likelihood(param)
            return lp + ll

        if run:
            hammer = EmceeHammer(
                log_probability=log_probability,
                initial_state=paramset.p0_init,
                initial_state_blob_mag=1e-4,
                nwalker=nwalker,
                nstep=nstep,
            )

            print(f"Starting to fit the spectrum at {r} au...")
            with multiprocessing.Pool(processes=nprocess) as pool:
                hammer.run(
                    progress=True,
                    save=True,
                    savefilename=fitdir + "fit.h5",
                    pool=pool,
                    name=f"{int(r)}au",
                )
            print(f"Done.")

        else:
            hammer = EmceeHammer()
            hammer.load_backend(filename=filename, name=f"{int(r)}au")
