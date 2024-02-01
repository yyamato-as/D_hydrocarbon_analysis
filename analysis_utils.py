import numpy as np
import astropy.units as u
import astropy.constants as ac
import h5py
import matplotlib.pyplot as plt
from qdisk.classes import FitsImage
import casatasks
import casatools
import os
from astropy.modeling.models import Lorentz1D, Gaussian1D
from astropy.convolution import (
    convolve_fft,
    Gaussian1DKernel,
    Box1DKernel,
    Model1DKernel,
)
from astropy.coordinates import SkyCoord
from astropy.io.misc import hdf5
import json

ia = casatools.image()

m_p = ac.m_p.cgs.value
mH2 = 2.016 * ac.u.cgs.value
c = ac.c.cgs.value
k_B = ac.k_B.cgs.value
h = ac.h.cgs.value
ckms = ac.c.to(u.km / u.s).value

arcsec = np.pi / 180.0 / 3600.0  # in rad


class DiskProperty:
    def __init__(self, PA, incl, Mstar, vsys, distance):
        self.PA = PA
        self.incl = incl
        self.Mstar = Mstar
        self.vsys = vsys
        self.distance = distance


line_label = {
    "DCO+_4-3": "DCO$^+$ $J=4$--3",
    "DCN_4-3": "DCN $J=4$--3",
    "C2D_4-3_fs1": "C$_2$D $N=4$--3, $J=9/2$--7/2",
    "C2D_4-3_fs2": "C$_2$D $N=4$--3, $J=7/2$--5/2",
    "C2H_1-0_fs1": "C$_2$H $N=1$--0, $J=3/2$--1/2",
    "C2H_3-2_fs1": "C$_2$H $N=3$--2, $J=7/2$--5/2",
    "C2H_3-2_fs2": "C$_2$H $N=3$--2, $J=5/2$--3/2",
    # "C2H_3-2_hfs1": "C$_2$H $N=4$--3, $J=7/2$--5/2",
}

specie_label = {"C2H": "C$_2$H", "C2D": "C$_2$D"}

color = {
    "C2D_4-3_fs1": "tab:orange",
    "C2D_4-3_fs2": "tab:orange",
    "C2H_1-0_fs1": "lightskyblue",
    "C2H_3-2_fs1": "deepskyblue",
    "C2H_3-2_fs2": "deepskyblue",
}


def decorate_broken_axis(axes, delta=0.1):
    d = 0.5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(
        marker=[(-d, -1), (d, 1)],
        markersize=10,
        linestyle="none",
        color="k",
        mec="k",
        mew=1,
        clip_on=False,
    )
    for i, ax in enumerate(axes):
        # xmin, xmax = ax.get_xlim()
        # ymin, ymax = ax.get_ylim()
        if i != len(axes) - 1:
            ax.spines.right.set_visible(False)
            ax.tick_params(right=False)
            ax.plot(
                [1.0 + delta, 1.0 + delta], [0.0, 1.0], transform=ax.transAxes, **kwargs
            )
        if i != 0:
            ax.spines.left.set_visible(False)
            ax.tick_params(left=False)
            ax.plot(
                [0.0 - delta, 0.0 - delta], [0.0, 1.0], transform=ax.transAxes, **kwargs
            )


def load_radial_profile(filename):
    print("Loading ", filename)
    with h5py.File(filename, "r") as f:
        r = f["radius"][...]
        print("radius unit: ", f["radius"].attrs["unit"])
        I = f["intensity"][...]
        print("intensity unit: ", f["intensity"].attrs["unit"])
        dI = f["scatter"][...]
        assert f["scatter"].attrs["unit"] == f["intensity"].attrs["unit"]
    return r, I, dI


def save_radial_profile(filename, r, I, dI=None, r_unit=None, I_unit=None):
    print("Saving to ", filename)
    with h5py.File(filename, "w") as f:
        dataset = f.create_dataset("radius", data=r)
        dataset.attrs["unit"] = r_unit
        dataset = f.create_dataset("intensity", data=I)
        dataset.attrs["unit"] = I_unit
        dataset = f.create_dataset("scatter", data=dI)
        dataset.attrs["unit"] = I_unit


def load_line_table(tablefilename):
    table = hdf5.read_table_hdf5(tablefilename)
    return table


def save_line_table(table, tablefilename):
    hdf5.write_table_hdf5(table, tablefilename, overwrite=True)


def load_json_file(filename):
    if not os.path.exists(filename):
        print(
            f"No such file: {filename}. Instead will be newly created and return an empty dictionary."
        )
        return dict()
    else:
        with open(filename, "r") as f:
            data = json.load(f)
        return data


def dump_json_file(filename, data):
    if os.path.exists(filename):
        print(f"{filename} already exists. It will be overwritten.")
    with open(filename, "w") as f:
        json.dump(data, f, indent=4)


from matplotlib import colors

cmap = {}
freeze = np.loadtxt("/home/yamato/Project/qDisk/cmap_freeze.txt")
freeze /= 255.0
cmap["M0"] = colors.ListedColormap(freeze, name="freeze")
cmap["M1"] = "RdBu_r"


def flatten(d, parent_key="", sep="_"):
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)


def nest(d, sep="_"):
    nested = {}
    for k, v in d.items():
        context = nested
        for subkey in k.split(sep)[:-1]:
            if subkey not in context:
                context[subkey] = {}
            context = context[subkey]
        context[k.split(sep)[-1]] = v
    return nested


class SpectroscopicData:
    def __init__(self, specie, nu0, Aul, Eu, gu, Q):
        self.specie = specie
        self.nu0 = nu0
        self.Aul = Aul
        self.Eu = Eu
        self.gu = gu
        self.Q = Q


def shift_freq(nu0, vsys):
    return nu0 * (1.0 - vsys / ckms)


def unshift_freq(nu, vsys):
    return nu / (1.0 - vsys / ckms)


def freq2vel(nu, nu0):
    return (1.0 - nu / nu0) * ckms


def vel2freq(v, nu0):
    return (1.0 - v / ckms) * nu0


def dv2dnu(dv, nu0):
    """convert velocity shift to frequency shift with a given reference frequency nu0

    Parameters
    ----------
    dv : float or ndarray
        velocity shift(s) in km/s
    nu0 : float or ndarray
        reference frequency(s) in Hz

    Returns
    -------
    float or ndarray
        frequency shift in Hz
    """
    return -dv / ckms * nu0


def dv2dnu_abs(dv, nu0):
    """convert velocity width to frequency width (absolute value) with a given reference frequency nu0

    Parameters
    ----------
    dv : float or ndarray
        velocity width(s) in km/s
    nu0 : float or ndarray
        reference frequency(s) in Hz

    Returns
    -------
    float or ndarray
        frequency width in Hz
    """
    return np.abs(-dv / ckms * nu0)


def dnu2dv(dnu, nu0):
    """convert frequency shift to velocity shift with a given reference frequency nu0

    Parameters
    ----------
    dnu : float or ndarray
        frequency shift(s) in Hz
    nu0 : float or ndarray
        reference frequency(s) in Hz

    Returns
    -------
    float or ndarray
        frequency shift in Hz
    """
    return -dnu / nu0 * ckms


def dnu2dv_abs(dnu, nu0):
    """convert frequency width to velocity width (absolute value) with a given reference frequency nu0

    Parameters
    ----------
    dnu : float or ndarray
        frequency shift(s) in Hz
    nu0 : float or ndarray
        reference frequency(s) in Hz

    Returns
    -------
    float or ndarray
        frequency shift in Hz
    """
    return np.abs(-dnu / nu0 * ckms)


def subimage(
    fitsimage,
    outfile,
    chans=None,
    nu0=None,
    vsys=0.0,
    vlim=None,
    nulim=None,
    xlim=None,
    ylim=None,
):
    if chans is None:
        imagecube = FitsImage(fitsimage, skipdata=True)
        if vlim is not None and nulim is not None:
            print(
                "WARNING: Both vlim and nulim is provided. Going to use nulim over vlim to select the channel..."
            )
        if nulim is not None:
            numin, numax = nulim
            if numin > numax:
                raise ValueError(
                    "nulim must be a 2-component tuple whose first component is smaller than the second."
                )
            elif numin == numax:
                index = np.argmin(np.abs(imagecube.nu - numin))
                chans = f"{index}~{index}"
            else:
                indices = np.where((imagecube.nu >= numin) & (imagecube.nu <= numax))
                chans = f"{str(np.min(indices))}~{str(np.max(indices))}"
        elif vlim is not None:
            vmin, vmax = vlim
            vmin += vsys
            vmax += vsys
            if vmin > vmax:
                raise ValueError(
                    "vlim must be a 2-component tuple whose first component is smaller than the second."
                )
            if nu0 is not None:
                if isinstance(nu0, (list, np.ndarray)):
                    mask = np.zeros(imagecube.nu.shape, dtype=bool)
                    for nu in nu0:
                        v = (1 - imagecube.nu / nu) * ckms
                        m = (v >= vmin) & (v <= vmax)
                        mask = mask | m
                else:
                    v = (1.0 - imagecube.nu / nu0) * ckms
                    mask = (v >= vmin) & (v <= vmax)
                indices = np.where(mask)
                chans = f"{str(np.min(indices))}~{str(np.max(indices))}"
            else:
                if vmin == vmax:
                    index = np.argmin(np.abs(imagecube.v - vmin))
                    chans = f"{index}~{index}"
                else:
                    indices = np.where((imagecube.v >= vmin) & (imagecube.v <= vmax))
                    chans = f"{str(np.min(indices))}~{str(np.max(indices))}"
        else:
            chans = ""
    if xlim is not None:
        _xmin, _xmax = xlim
        if _xmin >= _xmax:
            raise ValueError(
                "xlim must be a 2-component tuple whose first component is smaller than the second."
            )
        xmax = np.argmin(abs(imagecube.x - _xmin))
        xmin = np.argmin(abs(imagecube.x - _xmax))  # note that x is in decending order
    else:
        xmin = 0
        xmax = imagecube.x.size - 1
    if ylim is not None:
        _ymin, _ymax = ylim
        if _ymin >= _ymax:
            raise ValueError(
                "ylim must be a 2-component tuple whose first component is smaller than the second."
            )
        ymin = np.argmin(abs(imagecube.y - _ymin))
        ymax = np.argmin(abs(imagecube.y - _ymax))
    else:
        ymin = 0
        ymax = imagecube.x.size - 1
    region = f"box[[{xmin}pix, {ymin}pix], [{xmax}pix, {ymax}pix]]"

    tmpimage = "tmp.image"
    casatasks.imsubimage(
        imagename=fitsimage,
        outfile=tmpimage,
        chans=chans,
        region=region,
        overwrite=True,
    )
    if nu0 is not None:
        if isinstance(nu0, (list, np.ndarray)):
            nu0 = nu0[int(len(nu0) / 2)]
        ia.open(tmpimage)
        csys = ia.coordsys()
        csys.setrestfrequency(value=str(nu0) + "Hz")
        ia.setcoordsys(csys.torecord())
        ia.close()
    casatasks.exportfits(
        imagename=tmpimage, fitsimage=outfile, overwrite=True, dropstokes=True
    )
    os.system("rm -r " + tmpimage)
    return


class Spectrum:
    def __init__(
        self,
        filename=None,
        nu=None,
        v=None,
        channel=None,
        I=None,
        dI=None,
        unit=None,
        nu0=None,
        beam=None,
        **kwargs,
    ):
        """_summary_

        Parameters
        ----------
        filename : str, optional
            hdf5 filename which contains the spectrum data, by default None
        nu : array_like, optional
            frequency axis in Hz, by default None
        v : array_like, optional
            velocity axis in km/s, by default None
        channel : array_like, optional
            channel index, by default None
        I : array_like, optional
            spectrum in an arbitrary unit, by default None
        dI : array_like, optional
            uncertainty of the spectrum in the same unit as that of spectrum, by default None
        nu0 : float, optional
            rest frequency (or representative frequency of the spectrum) in Hz, by default None
        beam : tuple, optional
            beam in (bmaj, bmin, bpa) in arcsec and deg, by default None
        """
        self.filename = filename
        self._nu = nu
        self._v = v
        self._channel = channel
        self._I = I
        self._dI = dI
        self._unit = unit
        self._nu0 = nu0
        self._beam = beam
        self._kwargs = kwargs
        for key in kwargs:
            setattr(self, "_" + key, kwargs[key])

        self.restore()

    def restore(self):
        if self.filename is not None:
            with h5py.File(self.filename, "r") as f:
                # spectrum data
                data = f["data"]
                self.v = data["velocity"][...] if "velocity" in data.keys() else None
                self.nu = data["frequency"][...] if "frequency" in data.keys() else None
                self.channel = (
                    data["channel"][...] if "channel" in data.keys() else None
                )
                self.I = data["spectrum"][...] if "spectrum" in data.keys() else None
                self.dI = data["scatter"][...] if "scatter" in data.keys() else None

                # metadata
                self.unit = data["spectrum"].attrs["unit"]
                md = f["metadata"]
                for key in md.attrs.keys():
                    setattr(self, key, md.attrs[key])
                # self.nu0 = md.attrs["nu0"]
                # self.beam = md.attrs["beam"]
        else:
            self.v = (
                (1 - self._nu / self._nu0) * ckms
                if (self._nu is not None) and (self._nu0 is not None)
                else self._v
            )
            self.nu = (
                (1 - self._v / ckms) * self._nu0
                if (self._v is not None) and (self._nu0 is not None)
                else self._nu
            )
            self.channel = (
                self._channel if self._channel is not None else np.arange(self._I.size)
            )
            self.I = self._I
            self.dI = self._dI
            self.unit = self._unit
            self.nu0 = self._nu0
            self.beam = self._beam
            for key in self._kwargs.keys():
                setattr(self, key, getattr(self, "_" + key))

    def _set_velocity_axis(self, nu0=None):
        if nu0 is None:
            pass
        else:
            self.v = (1 - self.nu / nu0) * ckms

    def split(self, nu0=None, vrange=None, nurange=None, vsys=0.0):
        if vrange is not None:
            vstart, vend = vrange
            vstart += vsys
            vend += vsys

            if nu0 is not None and isinstance(nu0, (list, np.ndarray)):
                mask = np.zeros(self.nu.shape, dtype=bool)
                for nu in nu0:
                    v = (1 - self.nu / nu) * ckms
                    m = (v >= vstart) & (v <= vend)
                    mask = mask | m
                self.nu0 = nu0[int(len(nu0) / 2)]
                self.v = (1 - self.nu / self.nu0) * ckms
            else:
                self.nu0 = nu0 if nu0 is not None else self.nu0
                self.v = (1 - self.nu / self.nu0) * ckms
                mask = (self.v >= vstart) & (self.v <= vend)

        elif nurange is not None:
            nustart, nuend = nurange
            mask = (self.nu >= nustart) & (self.nu <= nuend)

        else:
            mask = np.ones_like(self.I, dtype=bool)

        self.v = self.v[mask] if self.v is not None else self.v
        self.nu = self.nu[mask] if self.nu is not None else self.nu
        self.channel = self.channel[mask] if self.channel is not None else self.channel
        self.I = self.I[mask] if self.I is not None else self.I
        self.dI = self.dI[mask] if self.dI is not None else self.dI

    def plot(
        self,
        ax=None,
        axis="velocity",
        scaling=1.0,
        error_scaling=1.0,
        nu0=None,
        vrange=None,
        nurange=None,
        vsys=0.0,
        indicate_loc=False,
        **errorbar_kwargs,
    ):
        self.restore()
        self.split(nu0=nu0, vrange=vrange, nurange=nurange, vsys=vsys)
        self._set_velocity_axis(nu0=np.nanmedian(nu0) if nu0 is not None else self.nu0)

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        if "freq" in axis:
            x = self.nu * 1e-9
            if indicate_loc:
                dnu = (1.0 - vsys / ckms) * self.nu0
                x0 = (np.array(nu0) + dnu) * 1e-9
        elif "chan" in axis:
            x = self.channel
            if indicate_loc:
                raise NotImplementedError
        else:
            x = self.v
            if indicate_loc:
                dv = (1.0 - np.array(nu0) / self.nu0) * ckms
                x0 = vsys + dv

        y = self.I * scaling
        yerr = self.dI * scaling * error_scaling if self.dI is not None else self.dI

        self.errorbar = ax.errorbar(x, y, yerr=yerr, **errorbar_kwargs)
        ax.set(xlim=(x.min(), x.max()))

        if indicate_loc:
            for loc in x0:
                ax.axvline(x=loc, color="grey", lw=0.8, alpha=0.5, ls="dotted")

        return fig, ax

    def save_to_hdf5(self, savefilename):
        self.quantities_to_save = ["v", "nu", "channel", "I", "dI"]
        self.quantities_names = [
            "velocity",
            "frequency",
            "channel",
            "spectrum",
            "scatter",
        ]
        self.quantities_units = ["km/s", "Hz", "", self.unit, self.unit]
        with h5py.File(savefilename, "w") as f:
            g = f.create_group("data")
            for quantity, name, unit in zip(
                self.quantities_to_save, self.quantities_names, self.quantities_units
            ):
                q = getattr(self, quantity)
                if q is not None:
                    dataset = g.create_dataset(name, data=q)
                    dataset.attrs["unit"] = unit if unit is not None else ""
                else:
                    self.quantities_to_save.remove(quantity)
                    self.quantities_names.remove(name)
                    self.quantities_units.remove(unit)

            # dataset = g.create_dataset("velocity", data=self.v)
            # dataset.attrs["unit"] = "km/s"
            # dataset = g.create_dataset("frequency", data=self.nu)
            # dataset.attrs["unit"] = "Hz"
            # dataset = g.create_dataset("channel", data=self.channel)
            # dataset.attrs["unit"] = ""
            # dataset = g.create_dataset("spectrum", data=self.I)
            # dataset.attrs["unit"] = self.unit
            # dataset = g.create_dataset("scatter", data=self.dI)
            # dataset.attrs["unit"] = self.unit

            g_meta = f.create_group("metadata")
            if self.beam is not None:
                g_meta.attrs["beam"] = self.beam
            if self.nu0 is not None:
                g_meta.attrs["nu0"] = self.nu0
            for key in self._kwargs.keys():
                if getattr(self, key) is not None:
                    g_meta.attrs[key] = getattr(self, key)


def get_line_cluster(restfreqs, dv=10):
    array = np.sort(restfreqs)
    breakpoints = []
    seed = array[0]
    for i in range(array.size):
        dnu = dv2dnu_abs(dv, seed)
        if (array[i] - seed) > dnu:
            breakpoints.append(i)
        seed = array[i]
    return np.split(array, breakpoints)


# group the restfrequency
def group_frequencies_old(table, dv=10):
    freq_list = []
    spw_list = []
    nu_prev = table["Rest Frequency (GHz)"][0]
    spw_prev = table["spw"][0]
    nu_list = []
    for i in range(len(table)):
        nu = table["Rest Frequency (GHz)"][i]
        spw = table["spw"][i]
        dnu = np.abs(dv / ckms * nu)
        if np.abs(nu - nu_prev) <= dnu:
            nu_list.append(nu * 1e9)
            nu_prev = nu
            spw_prev = spw
        else:
            freq_list.append(nu_list)
            spw_list.append(spw_prev)
            nu_list = [nu * 1e9]
            nu_prev = nu
            spw_prev = spw
        if i == len(table) - 1:
            freq_list.append(nu_list)
            spw_list.append(spw)
    return spw_list, freq_list


def multiplot(
    npanel=None,
    ncols=None,
    nrows=None,
    figsize=(4, 3),
    xlabel=None,
    ylabel=None,
    sharex=False,
    sharey=False,
    max_figsize=(15, None),
):
    if (npanel is not None) and (ncols is None) and (nrows is None):
        nrows = int(npanel**0.5)
        ncols = int(npanel / nrows * 0.999) + 1

        _width, _height = figsize
        width = ncols * _width
        height = nrows * _height

        width_max, height_max = max_figsize

        if (width_max is not None) and width > width_max:
            ncols = int(width_max / _width)
            nrows = int(npanel / ncols * 0.999) + 1
        if (height_max is not None) and (height > height_max):
            nrows = int(height_max / _height)
            ncols = int(npanel / nrows * 0.999) + 1

    # ncols = 5
    # nrows = int(npanel/ncols*0.999) + 1
    fig, axes = plt.subplots(
        figsize=(ncols * figsize[0], nrows * figsize[1]),
        nrows=nrows,
        ncols=ncols,
        sharex=sharex,
        sharey=sharey,
        constrained_layout=True,
    )
    axes.flatten()[int((nrows - 1) * ncols)].set(xlabel=xlabel, ylabel=ylabel)

    if npanel is not None:
        for j in range(npanel, len(axes.flatten())):
            ax = axes.flatten()[j]
            ax.set_axis_off()

    return fig, axes


def spatially_integrate(spectrum, source_size, beam):
    FWHM = np.sqrt(beam**2 + source_size**2)
    Omega = get_beam_solid_angle((FWHM, FWHM))
    spectrum *= Omega
    return spectrum


def convolve_Lorentzian(x, y, gamma):
    dx = np.diff(x).mean()
    kernel = Model1DKernel(
        Lorentz1D(amplitude=1.0 / (np.pi * gamma), x_0=0.0, fwhm=2 * gamma / dx),
        x_size=x.size,
    )
    kernel = Model1DKernel(
        Lorentz1D(amplitude=1.0, x_0=0.0, fwhm=2 * gamma / dx), x_size=x.size
    )
    model = convolve_fft(y, kernel, normalize_kernel=True)
    return model


def convolve_Gaussian(x, y, sigma):
    dx = np.diff(x).mean()
    kernel = Gaussian1DKernel(stddev=sigma / dx, x_size=x.size)
    convolved = convolve_fft(y, kernel, normalize_kernel=True)
    return convolved


def convolve_boxcar(x, y, w):
    dx = np.diff(x).mean()
    kernel = Box1DKernel(w / dx)
    convolved = convolve_fft(y, kernel, normalize_kernel=True)
    return convolved


### deprecated
def load_spectrum(filename, xaxis="frequency"):
    with h5py.File(filename, "r") as f:
        dataset = f["spectrum"]
        if xaxis == "velocity":
            x = dataset["velocity"][...]
        else:
            x = dataset["frequency"][...]
        try:
            I = dataset["averaged_spectrum"][...]
        except KeyError:
            I = dataset["integrated_spectrum"][...]
        dI = dataset["scatter"][...]
    return x, I, dI


def get_beam_solid_angle(beam):
    """Calculate the beam solid angle.

    Parameters
    ----------
    beam : tuple
        A tuple of beam size, i.e., (bmaj, bmin) in arcsec.

    Returns
    -------
    float
        Beam solid angle in steradian.
    """
    return np.multiply(*beam) * arcsec**2 * np.pi / (4 * np.log(2))


def jypb_to_jypsr(I, beam):
    """Convert intensity in Jy / beam to Jy / sr.

    Parameters
    ----------
    I : float or array_like
        Intensity in Jy / beam.
    beam : tuple
        A tuple of beam size, i.e., (bmaj, bmin) in arcsec.

    Returns
    -------
    float or array_like
        Intensity in Jy / sr.
    """
    Omega_beam = get_beam_solid_angle(beam)
    return I / Omega_beam


def jypsr_to_jypb(I, beam):
    """Convert intensity in Jy / beam to Jy / sr.

    Parameters
    ----------
    I : float or array_like
        Intensity in Jy / beam.
    beam : tuple
        A tuple of beam size, i.e., (bmaj, bmin) in arcsec.

    Returns
    -------
    float or array_like
        Intensity in Jy / sr.
    """
    Omega_beam = get_beam_solid_angle(beam)
    return I * Omega_beam


def jypsr_to_K_RJ(I, nu):
    """Convert intensity in Jy / sr to K using RJ approx.

    Parameters
    ----------
    I : float or array_like
        Intensity in Jy / sr.
    nu : float or array_like
        Observing frequency in Hz.

    Returns
    -------
    float or array_like
        Intensity in K in R-J approx.
    """
    I *= 1e-23  # in cgs
    return c**2 / (2 * k_B * nu**2) * I


def cgs_to_jypb(I, beam):
    """Convert intensity in Jy / beam to cgs unit. Jy = 1e-23 erg/s/cm2/Hz

    Parameters
    ----------
    I : float or array_like
        Intensity in Jy / beam.
    beam : tuple
        A tuple of beam size, i.e., (bmaj, bmin) in arcsec.

    Returns
    -------
    float or array_like
        Intensity in erg s-1 cm-2 Hz-1 sr-1.
    """
    return jypsr_to_jypb(I, beam) / 1e-23


def jypb_to_cgs(I, beam):
    """Convert intensity in Jy / beam to cgs unit. Jy = 1e-23 erg/s/cm2/Hz

    Parameters
    ----------
    I : float or array_like
        Intensity in Jy / beam.
    beam : tuple
        A tuple of beam size, i.e., (bmaj, bmin) in arcsec.

    Returns
    -------
    float or array_like
        Intensity in erg s-1 cm-2 Hz-1 sr-1.
    """
    return jypb_to_jypsr(I, beam) * 1e-23


def jypb_to_K_RJ(I, nu, beam):
    """Convert intensity in Jy / beam to birghtness temeprature in Kelvin using RJ approximation.

    Parameters
    ----------
    I : float or array_like
        Intensity in Jy / beam.
    nu : float or array_like
        Observing frequency in Hz.
    beam : tuple
        A tuple of beam size, i.e., (bmaj, bmin) in arcsec.

    Returns
    -------
    float or array_like
        Brightness temperature in RJ approximation.
    """
    I = jypb_to_cgs(I, beam)
    return c**2 / (2 * k_B * nu**2) * I


def jypb_to_K(I, nu, beam):
    """Convert intensity in Jy /beam to brightness temperature in Kelvin using full planck function.

    Parameters
    ----------
    I : float or array_like
        Intenisty in Jy /beam.
    nu : flaot or array_like
        Observing frequency in Hz.
    beam : tuple
        A tuple of beam size, i.e., (bmaj, bmin) in arcsec.

    Returns
    -------
    float or array_like
        Brightness temperature.
    """
    T = np.abs(jypb_to_cgs(I, beam))
    T = h * nu / k_B / np.log(1 + 2 * h * nu**3 / (c**2 * T))

    if isinstance(I, np.ndarray):
        return np.where(I >= 0.0, T, -T)
    elif isinstance(I, float):
        return T if I >= 0.0 else -T


def K_to_jypb(T, nu, beam):
    """Convert intensity in Jy /beam to brightness temperature in Kelvin using full planck function.

    Parameters
    ----------
    I : float or array_like
        Intenisty in Jy /beam.
    nu : flaot or array_like
        Observing frequency in Hz.
    beam : tuple
        A tuple of beam size, i.e., (bmaj, bmin) in arcsec.

    Returns
    -------
    float or array_like
        Brightness temperature.
    """
    I = 2 * h * nu**3 / c**2 / (np.exp(h * nu / (k_B * np.abs(T))) - 1)
    I = cgs_to_jypb(I, beam)

    if isinstance(T, np.ndarray):
        return np.where(I >= 0.0, I, -I)
    elif isinstance(T, float):
        return I if T >= 0.0 else -I


def K_to_jypb_RJ(T, nu, beam):
    """Convert intensity in Jy /beam to brightness temperature in Kelvin using full planck function.

    Parameters
    ----------
    I : float or array_like
        Intenisty in Jy /beam.
    nu : flaot or array_like
        Observing frequency in Hz.
    beam : tuple
        A tuple of beam size, i.e., (bmaj, bmin) in arcsec.

    Returns
    -------
    float or array_like
        Brightness temperature.
    """
    I = 2 * k_B * nu**2 / c**2 * T
    I = cgs_to_jypb(I, beam)

    return I


def jypb_to_K_astropy(I, nu, beam):
    """Convert intensity in Jy /beam to brightness temperature in Kelvin using RJ approximation implemented in astropy.

    Parameters
    ----------
    I : float or array_like
        Intenisty in Jy /beam.
    nu : flaot or array_like
        Observing frequency in Hz.
    beam : tuple
        A tuple of beam size, i.e., (bmaj, bmin) in arcsec.

    Returns
    -------
    float or array_like
        Brightness temperature.
    """
    I *= u.Jy / u.beam
    nu *= u.Hz
    Omega_beam = np.multiply(*beam) * u.arcsec**2 * np.pi / (4 * np.log(2))
    return I.to(
        u.K, equivalencies=u.brightness_temperature(nu, beam_area=Omega_beam)
    ).value


def sigma_to_FWHM(sigma):
    return sigma * np.sqrt(8 * np.log(2))


def FWHM_to_sigma(FWHM):
    return FWHM / np.sqrt(8 * np.log(2))


def populate_imfit_result(imfitdict):
    if not imfitdict["converged"]:
        raise ValueError("The fit has not converged.")


def get_casa_ellipse_region(dir, maj, min, pa):
    ra = dir.split()[0].replace("h", ":").replace("m", ":").replace("s", "")
    dec = dir.split()[1].replace("d", ".").replace("m", ".").replace("s", "")

    region = "ellipse[[%s, %s], [%.1farcsec, %.1farcsec], %.1fdeg]" % (
        ra,
        dec,
        maj,
        min,
        pa,
    )

    return region


def get_source_position(imfitdict, deconvolved=True, unit=None, frame="icrs"):
    result = imfitdict["deconvolved"] if deconvolved else imfitdict["results"]
    ncomp = result["nelements"]

    position = {}
    for i in range(ncomp):
        ra = result["component0"]["shape"]["direction"]["m0"]["value"]
        ra_unit = u.Unit(result["component0"]["shape"]["direction"]["m0"]["unit"])
        dec = result["component0"]["shape"]["direction"]["m1"]["value"]
        dec_unit = u.Unit(result["component0"]["shape"]["direction"]["m1"]["unit"])
        refer = result["component0"]["shape"]["direction"]["refer"]

        if frame == "J2000":
            c = SkyCoord(
                ra=ra * ra_unit, dec=dec * dec_unit, frame="fk5", equinox=refer
            )
        elif frame == "icrs":
            c = SkyCoord(ra=ra * ra_unit, dec=dec * dec_unit, frame="icrs")

        if unit is None:
            position["component{:d}".format(i)] = c.to_string("hmsdms")
        else:
            position["component{:d}".format(i)] = {
                "ra": c.ra.to(u.Unit(unit)),
                "dec": c.dec.to(u.Unit(unit)),
            }

    return position


def get_source_size(imfitdict, deconvolved=True, unit="arcsec"):
    result = imfitdict["deconvolved"] if deconvolved else imfitdict["results"]
    ncomp = result["nelements"]

    source_size = {}
    for i in range(ncomp):
        resultdict = result["component{:d}".format(i)]["shape"]

        source_size["component{:d}".format(i)] = {}
        for key in ["majoraxis", "minoraxis", "positionangle"]:
            unit = "deg" if key == "positionangle" else unit
            val = (
                (resultdict[key]["value"] * u.Unit(resultdict[key]["unit"]))
                .to(unit)
                .value
            )
            err = (
                (
                    resultdict[key + "error"]["value"]
                    * u.Unit(resultdict[key + "error"]["unit"])
                )
                .to(unit)
                .value
            )
            source_size["component{:d}".format(i)][key] = (val, err)

    return source_size
