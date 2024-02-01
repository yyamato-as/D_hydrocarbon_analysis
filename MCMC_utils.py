import corner
import emcee
from multiprocessing import Pool
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import numpy as np

def setup_params(param_dict):

	param_name = [
		key for key in param_dict.keys() if not param_dict[key]["fixed"]
	]
	fixed_param_name = [
		key for key in param_dict.keys() if param_dict[key]["fixed"]
	]
	bound = [
		param_dict[key]["bound"]
		for key in param_dict.keys()
		if not param_dict[key]["fixed"]
	]
	initial_state = [
		param_dict[key]["p0"]
		for key in param_dict.keys()
		if not param_dict[key]["fixed"]
	]

	return param_name, fixed_param_name, bound, initial_state


class Parameter:
	def __init__(self, name, value, bound=None, free=True, label=None, sample=None):
		self.name = name
		self.value = value
		self.bound = bound
		self.free = free
		self.label = label
		self.sample = sample
		if self.free and (self.bound is None or self.label is None):
			raise ValueError("Provide both bound and label for a free parameter.")
		
	def set_value(self, value):
		self.value = value
	
	def get_value(self):
		return self.value


class ParameterSet:
	def __init__(self, param_list):
		self.param_list = param_list
		for param in self.param_list:
			setattr(self, param.name, param)
		self.param_name = self.get_param_name()
		self.free_param_name = self.get_free_param_name()
		self.fixed_param_name = self.get_fixed_param_name()
		self.bound = self.get_bound()
		self.p0_init = self.get_initial_value()
		self.p0_fixed = self.get_fixed_value()
		self.label = self.get_free_param_label()
	
	def get_values(self):
		return [getattr(self, name) for name in self.param_name]

	def get_param_name(self):
		param_name = [param.name for param in self.param_list]
		return param_name

	def get_free_param_name(self):
		free_param_name = [param.name for param in self.param_list if param.free]
		return free_param_name

	def get_fixed_param_name(self):
		fixed_param_name = [param.name for param in self.param_list if not param.free]
		return fixed_param_name

	def get_bound(self):
		bound = [param.bound for param in self.param_list if param.free]
		return bound

	def get_initial_value(self):
		p0_init = [param.value for param in self.param_list if param.free]
		return p0_init

	def get_fixed_value(self):
		p0_fixed = [param.value for param in self.param_list if not param.free]
		return p0_fixed

	def get_free_param_label(self):
		label = [param.label for param in self.param_list if param.free]
		return label



class EmceeHammer:
	def __init__(
		self,
		params=None,
		log_probability=None,
		initial_state=None,
		nwalker=200,
		nstep=500,
		initial_state_blob_mag=1e-4,
	):

		self.params = params
		self.log_probability = log_probability
		self.initial_state = initial_state
		self.nwalker = nwalker
		self.nstep = nstep

		# set the initial position of each walker
		if self.initial_state is not None:
			self.initial_state = np.atleast_1d(self.initial_state)
			self.ndim = self.initial_state.size
			self.p0 = self.initial_state + initial_state_blob_mag * np.random.randn(self.nwalker, self.ndim)

	def load_backend(self, filename="test.h5", name="mcmc"):
		self.sampler = emcee.backends.HDFBackend(filename, name=name)

	def run(self, pool=None, progress=True, blobs_dtype=None, save=True, savefilename="test.h5", name="mcmc"):
		
		backend = None

		if save:
			backend = emcee.backends.HDFBackend(savefilename, name=name)
			backend.reset(self.nwalker, self.ndim)

		# set sampler
		self.sampler = emcee.EnsembleSampler(
			self.nwalker,
			self.ndim,
			self.log_probability,
			pool=pool,
			blobs_dtype=blobs_dtype,
			backend=backend
		)

		# run
		print(
			"starting to run the MCMC sampling with: \n \t initial state:",
			self.initial_state,
			"\n \t number of walkers:",
			self.nwalker,
			"\n \t number of steps:",
			self.nstep,
		)
		self.sampler.run_mcmc(self.p0, int(self.nstep), progress=progress)

		self.full_sample = self.get_sample(thin=1, nburnin=0)
		if self.params is not None:
			for i, name in enumerate(self.params.free_param_name):
				p = getattr(self.params, name)
				setattr(p, "sample", self.full_sample[:, :, i])


	@property
	def chain_length(self):
		return self.sampler.get_chain().shape[0]

	def get_sample(self, thin=1, nburnin=100, flat=False):
		sample = self.sampler.get_chain(thin=thin, discard=nburnin, flat=flat)
		return sample

	def get_flat_sample(self, thin=1, nburnin=100):
		return self.get_sample(thin=thin, nburnin=nburnin, flat=True)

	def get_blobs(self, thin=1, nburnin=100, flat=False):
		blobs = self.sampler.get_blobs(thin=thin, discard=nburnin, flat=flat)
		return blobs

	def get_flat_blobs(self, thin=1, nburnin=100):
		return self.get_blobs(thin=thin, nburnin=nburnin, flat=True)

	def get_log_prob(self, thin=1, nburnin=100, flat=False):
		return self.sampler.get_log_prob(thin=thin, discard=nburnin, flat=flat)

	def _concat_multiple_sample_arrays(*args, axis=0):
		return np.concatenate(args, axis=axis)

	def get_MAP_params(self, thin=1, nburnin=100, include_blobs=False):
		MAP_params = self.get_flat_sample(thin=thin, nburnin=nburnin)[np.argmax(self.get_log_prob(thin=thin, nburnin=nburnin, flat=True))]
		return MAP_params

	def get_random_sample_params(self, thin=1, nburnin=100, nsample=50):
		# get flat sample
		full_sample = self.get_flat_sample(thin=thin, nburnin=nburnin)

		# generate random number
		rng = np.random.default_rng(20010714)
		rints = rng.integers(low=0, high=full_sample.shape[0], size=nsample)

		params = np.empty((nsample, full_sample.shape[1]))
		for i, n in enumerate(rints):
			params[i] = full_sample[n]
		
		return params



		
	def plot_corner(
		self, thin=1, nburnin=100, include_blobs=False, labels=None, return_fig=False, **kwargs
	):
		sample = self.get_flat_sample(thin=thin, nburnin=nburnin)
		if include_blobs:
			blobs = self.get_flat_blobs(thin=thin, nburnin=nburnin)
			sample = self._concat_multiple_sample_arrays(sample, blobs, axis=1)

		quantiles = kwargs.pop("quantiles", [0.16, 0.5, 0.84])
		show_titles = kwargs.pop("show_titles", True)
		bins = kwargs.pop("bins", 20)

		fig = corner.corner(
			sample,
			labels=labels,
			quantiles=quantiles,
			show_titles=show_titles,
			bins=bins,
			**kwargs
		)
		if return_fig:
			return fig
		
	def plot_walker(self, nburnin=100, labels=None, histogram=True):

		sample = self.get_sample(thin=1, nburnin=0, flat=False).transpose(2, 0, 1)
		# Cycle through the plots.

		figset = []

		for i, s in enumerate(sample):
			fig, ax = plt.subplots()
			for walker in s.T:
				ax.plot(walker, alpha=0.1, color="k")
			ax.set_xlabel("Step number")
			if labels is not None:
				ax.set_ylabel(labels[i])
			if nburnin is not None:
				ax.axvline(nburnin, ls="dotted", color="tab:blue")
			ax.set_xlim(0, s.shape[0])

			# Include the histogram.

			if histogram:
				fig.set_size_inches(
					1.37 * fig.get_figwidth(), fig.get_figheight(), forward=True
				)
				ax_divider = make_axes_locatable(ax)
				bins = np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 50)
				hist, _ = np.histogram(s[nburnin:].flatten(), bins=bins, density=True)
				bins = np.average([bins[1:], bins[:-1]], axis=0)
				ax1 = ax_divider.append_axes("right", size="35%", pad="2%")
				ax1.fill_betweenx(
					bins, hist, np.zeros(bins.size), step="mid", color="darkgray", lw=0.0
				)
				ax1.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1])
				ax1.set_xlim(0, ax1.get_xlim()[1])
				ax1.set_yticklabels([])
				ax1.set_xticklabels([])
				ax1.tick_params(which="both", left=0, bottom=0, top=0, right=0)
				ax1.spines["right"].set_visible(False)
				ax1.spines["bottom"].set_visible(False)
				ax1.spines["top"].set_visible(False)

				# get percentile
				q = np.nanpercentile(s[nburnin:].flatten(), [16, 50, 84])
				for val in q:
					ax1.axhline(val, ls="dashed", color="black")
				text = (
					labels[i]
					+ r"$ = {:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$".format(
						q[1], np.diff(q)[0], np.diff(q)[1]
					)
					if labels is not None
					else r"${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$".format(
						q[1], np.diff(q)[1], np.diff(q)[0]
					)
				)
				ax1.text(0.5, 1.0, text, transform=ax1.transAxes, ha="center", va="top")

			figset.append(fig)

		return figset



def condition(p, b):
	if (p >= b[0]) and (p <= b[1]):
		return True
	return False


def log_prior(param, bound):
	for p, b in zip(param, bound):
		if not condition(p, b):
			return -np.inf
	return 0.0


def emcee_run_wrapper(
	log_probability,
	initial_state,
	nwalker=200,
	nstep=500,
	nburnin=100,
	initial_blob_mag=1e-4,
	get_sample=True,
	get_blobs=False,
	discard=False,
	pool=None,
	blobs_dtype=None,
	nthread=None
):
	# set dimension and initial guesses
	ndim = len(initial_state)
	p0 = initial_state + initial_blob_mag * np.random.randn(nwalker, ndim)

	# set smapler
	sampler = emcee.EnsembleSampler(
		nwalker, ndim, log_probability, pool=pool, blobs_dtype=blobs_dtype, threads=nthread,
	)

	# run
	print(
		"starting to run the MCMC sampling with: \n \t initial state:",
		initial_state,
		"\n \t number of walkers:",
		nwalker,
		"\n \t number of steps:",
		nstep + nburnin,
		"including",
		nburnin,
		"steps as burn in",
	)
	sampler.run_mcmc(p0, int(nstep + nburnin), progress=True)

	if get_sample:
		if discard:
			discard = nburnin
		else:
			discard = 0
		sample = sampler.get_chain(discard=discard, flat=False)
		if get_blobs:
			blobs = sampler.get_blobs(discard=discard, flat=False)
			return sampler, sample, blobs
		return sampler, sample

	return sampler


def set_param(value, bound, fixed, unit):
	return {"value": value, "bound": bound, "fixed": fixed, "unit": unit}


def get_free_param_list(param_dict):
	free_param_list = []
	for param_name in param_dict.keys():
		if "fixed" in param_dict[param_name] and not param_dict[param_name]["fixed"]:
			free_param_list.append(param_name)
		elif "fixed" not in param_dict[param_name]:
			for trans in param_dict[param_name].keys():
				if not param_dict[param_name][trans]["fixed"]:
					free_param_list.append(param_name + "_" + trans)
	return free_param_list


def get_input_param(param, param_dict, free_param_list):
	input_param_dict = {}
	for param_name in param_dict.keys():
		if "value" in param_dict[param_name]:
			if param_name in free_param_list:
				input_param_dict[param_name] = (
					param[free_param_list.index(param_name)]
					* param_dict[param_name]["unit"]
				)
			else:
				input_param_dict[param_name] = (
					param_dict[param_name]["value"] * param_dict[param_name]["unit"]
				)
		else:
			input_param_dict[param_name] = {}
			for trans in param_dict[param_name].keys():
				if param_name + "_" + trans in free_param_list:
					input_param_dict[param_name][trans] = (
						param[free_param_list.index(param_name + "_" + trans)]
						* param_dict[param_name][trans]["unit"]
					)
				else:
					input_param_dict[param_name][trans] = (
						param_dict[param_name][trans]["value"]
						* param_dict[param_name][trans]["unit"]
					)
	return input_param_dict


def get_bound(param_dict):
	bound = []
	for param_name in param_dict.keys():
		if "fixed" in param_dict[param_name] and not param_dict[param_name]["fixed"]:
			bound.append(param_dict[param_name]["bound"])
		elif "fixed" not in param_dict[param_name]:
			for trans in param_dict[param_name].keys():
				if not param_dict[param_name][trans]["fixed"]:
					bound.append(param_dict[param_name][trans]["bound"])
	return bound


def get_initial_guess(param_dict):
	initial_guess = []
	for param_name in param_dict.keys():
		if "fixed" in param_dict[param_name] and not param_dict[param_name]["fixed"]:
			initial_guess.append(param_dict[param_name]["value"])
		elif "fixed" not in param_dict[param_name]:
			for trans in param_dict[param_name].keys():
				if not param_dict[param_name][trans]["fixed"]:
					initial_guess.append(param_dict[param_name][trans]["value"])
	return initial_guess


def plot_corner(
	sample, blob=None, blob_labels=None, labels=None, nburnin=100, **kwargs
):
	sample = sample[nburnin:].reshape(-1, sample.shape[-1])
	if blob is not None:
		if blob.ndim < 3:
			blob = np.array([blob[key] for key in blob.dtype.names]).transpose(
				(1, 2, 0)
			)
		blob = blob[nburnin:].reshape(-1, blob.shape[-1])
		sample = np.concatenate((sample, blob), axis=1)
		if labels is not None and blob_labels is not None:
			labels += blob_labels
	corner_fig = corner.corner(
		sample,
		labels=labels,
		quantiles=[0.16, 0.5, 0.84],
		show_titles=True,
		bins=20,
		**kwargs
	)
	return corner_fig


def plot_chain(sample, labels=None):
	n_param = sample.shape[-1]
	chain_fig, axes = plt.subplots(n_param, sharex=True)
	for i in range(n_param):
		ax = axes[i]
		ax.plot(sample[:, :, i], "k", alpha=0.3)
		ax.set_xlim(0, len(sample))
		if labels is not None:
			ax.set_ylabel(labels[i])
		# ax.yaxis.set_label_coords(-0.1, 0.5)

	axes[-1].set_xlabel("Step number")
	return chain_fig


# from rich teague eddy helper_function.py
def plot_walker(sample, nburnin=None, labels=None, histogram=True):
	#     # Check the length of the label list.

	# 	if labels is not None:
	# 		if sample.shape[0] != len(labels):
	# 			raise ValueError("Incorrect number of labels.")

	sample = sample.transpose((2, 0, 1))
	# Cycle through the plots.

	figset = []

	for i, s in enumerate(sample):
		fig, ax = plt.subplots()
		for walker in s.T:
			ax.plot(walker, alpha=0.1, color="k")
		ax.set_xlabel("Step number")
		if labels is not None:
			ax.set_ylabel(labels[i])
		if nburnin is not None:
			ax.axvline(nburnin, ls="dotted", color="tab:blue")
		ax.set_xlim(0, s.shape[0])

		# Include the histogram.

		if histogram:
			fig.set_size_inches(
				1.37 * fig.get_figwidth(), fig.get_figheight(), forward=True
			)
			ax_divider = make_axes_locatable(ax)
			bins = np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 50)
			hist, _ = np.histogram(s[nburnin:].flatten(), bins=bins, density=True)
			bins = np.average([bins[1:], bins[:-1]], axis=0)
			ax1 = ax_divider.append_axes("right", size="35%", pad="2%")
			ax1.fill_betweenx(
				bins, hist, np.zeros(bins.size), step="mid", color="darkgray", lw=0.0
			)
			ax1.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1])
			ax1.set_xlim(0, ax1.get_xlim()[1])
			ax1.set_yticklabels([])
			ax1.set_xticklabels([])
			ax1.tick_params(which="both", left=0, bottom=0, top=0, right=0)
			ax1.spines["right"].set_visible(False)
			ax1.spines["bottom"].set_visible(False)
			ax1.spines["top"].set_visible(False)

			# get percentile
			q = np.percentile(s[nburnin:].flatten(), [16, 50, 84])
			for val in q:
				ax1.axhline(val, ls="dashed", color="black")
			text = (
				labels[i]
				+ r"$ = {:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$".format(
					q[1], np.diff(q)[0], np.diff(q)[1]
				)
				if labels is not None
				else r"${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$".format(
					q[1], np.diff(q)[1], np.diff(q)[0]
				)
			)
			ax1.text(0.5, 1.0, text, transform=ax1.transAxes, ha="center", va="top")

		figset.append(fig)

	return figset
