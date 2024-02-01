import matplotlib.pyplot as plt


def set_plot_style(style="paper"):
    """Set the plot style including appearance etc.

    Parameters
    ----------
    style : str, optional
        plot style, "paper", "slides", or "japanese", by default "paper"
    """

    if style == "paper":
        plt.rcParams.update(
            {
                "text.usetex": True,
                "font.family": "serif",
                # "font.sans-serif": ["Times"],
                "xtick.top": True,
                "ytick.right": True,
                "xtick.direction": "out",
                "ytick.direction": "out",
                "xtick.color": "black",
                "xtick.labelcolor": "black",
                "xtick.labelsize": 8,
                "xtick.major.size": 3.0,
                "ytick.color": "black",
                "ytick.labelcolor": "black",
                "ytick.labelsize": 8,
                "ytick.major.size": 3.0,
            }
        )