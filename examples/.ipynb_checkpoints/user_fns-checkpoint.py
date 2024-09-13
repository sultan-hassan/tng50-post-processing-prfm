from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
def generate_kde(x, y, ax, bin_num, leg_x, co1, ls1, fw, v2, v1):
    """
    Plots a 1D histogram with a kernel density estimation (KDE) overlay.
    
    Parameters:
    -----------
    x : array-like
        Input data for KDE. Can contain NaNs or infinite values.
    
    y : array-like
        Weights for `x`, used only if `fw` is True.
    
    ax : matplotlib.axes.Axes
        Axis object where the plot is drawn.
    
    bin_num : int
        Number of bins for KDE evaluation.
    
    leg_x : str
        Legend label for the plot.
    
    co1 : str or color
        Color for the KDE plot line.
    
    ls1 : str
        Line style for the KDE plot line.
    
    fw : bool
        Flag indicating whether to use weights for KDE.
    
    v2 : float
        Lower bound of the KDE range.
    
    v1 : float
        Upper bound of the KDE range.
    
    Returns:
    --------
    None
    """
    # Remove non-finite values from x (and y if weighted KDE is used)
    mask = np.isfinite(x) & (np.isfinite(y) if fw else True)
    x, y = x[mask], y[mask]

    # Create the kernel density estimation (weighted if fw=True)
    kernel_1 = stats.gaussian_kde(x, weights=y if fw else None)

    # Plot the KDE on the provided axis
    bins = np.linspace(v2, v1, bin_num)
    ax.plot(bins, kernel_1(bins), color=co1, lw=2, ls=ls1, label=leg_x)