import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from scipy import stats

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
    

fw = 1 #flag to weight or not with SFR
nd = 200

line_style    = ['-','--']
xaxis_lims    = np.array([ [0.0,4.0], [-3.0,5.0], [-5.0,3.0], [-2.0,1.0], [-5.0,3.0], [-4.0,2.0]  ])

fig, ax = plt.subplots(2,3,  figsize=(13,7))

ind = 0
for i in [0,2]:
    df = pd.read_pickle('../all_data_z'+str(i)+'.pkl')

    Sigma_g    = np.hstack([np.hstack(x) for x in df['Sigma_g']])*1e-6 # convert from Msun kpc^2 to Msun pc^2
    Sigma_SFR  = np.hstack([np.hstack(x) for x in df['Sigma_SFR']])
    H_star     = np.hstack([np.hstack(x) for x in df['H_*']])
    Sigma_star = np.hstack([np.hstack(x) for x in df['Sigma_*']])*1e-6 # convert from Msun kpc^2 to Msun pc^2
    rho_star   = np.hstack([np.hstack(x) for x in df['rho_*']])*1e-9   # convert from Msun kpc^3 to Msun pc^3
    rho_dm     = np.hstack([np.hstack(x) for x in df['rho_dm']])*1e-9  # convert from Msun kpc^3 to Msun pc^3
    
    generate_kde(np.log10(Sigma_g)    ,Sigma_SFR, ax[0][0], nd, r'$z=$'+str(i),'k', line_style[ind] ,fw, xaxis_lims[0,0]-2,xaxis_lims[0,1]+2)
    generate_kde(np.log10(H_star)     ,Sigma_SFR, ax[1][0], nd, r'$z=$'+str(i),'k', line_style[ind] ,fw, xaxis_lims[1,0]-2,xaxis_lims[1,1]+2)
    generate_kde(np.log10(Sigma_star) ,Sigma_SFR, ax[0][1], nd, r'$z=$'+str(i),'k', line_style[ind] ,fw, xaxis_lims[2,0]-2,xaxis_lims[2,1]+2)
    generate_kde(np.log10(rho_star)   ,Sigma_SFR, ax[1][1], nd, r'$z=$'+str(i),'k', line_style[ind] ,fw, xaxis_lims[4,0]-2,xaxis_lims[4,1]+2)
    generate_kde(np.log10(rho_dm)     ,Sigma_SFR, ax[1][2], nd, r'$z=$'+str(i),'k', line_style[ind] ,fw, xaxis_lims[5,0]-2,xaxis_lims[5,1]+2)
    generate_kde(np.log10(Sigma_SFR)  ,Sigma_SFR, ax[0][2], nd, r'$z=$'+str(i),'k', line_style[ind] ,fw, xaxis_lims[3,0]-2,xaxis_lims[3,1]+2)
    
    ind += 1


label_axis  = [ r'$\log_{10}\, \Sigma_{g}\, [M_{\odot}\, {\rm pc^{-2}}]$', r'$\log_{10}\, \Sigma_{\star}\, [M_{\odot}\, {\rm pc^{-2}}]$', r'$\log_{10} \Sigma_{\rm SFR}\, [M_{\odot}\, {\rm yr^{-1}\, kpc^{-2}}]$', r'$\log_{10}\,  H_{\star}\, [{\rm kpc}]$', r'$\log_{10}\,  \rho_{\star}\, [M_{\odot}\, {\rm pc^{-3}}]$', r'$\log_{10}\,  \rho_{d}\, [M_{\odot}\, {\rm pc^{-3}}]$']


xaxis_ticks = [np.arange(xaxis_lims[0,0],xaxis_lims[0,1]+1,1), np.arange(xaxis_lims[1,0],xaxis_lims[1,1]+2,2),np.arange(xaxis_lims[2,0],xaxis_lims[2,1]+2,2), np.arange(xaxis_lims[3,0],xaxis_lims[3,1]+1,1),
               np.arange(xaxis_lims[4,0],xaxis_lims[4,1]+2,2), np.arange(xaxis_lims[5,0],xaxis_lims[5,1]+2,2)]



ind=0
for j in range(2):
    for k in range(3):
        ax[j,k].tick_params(labelsize=18)  
        ax[j,k].set_xlabel(label_axis[ind], fontsize=15)
        ax[j,k].set_xlim([xaxis_lims[ind][0],xaxis_lims[ind][1]])
        ax[j,k].set_xticks(xaxis_ticks[ind])
        ax[j,k].set_ylim(bottom=0)
        ind+=1

ax[0][0].legend(loc='upper right', fancybox=True,shadow=True, fontsize=14)
plt.subplots_adjust(wspace=0.2, hspace=0.3)
fig.supylabel(r'$ \frac{1}{\rm SFR_{\rm tot}}\frac{\rm d\, SFR}{{\rm d\, log_{10}\,} x} $', fontsize=25, y=0.5, x=0.05)
fig.savefig('./plots/figure_2.pdf', bbox_inches='tight')
