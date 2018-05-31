import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
import seaborn as sns
from formatter import FormatParams
from model import Model, calc_sys, calc_astro

sns.set_context("talk")
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})
matplotlib.rcParams.update({'lines.markeredgewidth':0.3})
matplotlib.rcParams.update({'axes.formatter.useoffset':False})

def plot_raw(data):
    palette = sns.color_palette("husl", data.nvisit)
    for i in range(data.nvisit): 	
        ind = data.vis_num==i
        plt.subplot((data.nvisit)*100+10+i+1)
        plt.plot(data.t_vis[ind]*24., data.flux[ind], marker='o', \
                markersize=4.5, linestyle="none", \
                label = "Visit {0}".format(i), color= palette[i])
        plt.xlim(((data.t_vis.min()-0.02)*24., (data.t_vis.max()+0.05)*24.))
        plt.ylim((0.998*data.flux.min(), 1.002*data.flux.max()))
        plt.legend()
    plt.xlabel("Time after visit start (hours)")
    plt.ylabel("Flux (e-)")
    plt.tight_layout()
    plt.show()	


def plot_fit(params, data, flags, obs_par, plot_sys=False):
    p = FormatParams(params, data)
    m = Model(params,data, flags)
    sns.set_palette("muted")
    palette = sns.color_palette("muted", m.data.nvisit)

    #calculate a range of times at higher resolution to make model look nice
    phase_hr = np.linspace(m.phase.min()-0.05, m.phase.max()+0.05, 1000)
    t_hr = phase_hr*p.per[0]+p.t0[0]

    #plot data
    plt.subplot(211)
    #plot best fit model from first visit
    plt.plot(phase_hr, calc_astro(t_hr, params, data, visit = 0))

    #plot normalized data
    for i in range(m.data.nvisit):
        ind = m.data.vis_num == i
        plt.plot(m.phase[ind], m.norm_flux[ind], 'ob')

    #add labels/set axes
    xlo, xhi = np.min(m.phase)*0.9, np.max(m.phase)*1.1
    plt.xlim(xlo,xhi)
    plt.ylabel("Relative Flux")

    #annotate plot with fit diagnostics
    ax = plt.gca()
    ax.text(0.85, 0.29, \
        '$\chi^2_\\nu$:    ' + '{0:0.2f}'.format(int(m.chi2red)) + '\n' \
        + 'obs. rms:  ' + '{0:0d}'.format(int(m.rms)) + '\n' \
        + 'exp. rms:  ' + '{0:0d}'.format(int(m.rms_predicted)), \
            verticalalignment='top',horizontalalignment='left', 
            transform=ax.transAxes, fontsize = 12) 
    
    #plot fit residuals
    plt.subplot(212)
    plt.axhline(0, zorder=1, color='0.2', linestyle='dashed')

    for i in range(m.data.nvisit):
        ind = m.data.vis_num == i
        plt.plot(m.phase[ind], 1.0e6*m.norm_resid[ind], 'ob')

    #add labels/set axes
    plt.xlim(xlo,xhi)
    plt.ylabel("Residuals (ppm)")
    plt.xlabel("Orbital phase")
    
    plt.show()


