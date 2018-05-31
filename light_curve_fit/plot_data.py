import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
import seaborn as sns
from formatter import FormatParams
from light_curve_model import Model, get_elc, get_tlc, get_phaselc

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
        plt.plot(data.t_vis[ind]*24., data.flux[ind], marker='o', markersize=4.5, linestyle="none", label = "Visit {0}".format(i), color= palette[i])
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
    if plot_sys==False:
        for i in range(m.data.nvisit):
            ind = m.data.vis_num == i
            plt.subplot(211)
            if data.lc_type == "transit":
                temp_phase = np.copy(m.phase[ind])
                temp_phase[temp_phase>0.5] -= 1.0
                plt.plot(temp_phase, m.data_corr[ind], marker='o', linestyle="None", markersize=5)
            else: plt.plot(m.phase[ind], m.data_corr[ind], marker='o', linestyle="None", markersize=5)
            plt.subplot(212)
            if data.lc_type == "transit": plt.plot(temp_phase, 1.0e6*m.resid[ind]/m.data.flux[ind], marker='o', linestyle="None", markersize=5)
            else: plt.plot(m.phase[ind], 1.0e6*m.resid[ind]/m.data.flux[ind], marker='o', linestyle="None", markersize=5)
        plt.subplot(211)
        phase_hr = np.linspace(m.phase.min()-0.05, m.phase.max()+0.05, 1000)
        t_hr = phase_hr*p.per[0]+p.t0[0]

        if data.lc_type == "transit":
            transit_model_hr = get_tlc(t_hr, p, data, 0) 
            lc_hr = transit_model_hr
        elif data.lc_type == "eclipse":
            eclipse_model_hr = get_elc(t_hr, p, data, 0) 
            lc_hr = eclipse_model_hr
        elif data.lc_type == "phase_curve":
            transit_model_hr = get_tlc(t_hr, p, data, 0) 
            eclipse_model_hr = get_elc(t_hr, p, data, 0) 
            phase_model_hr = get_phaselc(t_hr, p, data, 0)
            lc_hr = transit_model_hr + (eclipse_model_hr-1.)*phase_model_hr
        else: assert False, "Unknown option; supported light curve types are 'transit', 'eclipse', and 'phase_curve'"
        
        if data.lc_type == "transit": 
            phase_hr[phase_hr > 0.5] -= 1.0
            ind = np.argsort(phase_hr)
            phase_hr = phase_hr[ind]
            lc_hr = lc_hr[ind]
                
        plt.plot(phase_hr, lc_hr, color='0.2', zorder=1)
        plt.ylabel("Relative flux")
        ax = plt.gca()
        ax.text(0,1,'obs, exp rms (ppm); chi:\n '+'{0:d}'.format(int(m.rms))+", "+'{0:d}'.format(int(m.rms_predicted))+", "+'{0:0.2f}'.format(m.chi2red),verticalalignment='top',horizontalalignment='left',transform=ax.transAxes, fontsize=10) #LK comment!!
        if flags['fit-white']: plt.title("Fit to white light curve")
        else: plt.title("Fit to {0:0.2f}".format(m.data.wavelength)+" micron channel")
        delta = 3.0*m.rms/1.0e6
        plt.ylim((lc_hr.min() - delta , lc_hr.max()+delta))
        if data.lc_type == "phase_curve": plt.ylim((1.0-delta, lc_hr.max()+delta))
        plt.xlim((phase_hr.min(), phase_hr.max()))
        if obs_par['lc_type'] == 'transit': plt.xlim(-0.03, 0.03)
        ax = plt.gca()
        ax.ticklabel_format(useOffset=False)

        plt.subplot(212)
        plt.axhline(0, zorder=1, color='0.2', linestyle='dashed')
        plt.ylabel("Residuals (ppm)")
        plt.xlabel("Orbital phase")
        plt.xlim((phase_hr.min(), phase_hr.max()))
        if obs_par['lc_type'] == 'transit': plt.xlim(-0.03, 0.03)
        plt.show()	#LK comment
        #plt.savefig("temp.pdf")

        """plt.hist(m.resid/m.data.flux*1e6, normed=True)

        xx = np.linspace(-2000, 2000, 100)
        plt.plot(xx, norm.pdf(xx, 0, int(m.rms_predicted)))

        plt.show()"""
        #make rms vs. bin size plot

    elif plot_sys == True:
        #offset = np.mean(m.sys)*0.001
        print "FIXME: this option isn't implemented for phase curves"
        offset = 0.
        for i in range(m.data.nvisit):
            ind = m.data.vis_num == i
            plt.subplot(211)
            plt.plot(m.phase[ind], m.data_normalized[ind]+offset*i, marker='o', linestyle="None", markersize=5, color=palette[i])
            plt.plot(m.phase[ind], m.sys[ind]*m.transit_model[ind]+offset*i, color=palette[i], zorder=1)
            plt.plot(m.phase[ind], m.sys[ind]+offset*i, color=palette[i], zorder=1, linestyle='dashed')
            plt.ylabel("Corrected flux")
            plt.title("Flux corrected for orbit-long trends and scan direction")
            plt.subplot(212)
            plt.plot(m.phase[ind], m.resid[ind], marker='o', linestyle="None", markersize=5, color=palette[i])
            plt.xlabel("Orbital phase")
            plt.ylabel("Residuals (e-)")
        plt.subplot(211)
        plt.plot(m.phase[0], m.sys[0]*m.transit_model[0], color='0.5', zorder=1, label="with transit")
        plt.plot(m.phase[0], m.sys[0], color='0.5', zorder=1, linestyle='dashed', label="transit removed")
        plt.legend()
        plt.show()
                
