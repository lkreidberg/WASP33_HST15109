import numpy as np
from formatter import FormatParams

class Model:
    """
    doc
    """
    def __init__(self, params, data, flags):
        p = FormatParams(params, data)
        #full model light curve (with systematics)
        self.lc = np.zeros(len(data.time))				
        #full model light curve (no systematics)
        self.lc_nosys = np.zeros(len(data.time))			
        self.transit_model = np.zeros(len(data.time))			#transit model (relative flux; no systematics)
        self.eclipse_model = np.zeros(len(data.time))			#eclipse model (relative flux; no systematics)
        self.phase_model = np.zeros(len(data.time))			#phase curve model (relative flux; no systematics; includes eclipse)
        self.phase_model_no_eclipse = np.zeros(len(data.time))		#phase variation model (relative flux; no systematics; no eclipse))
        self.data_corr = np.zeros(len(data.time))			#data with the odd/even effect and orbit-long ramps removed
        self.phase = np.zeros(len(data.time))				#orbital phase	(defined in model because it depends on the period and ephemeris)
        self.sys = np.zeros(len(data.time))				#systematics showing visit-long trends only
        self.all_sys = np.zeros(len(data.time))				#all systematics 
        self.data_normalized = np.zeros(len(data.time))

        #loop over each observation
        for visit in range(data.nvisit):
            ind = data.vis_num == visit
            #calculate astrophysical model (transit, eclipse, phase curve)
            for m in astrophysical_models:
                self.physical_model *= models[m, params     #FIXME

            #calculate instrument systematics model
            for m in systematic_models:

        """for i in range(data.nvisit):
            ind = data.vis_num == i
            self.phase[ind] = (data.time[ind] - p.t0[i])/p.per[i] - np.floor((data.time[ind] - p.t0[i])/p.per[i])
            if data.lc_type == "eclipse": 
                self.lc[ind] = get_elc(data.time[ind], p, data, i) 	
                self.lc_nosys[ind] = get_elc(data.time[ind], p, data, i) 	
            elif data.lc_type == "transit": 
                self.lc[ind] = get_tlc(data.time[ind], p, data, i) 
                self.lc_nosys[ind] = get_tlc(data.time[ind], p, data, i) 
            elif data.lc_type == "phase_curve":
                self.eclipse_model[ind] = get_elc(data.time[ind], p, data, i) 	
                self.transit_model[ind] = get_tlc(data.time[ind], p, data, i) 
                self.phase_model[ind] = get_phaselc(data.time[ind], p, data, i)
                self.lc[ind] = self.transit_model[ind] + (self.eclipse_model[ind] - 1.)*self.phase_model[ind]
                self.lc_nosys[ind] = self.transit_model[ind] + (self.eclipse_model[ind] - 1.)*self.phase_model[ind]
            else: assert False, "Unknown option; supported light curve types are 'transit', 'eclipse', and 'phase_curve'"
            self.phase_hr = np.linspace(self.phase.min()-0.05, self.phase.max()+0.05, 1000)
            self.t_hr = self.phase_hr*p.per[0]+p.t0[0]
            self.lc_hr = get_tlc(self.t_hr, p, data, 0)

            if flags['divide-white'] == False:
                S = np.ones_like(self.transit_model[ind])+p.scale[i]*data.scan_direction[ind]
                D = np.ones_like(self.transit_model[ind])+p.r3[i]*data.t_delay[ind]
                self.all_sys[ind] = data.flux[ind]/self.lc[ind]
                self.lc[ind] *= (p.c[i]*S + p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2)*(1.0-np.exp((-p.r1[i]*data.t_orb[ind]-p.r2[i]-D)))*(p.A1[i]*np.sin(2.*np.pi*(data.t_vis[ind] - p.phi1[i])/p.P1[i]))
                self.data_corr[ind] = data.flux[ind]/((p.c[i]*S + p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2)*(1.0-np.exp(-p.r1[i]*data.t_orb[ind]-p.r2[i]-D))*(p.A1[i]*np.sin(2.*np.pi*(data.t_vis[ind] - p.phi1[i])/p.P1[i])))
                self.sys[ind] = p.c[i]+ p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2
                self.data_normalized[ind] = data.flux[ind]/(1.-np.exp((-p.r1[i]*data.t_orb[ind]-p.r2[i]-D))) - self.transit_model[ind]*p.c[i]*(S-1.)
            else:
                self.lc[ind] *= (p.c[i] + p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2)*data.all_sys[ind]
                self.data_corr[ind] = data.flux[ind]/((p.c[i] + p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2)*data.all_sys[ind])"""


        self.resid = data.flux - self.lc 
        self.chi2 = np.sum((self.resid/data.err)**2)		
        self.chi2red = self.chi2/data.dof
        self.rms = 1.0e6*np.sqrt(np.mean((self.resid/data.flux)**2))
        self.rms_predicted = 1.0e6*np.sqrt(np.mean((data.err/data.flux)**2))
        self.data = data
        self.ln_likelihood = -0.5*(np.sum((self.resid/data.err)**2 + np.log(2.0*np.pi*(data.err)**2)))
        self.bic = -2.*self.ln_likelihood + data.nfree_param*np.log(data.npoints)

                


