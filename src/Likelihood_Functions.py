import numpy as np
import pandas as pd
import radvel
from scipy import spatial
import scipy
from scipy.linalg import cho_factor, cho_solve


def Calculate_Likelihood_RadVel(vals, **kwargs):

    '''
    Function that uses internal Radvel Stuff to calculate the likelihood of a model
    given the values for each parameter and the data

    vals: collection of arrays with the value of each parameter
    NPlanets: number of planets in the system
    filename: path to the RV data
	GP: Whether or not we're using a GP Likelihood

    '''
    try:
        priors = kwargs['priors']
    except KeyError:
        print('No priors detected')
    try:
        NPlanets = kwargs['NPlanets']
    except KeyError:
        NPlanets = 1
    try:
        filename = kwargs['filename']
    except KeyError:
        print('Pass a filename to the likelihood')
    try:
        GP = kwargs['GP']
    except KeyError:
        GP = False
    try:
        hparam_dict = kwargs['hparam_dict']
    except KeyError:
        hparam_dict = None
    try:
        kernel = kwargs['kernel']
    except KeyError:
        kernel = None
    try:
        val_dict = kwargs['val_dict']
    except KeyError:
        val_dict = None
    try:
        optimizing = kwargs['optimizing']
    except KeyError:
        optimizing = False

    #update the value dictionary to have the latest values of vals
    counter = 0
    for key in priors.keys(): #ensures they are the varying Parameters
        val_dict[key] = vals[counter]
        counter += 1



    data = pd.read_csv(filename,sep=' ')

    t = np.array(data.time)
    vel = np.array(data.mnvel)
    errvel = np.array(data.errvel)
    try:
        tel_arr = np.array(data.tel)
        telgrps = data.groupby('tel').groups
        instnames = telgrps.keys()
    except AttributeError:
        tel_arr = np.repeat('test',len(t))
        instnames = ['test']

    telnames = []
    for inst in instnames:
        telnames.append(inst)
    telnames = np.array(telnames)

    vary_dict = {}
    for key in val_dict.keys():
        if key in priors.keys():
            vary_dict[key] = True
        else:
            vary_dict[key] = False


    if hparam_dict is not None:
        hnames = np.array(list(hparam_dict.keys()))





    params = radvel.Parameters(NPlanets,basis='per tc e w k')


    for i in range(NPlanets):
        params['per'+str(i+1)] = radvel.Parameter(value=val_dict['per'+str(i+1)], vary=vary_dict['per'+str(i+1)]) #period in days
        params['tc'+str(i+1)] = radvel.Parameter(value=val_dict['tc'+str(i+1)], vary=vary_dict['tc'+str(i+1)]) #Time of conjunction BJD
        params['e'+str(i+1)] = radvel.Parameter(value=val_dict['e'+str(i+1)],vary=vary_dict['e'+str(i+1)])
        params['w'+str(i+1)] = radvel.Parameter(value=val_dict['w'+str(i+1)],vary=vary_dict['w'+str(i+1)])
        params['k'+str(i+1)] = radvel.Parameter(value=val_dict['k'+str(i+1)], vary=vary_dict['k'+str(i+1)]) #initial guess for K amplitude, m/s

    params['dvdt'] = radvel.Parameter(value=val_dict['dvdt'],vary=vary_dict['dvdt']) #linear trend term
    params['curv'] = radvel.Parameter(value=val_dict['curv'],vary=vary_dict['curv']) #curvature term

    if GP:
        for key in hparam_dict.keys():
            params[key] = radvel.Parameter(value=val_dict[key], vary=vary_dict[key])
        # if kernel == 'QP':
        #     params['gp_amp'] = radvel.Parameter(value=val_dict['gp_amp'],vary=vary_dict['gp_amp']) #GP amplitude
        #     params['gp_per'] = radvel.Parameter(value=val_dict['gp_per'],vary=vary_dict['gp_per']) #GP period
        #     params['gp_explength'] = radvel.Parameter(value=val_dict['gp_explength'],vary=vary_dict['gp_explength']) #GP explength
        #     params['gp_perlength'] = radvel.Parameter(value=val_dict['gp_perlength'],vary=vary_dict['gp_perlength']) #GP perlength
        # elif kernel == 'KJ2':
        #     params['gp_amp0'] = radvel.Parameter(value=val_dict['gp_amp0'],vary=vary_dict['gp_amp0']) #GP amplitude
        #     params['gp_amplambda'] = radvel.Parameter(value=val_dict['gp_amplambda'],vary=vary_dict['gp_amplambda']) #Wavelength Scaling
        #     params['gp_per'] = radvel.Parameter(value=val_dict['gp_per'],vary=vary_dict['gp_per']) #GP period
        #     params['gp_explength'] = radvel.Parameter(value=val_dict['gp_explength'],vary=vary_dict['gp_explength']) #GP explength
        #     params['gp_perlength'] = radvel.Parameter(value=val_dict['gp_perlength'],vary=vary_dict['gp_perlength']) #GP perlength
        # elif kernel == 'KJ1':
        #     params['gp_amp_hires_pre'] = radvel.Parameter(value=val_dict['gp_amp_hires_pre'],vary=vary_dict['gp_amp_hires_pre']) #GP amplitude
        #     params['gp_amp_hires_post'] = radvel.Parameter(value=val_dict['gp_amp_hires_post'],vary=vary_dict['gp_amp_hires_post']) #GP amplitude
        #     params['gp_amp_carmenes'] = radvel.Parameter(value=val_dict['gp_amp_carmenes'],vary=vary_dict['gp_amp_carmenes']) #GP amplitude
        #     params['gp_amp_hpf_pre'] = radvel.Parameter(value=val_dict['gp_amp_hpf_pre'],vary=vary_dict['gp_amp_hpf_pre']) #GP amplitude
        #     params['gp_amp_hpf_post'] = radvel.Parameter(value=val_dict['gp_amp_hpf_post'],vary=vary_dict['gp_amp_hpf_post']) #GP amplitude
        #     params['gp_amp_NEID'] = radvel.Parameter(value=val_dict['gp_amp_NEID'],vary=vary_dict['gp_amp_NEID']) #GP amplitude
        #     params['gp_amp_NEID_post'] = radvel.Parameter(value=val_dict['gp_amp_NEID_post'],vary=vary_dict['gp_amp_NEID_post']) #GP amplitude
        #     params['gp_amp_spirou'] = radvel.Parameter(value=val_dict['gp_amp_spirou'],vary=vary_dict['gp_amp_spirou']) #GP amplitude
        #     params['gp_per'] = radvel.Parameter(value=val_dict['gp_per'],vary=vary_dict['gp_per']) #GP period
        #     params['gp_explength'] = radvel.Parameter(value=val_dict['gp_explength'],vary=vary_dict['gp_explength']) #GP explength
        #     params['gp_perlength'] = radvel.Parameter(value=val_dict['gp_perlength'],vary=vary_dict['gp_perlength']) #GP perlength

    for tel_suffix in instnames:
        params['gamma_'+tel_suffix] = radvel.Parameter(value=val_dict['gamma_'+tel_suffix], vary=vary_dict['gamma_'+tel_suffix])
        params['jit_'+tel_suffix] = radvel.Parameter(value=val_dict['jit_'+tel_suffix], vary=vary_dict['jit_'+tel_suffix])

    model = radvel.model.RVModel(params) #combine all the parameters into a Keplerian model under the hood

    likes = [] #a list of likelihoods for each instrument



    def initialize(tel_suffix):

        # Instantiate a separate likelihood object for each instrument.
        # Each likelihood must use the same radvel.RVModel object.
        try:
             indices = telgrps[tel_suffix]
        except:
             print('WARNING: telgrps did not initialize properly')
             indices = np.ones(len(t), dtype=bool)
        if GP:
            like = radvel.likelihood.GPLikelihood(model, t[indices], vel[indices],
                                                  errvel[indices], hnames[tel], suffix='_'+tel_suffix,
                                                  kernel='QuasiPer'
                                                  )
        else:
            like = radvel.likelihood.RVLikelihood(model, t[indices], vel[indices],
                                                  errvel[indices], suffix='_'+tel_suffix,
                                                  )
        like.params['gamma_'+tel_suffix] = radvel.Parameter(value=val_dict['gamma_'+tel_suffix], vary=True)
        like.params['jit_'+tel_suffix] = radvel.Parameter(value=val_dict['jit_'+tel_suffix], vary=True)
        # Add in instrument parameters
        likes.append(like)

    def initialize_KJ2():
        like = radvel.likelihood.ChromaticLikelihood(model=model, t=t, vel=vel,
        errvel=errvel, suffix=telnames, hnames=hnames,
        kernel_name="Chromatic_2", tel=tel_arr, telnames=telnames
        )
        return like

    def initialize_KJ1():
        like = radvel.likelihood.ChromaticLikelihood(model=model, t=t, vel=vel,
        errvel=errvel, suffix=telnames, hnames=hnames,
        kernel_name="Chromatic_1", tel=tel_arr, telnames=telnames
        )
        return like


    #loop through each instrument and append the likelihood
    if GP == False:
        for tel in instnames:
            initialize(tel)
    elif kernel == 'KJ2':
        like = initialize_KJ2()
        likes.append(like)
    elif kernel == 'KJ1':
        like = initialize_KJ1()
        likes.append(like)

    #merge all the likelihoods into a composite likelihood object for calculations
    like = radvel.likelihood.CompositeLikelihood(likes)


    return like.logprob()



def LikelihoodXPrior(x0, **kwargs):
    '''
    This function is the RV Likelihood times the prior function
    '''
    loglikelihood = Calculate_Likelihood_RadVel(x0, **kwargs)


    #know we cycle through each prior and add its non-constant portion to the exponent

    try:
        priors = kwargs['priors']
    except KeyError:
        print('You need a prior')
    try:
        optimizing = kwargs['optimizing']
    except KeyError:
        optimizing = False
    counter = 0
    for key in priors.keys():
        #print('Parameter: {}'.format(key))
        con, noncon = Prior_Components(priors[key], x0[counter])
        counter += 1
        loglikelihood += np.log(con*noncon)
    if optimizing:
        return loglikelihood*-1
    else:

        return loglikelihood


import numpy as np
import pandas as pd
import radvel
from scipy import spatial


def Prior_Components(priors, vals):
    '''
    This is a function where each kind of prior is divided into its "constant" portion,
    and its non-constant portion. The Non-constant portions will eventually be multiplied by
    the likelihood before we use the Laplace Approximation, and the constant portions will move
    outside of the integral, and their logarithms are added to the Evidence at the end

    -----------------------------
    Uniform

    U(a,b) = 1/(b-a) if a < x < b
             0 otherwise

    -----------------------------
    Gaussian

    G(a,b) = 1/(b*sqrt(2 pi)) * exp(-0.5*(x - a)^2/(b^2))


    -----------------------------
    Truncated Jeffreys

    TJ(a, b) = 1/x * (1/log(b/a))


    -----------------------------
    Truncated Modified Jeffreys

    TMJ(a, b) = 1/(a + x) * 1/(log(1 + b/a))


    -----------------------------
    Truncated Rayleigh Distribution

    TRD(a, b) = 1/(a^2) * exp(-x^2/2a^2) / (1 - exp(-b^2/(2a^2)))


    -----------------------------
    '''
    prior_func = priors[0] #by design

    #print('Prior function: {}'.format(prior_func))
    #print('Parameter value: {}'.format(vals))

    if prior_func == 'Uniform':
        constant = 1/(priors[2] - priors[1])
        #print('Lower: {}'.format(priors[1]))
        #print('Upper: {}'.format(priors[2]))
        if (vals > priors[1] and vals < priors[2]):
            nonconstant = 1.
        else:
            nonconstant = 0.
    elif prior_func == 'Gaussian':
        constant = 1/(priors[2] * np.sqrt(np.pi*2))
        nonconstant =  np.exp(-0.5 * (vals - priors[1])**2/(priors[2]**2))
    elif prior_func == 'Jeffreys':
        constant = 1/np.log(priors[2]/priors[1])
        nonconstant = 1/vals
    else:
        print(priors[0])
        print(prior_func)
        print("Is your prior Listed above? If not, you should add its functional form")
        constant = 1.
        nonconstant = 1.

    #print('contribution: {}'.format(np.log(constant*nonconstant)))


    return constant, nonconstant
