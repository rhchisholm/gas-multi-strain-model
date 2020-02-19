# Imports
import elfi
import GPy
import os

import logging
logging.basicConfig(level=logging.INFO)

import numpy as np
import matplotlib.pyplot as plt
import scipy as sc

from glmnet import LogitNet
from scipy.stats import gaussian_kde
from itertools import combinations
import datetime
import time

import matlab.engine

# Set an arbitrary seed and a global random state to keep the randomly generated quantities the same between runs
#seed = 20180530
#seed = 20190517
seed = 20190605
np.random.seed(seed)

eng = matlab.engine.start_matlab()

# Global dict for BOLFIRE
global_params = {}

# Global arrays for results
global_param_arr = []
global_lambda_arr = []
global_intercept_arr = []
global_coef_arr = []

# Global marginal
global_marginal = None

# Parameters
#R0 = p_1, DI = p_2, sigma = p_3

n_m = 320
n_theta = 80

p_1_true = 2 
p_2_true = 10
p_3_true = 1

batch_size = 1
random_state = None


def GASmodel2(p_1, p_2, p_3, batch_size=1, random_state=None):
    random_state = random_state or np.random
    
    logxt = np.log(10)
    omegat = 0.1
    alphat = 3;
    burnint = 30;
    
    global global_params
    global_params = {'p_1': p_1, 'p_2': p_2, 'p_3': p_3}
    
    # Make inputs 2d arrays for numpy broadcasting with w
    p_1 = np.asanyarray(p_1).reshape((-1, 1))
    p_3 = np.asanyarray(p_3).reshape((-1, 1))
    omega = np.asanyarray(omegat).reshape((-1, 1))
    x = np.asanyarray(np.exp(logxt)).reshape((-1, 1))
    alpha = np.asanyarray(alphat).reshape((-1, 1))
    p_2 = np.asanyarray(p_2).reshape((-1, 1))
    burnin = np.asanyarray(burnint).reshape((-1, 1))
    
    random_state = random_state or np.random
    
    # MATLAB array initialized with Python's list
    p_1 = matlab.double(p_1.tolist())
    p_3 = matlab.double(p_3.tolist())
    omega = matlab.double(omega.tolist())
    x = matlab.double(x.tolist())
    alpha = matlab.double(alpha.tolist())
    p_2 = matlab.double(p_2.tolist())
    burnin = matlab.double(burnin.tolist())
    
    simM = eng.gas_model(p_1,p_2,p_3,omega,x,alpha,burnin)
    
    # Convert back to numpy array
    sim = np.atleast_2d(simM).reshape(-1,12)
    sim = sim[:,0:11]
    
    return sim

def GASmodel2_local(p_1, p_2, p_3, batch_size=1, random_state=None):
    random_state = random_state or np.random
    
    logxt = np.log(10)
    omegat = 0.1
    alphat = 3;
    burnint = 30;
    
    p_1 = p_1*np.ones(batch_size)
    p_2 = p_2*np.ones(batch_size)
    p_3 = p_3*np.ones(batch_size)
    
    # Make inputs 2d arrays for numpy broadcasting with w
    p_1 = np.asanyarray(p_1).reshape((-1, 1))
    p_3 = np.asanyarray(p_3).reshape((-1, 1))
    omega = np.asanyarray(omegat).reshape((-1, 1))
    x = np.asanyarray(np.exp(logxt)).reshape((-1, 1))
    alpha = np.asanyarray(alphat).reshape((-1, 1))
    p_2 = np.asanyarray(p_2).reshape((-1, 1))
    burnin = np.asanyarray(burnint).reshape((-1, 1))
    
    random_state = random_state or np.random
    
    # MATLAB array initialized with Python's list
    p_1 = matlab.double(p_1.tolist())
    p_3 = matlab.double(p_3.tolist())
    omega = matlab.double(omega.tolist())
    x = matlab.double(x.tolist())
    alpha = matlab.double(alpha.tolist())
    p_2 = matlab.double(p_2.tolist())
    burnin = matlab.double(burnin.tolist())
    
    simM = eng.gas_model(p_1,p_2,p_3,omega,x,alpha,burnin)
    
    # Convert back to numpy array
    sim = np.atleast_2d(simM).reshape(-1,12)
    sim = sim[:,0:11]
       
    return sim
    
# summary statistics functions

def avgprev(x):  
    y = x[0,0]
    z = y.reshape(1,1)
    return z

def avgdiv(x):    
    y = x[0,1]
    z = y.reshape(1,1)
    return z

def maxnum1str(x):    
    y = x[0,2]
    z = y.reshape(1,1)
    return z

def avgtimestrobs(x):   
    y = x[0,3]
    z = y.reshape(1,1)
    return z

def maxtimestrobs(x):    
    y = x[0,4]
    z = y.reshape(1,1)
    return z

def numstrobs(x):    
    y = x[0,5]
    z = y.reshape(1,1)
    return z

def varprev(x):   
    y = x[0,6]
    z = y.reshape(1,1)
    return z
    
def vardiv(x):   
    y = x[0,7]
    z = y.reshape(1,1)
    return z

def avgtimerepeat(x):
    y = x[0,8]
    z = y.reshape(1,1)
    return z

def vartimerepeat(x):
    y = x[0,9]
    z = y.reshape(1,1)
    return z
    
def npmi(x):
	y = x[0,10]
	z = y.reshape(1,1)
	return z
	
def divalli(x):
	y = x[0,11]
	z = y.reshape(1,1)
	return z
    
# Samples for LFIRE
def sample_from_likelihood(params, batch_size=1, random_state=None):
    random_state = random_state or np.random
    p_1_prior = params['p_1']
    p_2_prior = params['p_2']
    p_3_prior = params['p_3']
    sample = GASmodel2_local(p_1_prior, p_2_prior, p_3_prior, batch_size, random_state)
    return sample

def sample_from_marginal(batch_size=1, random_state=None):
    random_state = random_state or np.random
    p_1_prior = random_state.uniform(0, 5, batch_size)
    p_2_prior = random_state.uniform(5, 15, batch_size)
    #p_2_prior = random_state.uniform(15, 25, batch_size)
    p_3_prior = random_state.uniform(0.5, 1, batch_size)
    sample = GASmodel2_local(p_1_prior, p_2_prior, p_3_prior, batch_size, random_state)
    return sample


# Marginal data
p_1_prior = np.random.uniform(0, 5, n_m)
p_2_prior = np.random.uniform(5, 15, n_m)
#p_2_prior = np.random.uniform(15, 25, n_m)
p_3_prior = np.random.uniform(0.5, 1, n_m)

# Marginal data reduced to summary statistics
# Calculate global marginal data
global_marginal = sample_from_marginal(n_m, None)
save_gm_name='GM_syn_DI5_15_'+str(n_m)
#save_gm_name='GM_syn_DI15_25_'+str(n_m)
np.savez(save_gm_name,global_marginal=global_marginal)

# OR load global marginal data
#save_gm_name='GM_syn_DI5_15_320.npz'
#save_gm_name='GM_syn_DI15_25_320.npz'
#saved_marginal = np.load(save_gm_name)
#global_marginal = saved_marginal['global_marginal']
    
# Synthetic data: 
# Use pre-generated data with DI = 10 or 20, R0=2, sigma = 1 (rng=1)
y_obs = np.array([[5.8847,1.8307,26.0000,3.8889,9.0000,9.0000,29.0527,4.4326,1.8000,1.7000,-0.1976]]) # DI = 10yr
# y_obs = np.array([[10.2160571656385, 2.54025816073426, 22, 6, 12,14, 33.4501948260742, 2.63235705270612, 1.83333333333333, 1.76666666666667, -0.213817505572094]]) # DI = 20 yr
# OR regenerate synthetic data:
# random_state_obs = np.random.RandomState(0)
# y_obs = GASmodel2_local(p_1_true, p_2_true, p_3_true, random_state=random_state_obs)

# Build an elfi model
m = elfi.ElfiModel()

# Priors
p_1 = elfi.Prior('uniform', 0, 5, model=m)
p_2 = elfi.Prior('uniform',  5, 10, model=m)
#p_2 = elfi.Prior('uniform',  15, 10, model=m)
p_3 = elfi.Prior('uniform',  0.5, 1, model=m)

priors = [p_1, p_2, p_3]

# Simulator
Y = elfi.Simulator(GASmodel2, *priors, observed=y_obs)

# Summary statistics
sumstats = []

sumstats.append(elfi.Summary(avgprev, Y, name='AvgPrev'))
sumstats.append(elfi.Summary(avgdiv, Y, name='AvgDiv'))
sumstats.append(elfi.Summary(maxnum1str, Y, name='MaxNum1Str'))
sumstats.append(elfi.Summary(avgtimestrobs, Y, name='AvgTimeStrObs'))
sumstats.append(elfi.Summary(maxtimestrobs, Y, name='MaxTimeStrObs'))
sumstats.append(elfi.Summary(numstrobs, Y, name='NumStrObs'))
sumstats.append(elfi.Summary(varprev, Y, name='VarPrev'))
sumstats.append(elfi.Summary(vardiv, Y, name='VarDiv'))
sumstats.append(elfi.Summary(avgtimerepeat, Y, name='AvgTimeBtwRepeatObs'))
sumstats.append(elfi.Summary(vartimerepeat, Y, name='VarTimeBtwRepeatObs'))
sumstats.append(elfi.Summary(npmi, Y, name='NMPI'))

# LFIRE method
def LFIRE(X, Y, n_m=n_m, n_theta=n_theta, random_state=random_state):
    """ LFIRE method as a distance node.
    
    Parameters
    ----------
    X: np.ndarray
        Simulated summary statistics
    Y: np.ndarray
        Observed summary statistics
    n_m: int
        Number of simulations from marginal
    n_theta: int
        Number of simulations from likelihood
    random_state: np.random
        Random state for random generators
        
    Output
    ------
    res: np.ndarray
        Negative value of the posterior function
    
    """
    
    random_state = random_state or np.random
    
    logreg = LogitNet(alpha=1, n_splits=10,
                      n_jobs=-1, verbose=0)
    
    # Global variables
    global global_params
    params = global_params
    
    global global_marginal
    marginal = global_marginal
    
    # Generate training data
    X_sim = sample_from_likelihood(params,
                                   batch_size=n_theta-1,
                                   random_state=random_state)
    
    X_ = np.concatenate((X, X_sim, marginal), axis=0)
    Y_ = np.concatenate((np.ones(n_theta), -1*np.ones(n_m)))
    
    # Fit the model
    logreg.fit(X_, Y_)
    
    # (Unnormalized) log posterior value
    res = logreg.intercept_ +  np.sum(np.multiply(logreg.coef_, Y))
    
    # Store results
    global global_param_arr
    global global_lambda_arr
    global global_intercept_arr
    global global_coef_arr
    
    global_param_arr.append([params['p_1'][0], params['p_2'][0], params['p_3'][0]])
    global_lambda_arr.append(logreg.lambda_best_[0])
    global_intercept_arr.append(logreg.intercept_)
    global_coef_arr.append(logreg.coef_.ravel())
    
    # Negative posterior value
    return np.atleast_2d(-0.5 * np.exp(res))
    
# Distance (LFIRE)
lfire = elfi.Distance(LFIRE, *sumstats, name='LFIRE')

# Define polynomial mean function
class Polynomial(GPy.core.Mapping):
    
    def __init__(self, input_dim, output_dim, name='polymap'):
        super(Polynomial, self).__init__(input_dim=input_dim, output_dim=output_dim, name=name)
        
        self.A = GPy.core.Param('A', np.random.rand(self.input_dim, self.output_dim))
        self.B = GPy.core.Param('B', np.random.rand(self.input_dim, self.output_dim))
        self.C = GPy.core.Param('C', np.random.rand())
        
        self.link_parameter(self.A)
        self.link_parameter(self.B)
        self.link_parameter(self.C)
        
        self.A.constrain_positive()
    
    def f(self, X):
        XX = np.power(X, 2)
        return np.dot(XX, self.A) + np.dot(X, self.B) + self.C
    
    def update_gradients(self, dL_dF, X):
        XX = np.power(X, 2)
        self.A.gradient = np.dot(XX.T, dL_dF)
        self.B.gradient = np.dot(X.T, dL_dF)
        self.C.gradient = np.sum(dL_dF, axis=0)
        
# Define a custom target model
parameter_names = ['p_1', 'p_2','p_3']

bounds = {'p_1': (0, 5),
          'p_2': (5, 15),
          #'p_2': (15, 25),
          'p_3': (0.5, 1)}

mf = Polynomial(3, 1)

# Radial Basis Function kernel, aka squared-exponential, exponentiated quadratic or 
# Gaussian kernel: k(r) = sigma^2 exp(-1/2 r^2) 
kernel = GPy.kern.RBF(input_dim=len(parameter_names), ARD=True)

# Matern 5/2 kernel: 𝑘(𝑟)=𝜎^2(1+√5𝑟+5/3𝑟^2)exp(−√5𝑟)
#kernel = GPy.kern.Matern52(input_dim=len(parameter_names), ARD=True)

# The White kernel is used as part of a sum-kernel where it explains the
# noise-component of the signal
kernel += GPy.kern.White(input_dim=len(parameter_names))

target_model = elfi.GPyRegression(parameter_names=parameter_names,
                                  kernel=kernel,
                                  mean_function=mf,
                                  bounds=bounds)

# Bolfire with a custom target
acq_noise_var = np.power([2e-3, 2e-3, 2e-3], 2)

bolfire = elfi.BOLFI(lfire, batch_size=1, initial_evidence=50,
                     update_interval=10, target_model=target_model,
                     async_acq=False,
                     acq_noise_var=acq_noise_var)

                     
# Set parameters
step_size = 10
max_evidence = 500

# Run inference
for k in np.arange(bolfire.n_initial_evidence, max_evidence + step_size, step_size):
    
    print('n_evidence = {}'.format(k))
    
    print('Beginning bolfire fit at {}'.format(datetime.datetime.now()))
    elfi.set_client('multiprocessing')
    post = bolfire.fit(n_evidence=k, bar=True)
    print('Finished bolfire fit at {}'.format(datetime.datetime.now()))
    
    X = bolfire.target_model.X
    Y = bolfire.target_model.Y
    
    save_model_name1='FIT_syn_DI5_15_'+str(k)
    save_model_name2='PARAMS_syn_DI5_15_'+str(k) 
    #save_model_name1='FIT_syn_DI15_25_'+str(k)
    #save_model_name2='PARAMS_syn_DI15_25_'+str(k)   
    np.savez(save_model_name2,global_param_arr=global_param_arr,global_lambda_arr=global_lambda_arr,global_intercept_arr=global_intercept_arr,global_coef_arr=global_coef_arr)   
    np.savez(save_model_name1,X=X,Y=Y,param_array=bolfire.target_model._gp.param_array,bounds=bolfire.target_model.bounds,y_obs=y_obs)

print('Beginning posterior samples at {}'.format(datetime.datetime.now()))

num_iterations = 100000  
elfi.set_client('multiprocessing')   
result=bolfire.sample(num_iterations)

save_model_name3='SAMPLES_syn_DI5_15_'+str(num_iterations)
#save_model_name3='SAMPLES_syn_DI15_25_'+str(num_iterations)
np.savez(save_model_name,p_1=result.samples['p_1'],p_2=result.samples['p_2'],p_3=result.samples['p_3'])
   
eng.quit()