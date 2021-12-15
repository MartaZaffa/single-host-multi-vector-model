#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Marta Zaffaroni
"""
"""
This code 
simulates the temporal dynamics of the single host-multi vector system,
computes the epidemic equilibrium for the given parameters,
explores the changes in the interference function for changhing \nu 
show the values of "I" at the equilibrium for value of \pi > 0
"""
import numpy as np
from scipy.integrate import odeint
from myfunctions_library import ODE_solve, compute_equilibrium # import user-defined functions
from parameters import params # import the parameters defined by the users
import myplot_library # module inlcuding the user-defined plots

#############################################################################################################################
######### TEMPORAL DYNAMICS OF THE SYSTEM ###################################################################################
   
# Definition of initial states 
s_0 = params['N_P']-100 
i_0 = 100.0
if params['mu']>params['r']:
    x_r_0 = 0.0
    z_r_0 = 0.0
else:
    x_r_0 = 5.0
    z_r_0 = 5.0
x_t_0 = 5.0
z_t_0 = 1.0
    
state_0 = [s_0, i_0, x_r_0,  z_r_0, x_t_0, z_t_0]

# Definition of time horizon
t_in = 1
veg_seas = 500 # duration of the vegetative season (day)
step_t = 1

time_horiz = np.arange(t_in, t_in + veg_seas, step_t)

# Initialization of variables
state = np.zeros(( len(time_horiz),len(state_0))) 
# convert the dictionary "params" into an array to give as imput to odeint
params_array = list(params.values())

# Integration via the scipy's solver odeint     
state = odeint(ODE_solve, state_0, time_horiz, args = (params_array,))

k_1 = params['k']*(1-params['mu']/params['r'])

# Calculation of the basic reproduction number and of its components
# NB: parameter pigr should be equal to 0 to compute R0

k_ave = params['k']
if params['mu'] > params['r']:
    f_comp = 1
    g_comp = 1
else :
    f_comp = 1/(1 + (params['nu']/(k_ave**(1-params['beta'])) * k_1/(params['k']**(params['beta'])))**params['alpha']) # competition function for R_0 computation 
    g_comp = (1 + (params['nu']/(k_ave**(1-params['beta']))* k_1/(params['k']**(params['beta'])))**params['alpha'])
R0R = 1/(params['rho'] + \
     params ['theta'])*(params['Lambda_R']**2 * params['delta_R']*params['epsilon_R']*params['k']*max(0,(1-params['mu']/params['r'])))/\
    (params['mu'] + params['gamma'])
R0T = (1/(params['rho'] + \
     params ['theta']))*((params['Lambda_T']**2 * params['delta_T']*params['epsilon_T']*f_comp**2*params['N_T'])/\
    (g_comp*(params['tau']*g_comp + params['gamma'])))
R0 = R0R + R0T

#############################################################################################################################
######### COMPUTATION OF THE EPIDEMIC EQUILIBRIUM ###########################################################################

equilibrium = compute_equilibrium(params, R0R, R0T, f_comp, g_comp)

#############################################################################################################################
######## INTERFERENCE FUNCTION'S RESPONSE TO VARIATION OF PARAMETER nu ######################################################

N_R = np.arange(1.0,50.0e3,1.0)
h = np.arange(1.0,50.0e3,1.0)

a = [0.5, 1.0, 5.0]
nu = [6.0, 12.0, 18.0]

# sensitivity of f(Nc) to the value of alpha and nu, h is equal to \hat h
h_ref = 50e3

myplot_library.competition_function(N_R, h, a[0], a[1],a[2], nu[0], nu[1], nu[2], h_ref)

######### GRAPHICAL REPRESENTATION OF "I" VALUE AT ##########################################################################
######### THE EQUILIBRIUM FOR PARAMETER pigr > 0  ###########################################################################

params["pigr"] = 0.1
I_eqm = np.arange(0.0,1.1,0.001)
i_r = (params['Lambda_R']*params['delta_R']*params['k']*(1-params['mu']/params['r']))/(params['rho'] + params['theta'])
i_t = (params['Lambda_T']*params['delta_T']*params['N_T'])/(params['rho'] + params['theta'])
a_r = (params['Lambda_R']*params['epsilon_R'])/(params['rho'] + params['theta'])
a_t = (params['Lambda_T']*params['epsilon_T'])/(params['rho'] + params['theta'])
m = (params['gamma'] + params['mu'])/(params['rho'] + params['theta'])
u = (params['pigr']*params['tau'])/(params['rho'] + params['theta'])
h = params['gamma']/(params['rho'] + params['theta'])
e = params['tau']/(params['rho'] + params['theta'])
nu_a = (1-params['mu']/params['r'])/(params['k']**(params['beta']-1))*params['nu'] 

# N_R = 1 
b_2 = ((i_r*a_r/(a_r*I_eqm + m)) + (i_t*f_comp*u*g_comp/(I_eqm*g_comp*(a_t*f_comp*I_eqm+h+e*g_comp))) + 
       (i_t*a_t*f_comp**2/(g_comp*(a_t*f_comp*I_eqm + h + e*g_comp))))
v =  1/(1-I_eqm)

# N_R = 0
f_comp = 1
g_comp = 1

b_1 = ((i_t*f_comp*u*g_comp/(I_eqm*g_comp*(a_t*f_comp*I_eqm+h+e*g_comp))) + 
       (i_t*a_t*f_comp**2/(g_comp*(a_t*f_comp*I_eqm + h + e*g_comp))))

myplot_library.equilibrium_function(I_eqm, b_1, b_2, v)



   

