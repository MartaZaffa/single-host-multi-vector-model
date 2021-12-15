#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Marta Zaffaroni

This code allows to explore the responses of the basic reproduction number R0 to variations of selected
agricultural practices (modeled through plant hosting capacity (k), resident aphid 
mortality (mu), infected plant roguing rate (rho) under the direct and 
indirect interference scenarios. 

In the first part, the code allows to explore the response of R0 when one parameter
at a time is variated. In this par of code it is also possible to explore the influence of
interference strength parameter /nu on the result.
Furthermore, in this first part of code allows to analyse the response of the infection
incidence on plant when one of the aformentioned parameters is variated.
In the second part, the code allows to explore the response of R0 when two parameters
at a time are variated.
"""
from parameters import params
import numpy as np
from myfunctions_library import compute_equilibrium
import myplot_library

# Declare the parameter(s) that is (are) variated and define its (their) variation
# range

test_param = ['k'] # declare the name of the parameter(s) to vary
var_range = ([1e3,700e3]) # define the variation range limits
step = 500 

k_ave = params['k']

nu_values = [2.0,12.0,22.0] 
x_lim = ([150e3,700e3],[30e3,150e3],[10e3,80e3]) # limits for the x axis in the plot 
                                                 # (each plot corresponds to an interference strenght scenario)
sz_nu = len(nu_values)
   
if len(test_param) < 2: # change 1 parameter at time
    print('1 PARAMETER sensitivity analysis') 
    test_values = np.linspace(var_range[0], var_range[1], step)
   
    # variables initialization
    R0 = np.zeros((len(test_values), sz_nu))
    R0R = np.zeros((len(test_values), sz_nu))
    R0T = np.zeros((len(test_values), sz_nu))
    i_eq = np.zeros((len(test_values), sz_nu))
    Nc_eq = np.zeros((len(test_values), sz_nu))
    
    R0_2 = np.zeros((len(test_values), sz_nu))
    R0R_2 = np.zeros((len(test_values), sz_nu))
    R0T_2 = np.zeros((len(test_values), sz_nu))
    i_eq_2 = np.zeros((len(test_values), sz_nu))
    Nc_eq_2 = np.zeros((len(test_values), sz_nu))
    
    for idx1, v1 in enumerate(test_values):
        params[test_param[0]] = v1
        for col, nu in enumerate(nu_values):
           
            ################### interference function depends on the relative number of colonizer aphids ####################################################
            beta = 1  # DIRECT INTERFERENCE
            params['nu'] = nu*k_ave**(beta-1)
            if params['mu'] > params['r']:
                f_comp = 1
                g_comp = 1/f_comp
            else :
                f_comp = 1/(1 + (params['nu']*params['k']**(-beta)*params['k']*(1-params['mu']/params['r']))**params['alpha']) 
                g_comp = 1/f_comp
            R0R[idx1,col] = 1/(params['rho'] + \
                            params ['theta'])*(params['Lambda_R']**2 * params['delta_R']*params['epsilon_R']*params['k']*max(0,(1-params['mu']/params['r'])))/\
                            (params['mu'] + params['gamma'])
            R0T[idx1,col] = (1/(params['rho'] + \
                  params ['theta']))*((params['Lambda_T']**2 * params['delta_T']*params['epsilon_T']*f_comp**2*params['N_T'])/\
                 (g_comp*(params['tau']*g_comp + params['gamma'])))
            R0[idx1,col] = R0R[idx1,col] + R0T[idx1,col]
            equilibrium = compute_equilibrium(params, R0R[idx1,col], R0T[idx1,col], f_comp, g_comp)
            i_eq[idx1,col] = equilibrium[1] # infected plant at the equilibrium
            Nc_eq[idx1,col] = equilibrium[2] + equilibrium[3] # number of colonizer at the equilibrium
            
            ################### interference function depends on the absolute number of colonizer aphids ####################################################
            beta = 0 # INDIRECT  INTERFERENCE
            params['nu'] = nu*k_ave**(beta-1)
            if params['mu'] > params['r']:
                f_comp = 1
                g_comp = 1/f_comp
            else :
                f_comp = 1/(1 + (params['nu']*params['k']**(-beta)*params['k']*(1-params['mu']/params['r']))**params['alpha']) 
                g_comp = 1/f_comp
            R0R_2[idx1,col] = 1/(params['rho'] + \
                            params ['theta'])*(params['Lambda_R']**2 * params['delta_R']*params['epsilon_R']*params['k']*max(0,(1-params['mu']/params['r'])))/\
                            (params['mu'] + params['gamma'])
            R0T_2[idx1,col] = (1/(params['rho'] + \
                  params ['theta']))*((params['Lambda_T']**2 * params['delta_T']*params['epsilon_T']*f_comp**2*params['N_T'])/\
                 (g_comp*(params['tau']*g_comp + params['gamma'])))
            R0_2[idx1,col] = R0R_2[idx1,col] + R0T_2[idx1,col]
            equilibrium = compute_equilibrium(params, R0R_2[idx1,col], R0T_2[idx1,col], f_comp, g_comp)
            i_eq_2[idx1,col] = equilibrium[1] # infected plant at the equilibrium
            Nc_eq_2[idx1,col] = equilibrium[2] + equilibrium[3] # number of colonizer at the equilibrium
                           
    myplot_library.sa_1_param(i_eq, i_eq_2, R0, R0R, R0T, R0_2, R0R_2, R0T_2, test_param[0], test_values,x_lim)

else:
    print('2 PARAMETER sensitivity analysis')  # change 2 parameter at time
    
    beta =0 # beta = 0 (indirect interference), 
            # beta = 1 (direct interference)
    nu_value = 12 # set the interference strength scenario
    test_values = np.zeros((len(var_range), step))
    for idx, v in enumerate(var_range):
        test_values[idx] = np.linspace(v[0], v[1], step)
    
    # inizialisation of the variables
    R0 = np.zeros((step, step))
    R0R = np.zeros((step, step))
    R0T = np.zeros((step, step))
    
    if test_param[0] == 'rho':
        test_values_1 = test_values[0]
    else:
        test_values_1 = test_values[0]
    if test_param[1] == 'rho':
        test_values_2 = test_values[1]
    else:
        test_values_2 = test_values[1]
    
    for idx1, v1 in enumerate(test_values_1):
        params[test_param[0]] = v1
        for idx2, v2 in enumerate(test_values_2):
            params[test_param[1]] = v2
            if beta == 0:
                params['nu'] = nu_value/k_ave
            else :
                params['nu'] = nu_value
            if params['mu'] > params['r']:
                f_comp = 1
                g_comp = 1/f_comp
            else :
                f_comp = 1/(1 + (params['nu']*params['k']**(-beta)*params['k']*(1-params['mu']/params['r']))**params['alpha']) 
                g_comp = 1/f_comp 
                
            R0R[idx1,idx2] = 1/(params['rho'] + \
                                params ['theta'])*(params['Lambda_R']**2 * params['delta_R']*params['epsilon_R']*params['k']*max(0,(1-params['mu']/params['r'])))/\
                                (params['mu'] + params['gamma'])
            R0T[idx1,idx2] = (1/(params['rho'] + \
              params ['theta']))*((params['Lambda_T']**2 * params['delta_T']*params['epsilon_T']*f_comp**2*params['N_T'])/\
             (g_comp*(params['tau']*g_comp + params['gamma'])))
            R0[idx1,idx2] = R0R[idx1,idx2] + R0T[idx1,idx2]
         
    myplot_library.sa_2_param(R0, R0R, R0T, test_param, test_values_1, test_values_2) 

   