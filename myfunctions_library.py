#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Marta Zafffaroni

This library include the user defined functions
"""
import numpy as np
from sympy import symbols, solve

######################################################################################################################################################
# This function integrate the dynamical system of the single host-multi vector model

def ODE_solve(state_0, time, params):
    # assign the value to each parameter : NB = DO NOT CHANGE THE ORDER IN THE 
    # PARAMS DICTIONARY
    s, i, x_r, z_r, x_t, z_t = state_0

    N_P = params[0]
    N_T = params[1]
    Lambda_R = params[2]
    delta_R = params[3]
    Lambda_T = params[4]
    delta_T = params[5]
    alpha = np.round(params[6],0)
    nu = params[7]
    rho = params[8]
    theta = params[9]
    r = params[10]
    k = params[11]
    epsilon_R = params[12]
    mu =  params[13]
    gamma = params[14]
    pigr = params[15]
    epsilon_T = params[16]
    tau = params[17]
    beta = params[18]
    
    k_ave = k
    
    f = 1/(1 +  (nu/(k_ave**(1-beta))*(x_r + z_r)/(k**beta))**alpha)
    g = (1 + (nu/(k_ave**(1-beta))*(x_r + z_r)/(k**beta))**alpha)

 
    state_dot = [ (rho + theta)*i - (Lambda_R*delta_R*z_r*N_P + Lambda_T*delta_T*f*z_t*N_P)*s/N_P ,
                 (Lambda_R*delta_R*z_r*N_P + Lambda_T*delta_T*f*z_t*N_P)*s/N_P - (rho + theta)*i,
                 r*(x_r+z_r)*(1-(x_r+z_r)/k) - mu*x_r - Lambda_R*epsilon_R*i/N_P*x_r + gamma*z_r,
                 Lambda_R*epsilon_R*i/N_P*x_r - (gamma + mu)*z_r,
                 (1-pigr)*tau*N_T - Lambda_T*epsilon_T*f*i/N_P*x_t + gamma*z_t - g*tau*x_t,
                 pigr*tau*N_T + Lambda_T*epsilon_T*f*i/N_P*x_t - (gamma + g*tau)*z_t ]

    return state_dot
####################################################################################################################################################
# This function compute the epidemic equilibrium 

def compute_equilibrium(params, R0R, R0T, f_comp, g_comp):

    i_r = (params['Lambda_R']*params['delta_R']*params['k']*(1-params['mu']/params['r']))/(params['rho'] + params['theta'])
    i_t = (params['Lambda_T']*params['delta_T']*params['N_T'])/(params['rho'] + params['theta'])
    a_r = (params['Lambda_R']*params['epsilon_R'])/(params['rho'] + params['theta'])
    a_t = (params['Lambda_T']*params['epsilon_T'])/(params['rho'] + params['theta'])
    m = (params['gamma'] + params['mu'])/(params['rho'] + params['theta'])
    i = (params['pigr']*params['tau'])/(params['rho'] + params['theta'])
    h = params['gamma']/(params['rho'] + params['theta'])
    e = params['tau']/(params['rho'] + params['theta'])
 
    if params['pigr'] == 0.0:
        if params['r'] < params['mu']:
            k_eq = 0.0
            z_r_eq = 0.0
            n_t_eq = 1.0
            f_in = 1
            g_in = 1
            
            R0R = i_r*a_r/m
            R0T = i_t*a_t*f_in**2/(g_in*(h+e*g_in))
            
            if R0T < 1:
                i_eq = 0.0
                z_t_eq = 0.0
            else:
                i_eq = (R0T - 1)/(R0T + a_t/(h+e))
                z_t_eq =  (a_t*i_eq)/(a_t*i_eq + h + e)
        else:
            k_eq = 1.0
            n_t_eq = 1/g_comp
            
            if (R0T + R0R) < 1:
                i_eq = 0.0
                z_r_eq = 0.0
                z_t_eq = 0.0
            else:
                a2 = - (R0R * a_t/(h+e*g_comp) * f_comp + R0T * a_r/m + a_r*a_t*f_comp/(m*(e*g_comp + h)))
                a1 = (a_t/(e*g_comp + h)*f_comp*(R0R-1) + a_r/m*(R0T - 1) - (R0R + R0T))
                a0 = (R0R + R0T - 1)
                i_eq = (-a1 - np.sqrt(a1**2-4*a2*a0))/(2*a2) 
                z_r_eq = a_r*i_eq/(a_r*i_eq + m)
                z_t_eq = a_t*f_comp*i_eq/(g_comp*(a_t*f_comp*i_eq + h + e*g_comp))
    elif params['pigr'] > 0.0:
        print('pigr >1')
        if params['r'] < params['mu']:
            k_eq = 0.0
            z_r_eq = 0.0
            n_t_eq = 1.0
            
            I = symbols('I')
            root = (i_t * i/(I *(a_t*I+h+e))) + (i_t*a_t/(a_t*I+h+e)) - (1/(1-I))
            i_root = solve(root, I)
            i_root_new = [expr.evalf(chop=True) for expr in i_root] # need this to get rid of *I
            for idx in i_root_new:
                if idx >= 0 and idx <= 1:
                    i_eq = idx
            z_t_eq = (i+a_t*i_eq)/(a_t*i_eq+h+e)
        else:
            k_eq = 1.0
            n_t_eq = 1/g_comp
            
            I = symbols('I')
            root = (i_r*a_r/(a_r*I + m)) + (i_t*f_comp*i*g_comp/(I*g_comp*(a_t*f_comp*I+e*g_comp+h))) + (i_t*a_t*f_comp**2/(g_comp*(a_t*f_comp*I+e*g_comp+h))) - (1/(1-I))
            i_root = solve(root, I)
            i_root_new = [expr.evalf(chop=True) for expr in i_root] # need this to get rid of *I
            for idx in i_root_new:
                if idx >= 0 and idx <= 1:
                    i_eq = idx
            z_r_eq = a_r*i_eq/(a_r*i_eq + m)
            z_t_eq = (i*g_comp + a_t*f_comp*i_eq)/(g_comp*(a_t*f_comp*i_eq + e*g_comp + h))
        
    s_eq = 1.0 - i_eq
    x_r_eq = k_eq - z_r_eq
    x_t_eq = n_t_eq - z_t_eq
    equilibria = [s_eq, i_eq, x_r_eq, z_r_eq, x_t_eq, z_t_eq]
    return equilibria