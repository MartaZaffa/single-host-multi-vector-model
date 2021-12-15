#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Marta Zaffaroni

This library cointains the user defined functions used to made plots  
"""
import matplotlib.pyplot as plt # module for graphical representation
import numpy as np
from matplotlib import colors
import seaborn as sns

# Define tick size
plt.rc('xtick', labelsize = 12) 
plt.rc('ytick', labelsize = 12) 


    
##################################################################################################################################################
# Plot function used in dynamical system.py
def competition_function(N_C, h, a_low, a_ave, a_high, nu_low, nu_ave, nu_high, h_ref):
    
    f, ax1 = plt.subplots()
    ax1.plot(N_C, 1/(1+(nu_low*N_C/h_ref)**a_low),'--g',label = r'$\nu =2.0, \alpha = 0.5$')
    ax1.plot(N_C, 1/(1+(nu_high*N_C/h_ref)**a_low),'--r',label = r'$\nu = 22.0, \alpha = 0.5$')
    ax1.plot(N_C, 1/(1+(nu_low*N_C/h_ref)**a_high),'g',label = r'$\nu = 2.0, \alpha = 5.0$')
    ax1.plot(N_C, 1/(1+(nu_high*N_C/h_ref)**a_high),'r',label = r'$\nu = 22.0, \alpha = 5.0$')
    ax1.plot(N_C, 1/(1+(nu_ave*N_C/h_ref)**a_ave),'k',label = r'$\nu = 12.0, \alpha = 1.0$')
    ax1.set_xlabel('Colonizer aphids $N_R$', fontsize = 20)
    ax1.set_ylabel('Interference function $f(N_R)$', fontsize = 20)
    ax1.set_aspect(1./ax1.get_data_ratio()) # for having squared boxes
    ax1.legend()

############################################################################################################################    
# Plot function used in Sens_An.py and Sens_An_alpha.py

def sa_1_param(i_eq, i_eq_2, R0, R0C, R0M, R0_2, R0C_2, R0M_2, par_name, par_values,x_lim):
    
    if par_name == 'k':
        val = par_values
        xlab = '$\kappa$' 
    elif par_name == 'mu':
        val = par_values
        xlab = '$\mu$'
    elif par_name == 'rho':
        val = par_values
        xlab = r'$\rho$'
  
    small = 1
    big = 2
    
    # fig1 is a (3,1) subplot showing the response of R0 to variation in the selected 
    # parameter under three intereference strength scenarios
    
    fig1, axs = plt.subplots(np.size(R0,1),1)
    
    plt.rcParams['axes.titley'] = 1.0    # y is in axes-relative co-ordinates.
    plt.rcParams['axes.titlepad'] = -14  # pad is in points...
    
    ################# INTERMEDIATE INTERFERENCE ############################## 
    
    axs[1].plot(val, R0C[:,1], 'b--', label = '$R_0^C$', linewidth = small)
    axs[1].plot(val, R0M[:,1], 'r--', label = '$R_0^M$', linewidth = small)
    axs[1].plot(val, R0C[:,1] + R0M[:,1], 'g--', label = '$R_0$', linewidth = big)
    axs[1].plot(val, R0C_2[:,1], color ='b', label = '$R_0^C$', linewidth = small)
    axs[1].plot(val, R0M_2[:,1], color ='r', label = '$R_0^M$', linewidth = small)
    axs[1].plot(val, R0C_2[:,1] + R0M_2[:,1], color ='g', label = '$R_0$', linewidth = big)
    axs[1].set_xlim(x_lim[1])
    axs[1].set_title("intermediate interference")

    if par_name == 'k':
        axs[1].ticklabel_format(style='sci', axis='x', scilimits=(3,3))
        axs[1].set_ylim((0-0.05, 2+0.05))
        
    ################# WEAK INTERFERENCE #################################### 
        
    axs[0].plot(val, R0C[:,0], 'b--', label = '$R_0^C$', linewidth = small)
    axs[0].plot(val, R0M[:,0], 'r--', label = '$R_0^M$', linewidth = small)
    axs[0].plot(val, R0C[:,0] + R0M[:,0], 'g--', label = '$R_0$', linewidth = big)
    axs[0].plot(val, R0C_2[:,0],'b', label = '$R_0^C$', linewidth = small)
    axs[0].plot(val, R0M_2[:,0],'r', label = '$R_0^M$', linewidth = small)
    axs[0].plot(val, R0C_2[:,0] + R0M_2[:,0],'g', label = '$R_0$', linewidth = big)
    axs[0].set_ylim((0-0.05, 11+0.05))
    axs[0].set_xlim(x_lim[0])
    axs[0].set_title("weak interference")

    if par_name == 'k':
        axs[0].ticklabel_format(style='sci', axis='x', scilimits=(3,3))
        
    ################# STRONG INTERFERENCE ####################################
    
    axs[2].plot(val, R0C[:,2], 'b--', label = '$R_0^C$', linewidth = small)
    axs[2].plot(val, R0M[:,2], 'r--', label = '$R_0^M$', linewidth = small)
    axs[2].plot(val, R0C[:,2] + R0M[:,2], 'g--', label = '$R_0$', linewidth = big)
    axs[2].plot(val, R0C_2[:,2],'b', label = '$R_0^C$', linewidth = small)
    axs[2].plot(val, R0M_2[:,2],'r', label = '$R_0^M$', linewidth = small)
    axs[2].plot(val, R0C_2[:,2] + R0M_2[:,2],'g', label = '$R_0$', linewidth = big)
    axs[2].set_ylim((0-0.05, 1.5+0.05))
    axs[2].set_xlim(x_lim[2])
    axs[2].set_title("strong interference")

    if par_name == 'k':
        axs[0].ticklabel_format(style='sci', axis='x', scilimits=(3,3))
        axs[1].ticklabel_format(style='sci', axis='x', scilimits=(3,3))
        axs[2].ticklabel_format(style='sci', axis='x', scilimits=(3,3))
    
    axs[1].set_ylabel('$R_0$', fontsize = 12)
    axs[2].set_xlabel(xlab, fontsize = 12)
    axs[0].axhline(y = 1, c = "gray", linewidth = 0.5)
    axs[1].axhline(y = 1, c = "gray", linewidth = 0.5)
    axs[2].axhline(y = 1, c = "gray", linewidth = 0.5)
    axs[1].set_ylim((0-0.05, 2+0.05))
    axs[0].set_aspect(1./axs[0].get_data_ratio())  
    axs[1].set_aspect(1./axs[1].get_data_ratio())   
    axs[2].set_aspect(1./axs[2].get_data_ratio())   
    
    # fig 2 is a plot showing the response of the infection incidence in plant and
    # of R0 to variation in the selected parameter (for direct and indirect interference)
 
    fig2, ax = plt.subplots()

    ax.plot(val, i_eq_2[:,1], linewidth = big, color = 'k')
    ax.plot(val, i_eq[:,1],'k--', linewidth = big)
    ax.set_ylim((0-0.01, 1))
    ax.set_xlabel(xlab, fontsize = 26)
    ax.set_ylabel('Infection incidence in plant', fontsize = 15)
    ax.set_xlim(x_lim[1])
    ax.set_aspect(1./ax.get_data_ratio()) # for having squared boxes
    
    if par_name == 'k':
        ax.ticklabel_format(style='sci', axis='x', scilimits=(3,3))
        
    ax2 = fig2.add_axes(ax.get_position())
    ax2.set_facecolor("None")
    
    ax2.plot(val,R0[:,1],'g--', linewidth = big)
    ax2.plot(val,R0_2[:,1],linewidth = big, color = 'green')
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel('R0', color = 'green', fontsize = 15)
    
    ax2.set_ylim((0-0.05, 2+0.05))
    ax2.axhline(y = 1, c = "gray", linewidth = 0.5)
    ax2.tick_params(bottom=0, top=0, left=0, right=1, 
                labelbottom=0, labeltop=0, labelleft=0, labelright=1)
    ax2.set_xlim(x_lim[1])
    ax2.set_aspect(1./ax2.get_data_ratio()) # for having squared boxes
    
    if par_name == 'k':
        ax2.ticklabel_format(style='sci', axis='x', scilimits=(3,3))
    
    
 ###############################################################################################################################   
# Plot function used in Sens_An.py

def sa_2_param(R0, R0C, R0M, par_name, par_values_1, par_values_2):
    
    c2 = (128/225, 128/225, 128/225, 0)
    c1 = (0/225, 0/225, 0/225, 1)
    
    heat_cmap = 'seismic'
    
    if par_name[0] == 'k':
        val = par_values_1
        xlab = '$\kappa$' 
        arr = -1
    elif par_name[0] == 'mu':
        val = par_values_1
        xlab = '$\mu$'
        arr = 2
            
    if par_name[1] == 'k':
        val_2 = par_values_2
        xlab_2 = '$\kappa$' 
        arr_2 = 0
    elif par_name[1]  == 'mu':
        val_2 = par_values_2
        xlab_2 = '$\mu$'
        arr_2 = 2
    elif par_name[1]  == 'rho':
        val_2 = par_values_2
        arr_2 = 2
        xlab_2 = r'$\rho$'
    

    div_1 = 125
    div_2 = 100
    
    cmap = colors.ListedColormap([c1,c2])
    bounds=[0, 1, R0.max()]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    
    # fig is a single plot showing the response of R0 to variation of two 
    # selected parameter at time
    
    fig, axs2 = plt.subplots()
    x_tick = np.round(val, arr)
    y_tick = np.round(val_2, arr_2)

    axs2 = sns.heatmap(R0.transpose(), yticklabels = y_tick, xticklabels=x_tick, vmin = 0, vmax=4,  cbar = True,  cmap = heat_cmap) 
    axs2.tick_params(axis=u'both', which=u'both',length=0)
    plt.yticks(fontsize = 14, rotation = 0)   
    plt.xticks(fontsize = 14, rotation = 0)
    
    for ind, label in enumerate(axs2.get_xticklabels()):
        if ind % div_1 == 0:
            label.set_visible(True)   
        else:
            label.set_visible(False)
             
    for ind, label in enumerate(axs2.get_yticklabels()):
        if ind % div_2 == 0: 
            label.set_visible(True)
            
        else:
            label.set_visible(False)
           
    axs2.set_ylabel(xlab_2, fontsize = 20)
    axs2.set_xlabel(xlab, fontsize = 20)
    axs2.set_aspect(1./axs2.get_data_ratio()) # for having squared
    axs2 = sns.heatmap(R0.transpose(), yticklabels = y_tick, xticklabels = x_tick, cmap=cmap, norm=norm,cbar = False)
    plt.yticks(fontsize = 20, rotation = 0)   
    plt.xticks(fontsize = 20, rotation = 0)
    axs2.set_ylabel(xlab_2, fontsize = 20)
    axs2.set_xlabel(xlab, fontsize = 20)
      
##################################################################################################################################
# This plot is used in dynamical_system.py

def equilibrium_function(I_eqm, b_1, b_2, v):
    
    # fig is a (2,1) subplot showing the graphical representation of the equilibrium
    # value of "I" (infected plant) for fraction of viruliferous transient 
    # aphids entering the system >0 (pigr > 0) and for NR=0 (first plot) and 
    # NR>0 (second plot)
    
    fig,[ax1, ax2] = plt.subplots(1,2)
    ax1.plot(I_eqm, b_1, label = '$b_2(\hat I)$')
    ax1.plot(I_eqm, v, label = '$v(\hat I)$')
    ax1.legend()
    ax1.set_ylim([0,500])
    ax1.set_xlim([0.6,1.1])
    ax1.set_aspect(1./ax1.get_data_ratio()) # for having squared boxes
    ax1.set_xlabel('$\hat I$', fontsize = 30)
    ax1.set_title('$\hat N_R = 0$, $u > 0$', fontsize = 30)
    
    ax2.plot(I_eqm, b_2, label = '$b_2(\hat I)$')
    ax2.plot(I_eqm, v, label = '$v(\hat I)$')
    ax2.legend()
    ax2.set_ylim([0,50])
    ax2.set_xlim([0.6,1.1])
    ax2.set_aspect(1./ax2.get_data_ratio()) # for having squared boxes
    ax2.set_xlabel('$\hat I$', fontsize = 30)
    ax2.set_title('$\hat N_R = 1$, $u > 0$', fontsize = 30)
    
