#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Marta Zaffaroni

This module create a dictionary with all the model parameters
"""
params = ({"N_P" :720.0,  # 0 Total number of plant
           "N_T" : 500.0, # 1 Average number of transient aphid per plaant in absent of resident aphid
           "Lambda_R" : 0.0485, # 2 Number of plant visited by a resident aphid
           "delta_R" : 0.035, # 3 Probability of virus transmission from the resident aphid to the plant
           "Lambda_T" : 8.5, # 4 Number of plant visited by a transient aphid
           "delta_T" : 0.035, # 5 Probability of virus transmission from the transient aphid to the plant
           "alpha": 1.0, # 6 Visiting (=Emigration) interference curvature
           "nu": 12.0, # 7 Visiting (=Emigration) interference strength (for direct interference)
           "rho": 0.02, # 8 Infected plant roguing rate
           "theta": 0.003, # 9 Plant harvesting rate
           "r" : 0.21, # 10 Intrinsic growth rate of resident aphids
           "k" : 50.0e3, # 11 Plant hosting capacity
           "epsilon_R": 0.024, # 12 Probability of virus transmission from the plant to the resident aphid
           "mu" : 0.08, # 13 Mortality rate of resident aphids
           "gamma" : 4.0, # 14 Virus clearance rate in aphid vectors
           "pigr": 0.0, # 15 Fraction of viruliferous transient aphids entering the system
           "epsilon_T": 0.024, # 16 Probability of virus transmission from the plant to the transient aphid
           "tau": 0.5, # 17 Transient aphids emigration rate in absence of resident aphids
           "beta": 0.0, # 18 Parameter of the Visiting (= interference function)
           "nu_2":3.0}) # 19                                 
