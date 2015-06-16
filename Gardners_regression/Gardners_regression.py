# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 13:35:04 2015

@author: atorlucci
"""

# Gardners regression
# algorithm for finding the parameters alpha and beta that fit Gardners general
# expression for the relationship between p-wave velocity and density in the
# least squares sense.
#
# Gardner's Equation:
# rho = alpha*(Vp**beta)
# this is obviously hyperbolic, but we can make it linear by taking the 
# natural log
# ln(rho) = ln(alpha*(Vp**beta))
# ln(rho) = ln(alpha) + ln(Vp**beta)
# ln(rho) = ln(alpha) + beta*ln(Vp)
# D       = A         + B   *V
# where A is the intercept and B is the slope (or gradient)

# import 3rd party libraries
from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# import local library
from las import LASReader


# import las file (see creating a synthetic ipython notebook)
L30 = LASReader('L-30.las', null_subs=np.nan)

# load log data into pandas dataframe
L30_df = pd.DataFrame(L30.data2d, columns=L30.curves.names, index=L30.data['DEPTH'])

# Create calculated Logs: P-wave Velocity from Sonic and P-wave Normal Incidence Impedance
# Assumes the dt curve is in standard units microseconds per foot and we want to convert to feet per second
L30_df['PVEL'] = (1/L30_df['DT'])*1e6
L30_df['PNII'] = L30_df['RHOB']*L30_df['PVEL']

# create 1D arrays of density and p-wave velocity from pandas dataframe
rhob = np.array(L30_df['RHOB'])
pvel = np.array(L30_df['PVEL'])

# Determine the index range for which both values are not nan's
# starting index
for i in range(0,len(L30_df.index)+1):
    if np.isnan(rhob[i]) == False and np.isnan(pvel[i]) == False:
        first_index = i
        break

# ending index
for j in range(0,len(L30_df.index)+1):
    if np.isnan(rhob[-j]) == False and np.isnan(pvel[-j]) == False:
        last_index = len(L30_df.index) - j
        break

rhobX = rhob[first_index:last_index:1]
pvelX = pvel[first_index:last_index:1]

# Create a crossplot of pvel vs rhob
#plt.plot(pvelX, rhobX, 'bx')
#plt.show()

# create new 1D arrays of the natural log of density and p-wave velocity
D = np.log(rhobX)
V = np.log(pvelX)

# create a crossplot of V vs D
#plt.plot(V, D, 'bx')
#plt.show()

# Setup the regression for A and B using scipy.stats
# slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
B, A, r_value, p_value, std_err = stats.linregress(V,D)

D_reg = A + B * V

# plot the regression
plt.plot(V, D, 'bx')
plt.plot(V, D_reg, 'rx')
plt.show()

# calculate alpha and beta from A and B
#alpha = np.exp(A)
#beta = B

# Estimate the density (rhobCALC) from velocity
#rhobCalc = alpha*(pvelX**beta)

# crossplot D and V and draw a line from the predicted values
#plt.plot(pvelX, rhobX, 'bx')
#plt.plot(pvelX, rhobCalc, 'r')
#plt.show()

# REFERENCES
# Castagna, J. P., and M. M. Backus, 1993, Offset-dependent reflectivity - Theory
# and practice of AVO analysis: SEG Investigations in Geophysics No. 8.
#
# Gardner,  G. H., Gardner, L. W., Gregory, A. R., 1974, Formation velocity and 
# density - the diagnostic basics for stratigraphic traps. GEOPHYSICS: 39, 770-780.
#
# Gregory, A. R., 1977, Aspects of rock physics from laboratory and log data 
# that are important to seismic interpretation.  Seismic Stratigraphy - 
# applications to hydrocarbon exploration: AAPG Memoir 26.
#
# 
#

 