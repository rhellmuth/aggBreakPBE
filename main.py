#! /usr/bin/env python
# -*- coding: utf-8 -*-

from aggBreakPBE import *

import numpy as np
import matplotlib.pyplot as plt


#-Plot with LaTeX
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


#---------------------------Simulation Inputs----------------------------------#

abProps = {} # input dictionary of aggregation-breakup properties

abProps['saveFolder'] = "./test/"
abProps['jobName'] = "test"
abProps['jobDate'] =   str(time.localtime().tm_year) + '-' \
          + str(time.localtime().tm_mon) + '-' \
          + str(time.localtime().tm_mday)
abProps['jobTime'] =   str(time.localtime().tm_hour) + ':' \
          + str(time.localtime().tm_min) + ':' \
          + str(time.localtime().tm_sec)

# CMD parameters:
abProps['C_p'] = 3.0e14                   # (m-3) absolute monomer concentration
R_p = 1.16e-6                             # (m) radius of each monomer
abProps['R_p'] = R_p
abProps['V_p'] = 4./3.*np.pi * R_p**3.    # (m3) volume of each monomer
abProps['D_F'] = 2.2                      # fractal dimension of aggregates
abProps['v_max'] = 2e3  #maximum number of monomers in an aggregate

#-Specify the bin structure you want to use: isGridUniform ?
# True  ---> v = [1, 2, 3, ..., vmax]
# False ---> v = [1, 2, 4, 8, ..., 2**i <= v_max]

#abProps['isGridUniform'] = True
abProps['isGridUniform'] = False

# PBE parameters:
tmin, tmax, tstep = (0., 30., 0.05)


abProps['t'] = np.arange(tmin, tmax, tstep) #time vector
abProps['G'] = 3350.                         #shear rate >= 0 (s-1)

# Aggregation parameters:
abProps['eta']        = 0.30                #aggregation efficiency
abProps['aggPhysics'] = [
                         'shear',
                         #'Brownian',
                         #'Sorensenian',
                         ]
abProps['T']  = 310.                       # (K)temperature of the system
abProps['mu'] = 1.25e-3                    # (Pa s) plasma viscosity

# Breakup Parameters:
abProps['isBreakupOn'] = True
#abProps['isBreakupOn'] = False

abProps['fragModel'] = 'binary'  #aggregate splits in half parts
#abProps['fragModel'] = 'normal'  #splits in parts that follow a normal PDF
#abProps['fragModel'] = 'uniform' #splits in parts with the same probability
#abProps['fragModel'] = 'erosion' #aggregate always loses 1 single monomer

# regulates the standard deviation of the normal splitting mass distribution
abProps['lambd'] = 1000 

a = 3.7e9
b = 2.0
c = 3.0
Gstar = (a*R_p**c)**(-1./b)
#Gstar = 1000.

abProps['b'] = b
abProps['c'] = c
abProps['Gstar'] = Gstar

ab = aggBreak(abProps)
ab.solve()
ab.saveOutput()

#------------------------------plot results------------------------------------#


plt.plot(ab.t, 100. * ab.PA_ts, color='black')
plt.ylabel('Fraction of Aggregated Platelets, PA (\%)')
plt.xlabel('Time, t (s)')
plt.show()

plt.plot(ab.t, ab.v_mean_ts)
plt.ylabel('Mean Number of Platelets per Aggregate, \\textlangle v \\textrangle')
plt.xlabel('Time, t (s)')
plt.show()

for i in xrange(ab.n_bin):
  plt.plot(ab.t, 100.*ab.PPD_list_ts[:,i], label='v = ' +str(ab.v_list.astype(int)[i]))
plt.ylabel('Platelet Population Distribution, \\tilde{C} (\%)')
plt.xlabel('Time, t (s)')
plt.legend()
plt.show()

plt.loglog(ab.v_list/ab.v_mean_ts[-1], ab.PPD_list_ts[-1], "ko-")
plt.xlabel('Normalized Number of Platelets per Aggregate, v / \\textlangle v \\textrangle')
plt.ylabel('Platelet Population Distribution, \\tilde{C} (\%)')
plt.ylim(ymin=1e-6)
plt.xlim(1e-3, 1e3)
plt.grid()
#plt.savefig('./self-similarDistribution.pdf', format='pdf')
plt.show()

plt.plot(ab.v_list, 100.*ab.PPD_list_ts[-1], 'ko-')
plt.axis([1,6,0, 80])
plt.xlabel('Number of Platelets in Cluster')
plt.ylabel('Platelet Population Distribution, \\tilde{C} (\%)')
#plt.savefig('PPD.pdf', format='pdf')
plt.show()



