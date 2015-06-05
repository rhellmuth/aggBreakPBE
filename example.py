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

abProps['saveFolder'] = "./output/"
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
abProps['D_F'] = 2.2                      # fractal dimension of aggregates
abProps['v_max'] = 2e3  #maximum number of monomers in an aggregate

#-Specify the bin structure you want to use: isGridUniform ?
# True  ---> v = [1, 2, 3, ..., vmax]
# False ---> v = [1, 2, 4, 8, ..., 2**i <= v_max]

#abProps['isGridUniform'] = True
abProps['isGridUniform'] = False

# PBE parameters:
tmin, tmax, tstep = (0., 60., 0.01)


abProps['t'] = np.arange(tmin, tmax, tstep) #time vector
abProps['G'] = 1000.                         #shear rate >= 0 (s-1)

# Aggregation parameters:
abProps['eta']        = 0.30                #aggregation efficiency
abProps['aggPhysics'] = [
                         'shear',
                         #'Brownian',
                         #'Sorensenian',
                         ]
abProps['T']  = 310.                       # (K) blood temperature
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

a = 1e11
b = 2.0
c = 3.0
Gstar = (a*R_p**c)**(-1./b)
#Gstar = 1000.

abProps['b'] = b
abProps['c'] = c
abProps['Gstar'] = Gstar


#---------------------Case 1: using above dictionary---------------------------# 
ab1 = aggBreak(abProps)
ab1.printSetup()
ab1.solve()
ab1.processData()
#ab1.saveOutput()

plt.plot(ab1.t, 100. * ab1.PA_ts, color='black')
plt.ylabel('Fraction of Aggregated Platelets, PA (\%)')
plt.xlabel('Time, t (s)')
plt.show()

plt.plot(ab1.t, ab1.v_mean_ts)
plt.ylabel('Mean Number of Platelets per Aggregate, \\textlangle v \\textrangle')
plt.xlabel('Time, t (s)')
plt.show()

for i in xrange(ab1.n_bin):
  plt.plot(ab1.t, 100.*ab1.PPD_list_ts[:,i], label='v = ' +str(ab1.v_list.astype(int)[i]))
plt.ylabel('Platelet Population Distribution, \\tilde{C} (\%)')
plt.xlabel('Time, t (s)')
plt.legend()
plt.show()

plt.loglog(ab1.v_list/ab1.v_mean_ts[-1], ab1.PPD_list_ts[-1], "ko-")
plt.xlabel('Normalized Number of Platelets per Aggregate, v / \\textlangle v \\textrangle')
plt.ylabel('Platelet Population Distribution, \\tilde{C} (\%)')
plt.ylim(ymin=1e-6)
plt.xlim(1e-3, 1e3)
plt.grid()
plt.show()

plt.plot(ab1.v_list, 100.*ab1.PPD_list_ts[-1], 'ko-')
plt.axis([1,6,0, 80])
plt.xlabel('Number of Platelets in Cluster')
plt.ylabel('Platelet Population Distribution, \\tilde{C} (\%)')
plt.show()

#-----------------Case 2: same theta, different conditions---------------------#

abProps2 = abProps.copy()
abProps2['G'] = 500.
a = 2e11
abProps2['Gstar'] = (a*R_p**c)**(-1./b)

ab2 = aggBreak(abProps2)
ab2.printSetup()
ab2.solve()
ab2.processData()
#ab2.saveOutput()

lbl1 = "G = {0:3.0f}, G* = {1:3.0f}, $\\theta$ = {2:.3g}".format(ab1.G,
                                                                ab1.Gstar,
                                                                ab1.theta)
lbl2 = "G = {0:3.0f}, G* = {1:3.0f}, $\\theta$ = {2:.3g}".format(ab2.G,
                                                                ab2.Gstar,
                                                                ab2.theta)
                                                                
plt.plot(ab1.t, 100. * ab1.PA_ts, color='blue', label=lbl1)
plt.plot(ab2.t, 100. * ab2.PA_ts, color='red', label=lbl2)
plt.legend(loc=4)
plt.ylabel('Fraction of Aggregated Platelets, PA (\%)')
plt.xlabel('Time, t (s)')
plt.show()

plt.plot(ab1.t, ab1.v_mean_ts, color='blue', label=lbl1)
plt.plot(ab2.t, ab2.v_mean_ts, color='red', label=lbl2)
plt.legend(loc=4)
plt.ylabel('Mean Number of Platelets per Aggregate, \\textlangle v \\textrangle')
plt.xlabel('Time, t (s)')
plt.show()

del abProps2
del ab2
del lbl2

#----------------Case 3: same conditions of case 1 + Diffusion-----------------#

abProps3B = abProps.copy()
abProps3B['aggPhysics']= ['Brownian', 
                         'shear',
                         ]

ab3B = aggBreak(abProps3B)
ab3B.printSetup()
ab3B.solve()
ab3B.processData()

abProps3S = abProps.copy()
abProps3S['aggPhysics']= ['Sorensenian', 
                         'shear',
                         ]

ab3S = aggBreak(abProps3S)
ab3S.printSetup()
ab3S.solve()
ab3S.processData()

lbl1 = "Shear alone, G = {0:3.0f}, Pe = {1:1.3f}".format(ab1.G, ab1.Pe)
lbl3B = "Shear + Brownian, G = {0:3.0f}, Pe = {1:1.3f}".format(ab3B.G, ab3B.Pe)
lbl3S = "Shear + Sorensenian, G = {0:3.0f}, Pe = {1:1.3f}".format(ab3S.G,
                                                                ab3S.Pe)

plt.plot(ab1.t, 100. * ab1.PA_ts, color='blue', label=lbl1)
plt.plot(ab3B.t, 100. * ab3B.PA_ts, color='red', label=lbl3B)
plt.plot(ab3S.t, 100. * ab3S.PA_ts, color='green', label=lbl3S)
plt.legend(loc=4)
plt.ylabel('Fraction of Aggregated Platelets, PA (\%)')
plt.xlabel('Time, t (s)')
plt.show()

plt.plot(ab1.t,  ab1.v_mean_ts, color='blue', label=lbl1)
plt.plot(ab3B.t, ab3B.v_mean_ts, color='red', label=lbl3B)
plt.plot(ab3S.t, ab3S.v_mean_ts, color='green', label=lbl3S)
plt.legend(loc=4)
plt.ylabel('Mean Number of Platelets per Aggregate, \\textlangle v \\textrangle')
plt.xlabel('Time, t (s)')
plt.show()

del abProps3B
del ab3B
del lbl3B

#----------Case 4: same theta comparing effect of Sorensenian collision--------#

abProps4S = abProps3S.copy()
abProps4S['G']= 100.
a = 1005640024739.8469
abProps4S['Gstar'] = (a*R_p**c)**(-1./b)

ab4S = aggBreak(abProps4S)
ab4S.printSetup()
ab4S.solve()
ab4S.processData()

abProps4lt1 = abProps3S.copy()
abProps4lt1['G']= 500.
a = 2.001e10
abProps4lt1['Gstar'] = (a*R_p**c)**(-1./b)

ab4lt1 = aggBreak(abProps4lt1)
ab4lt1.printSetup()
ab4lt1.solve()
ab4lt1.processData()

abProps4lt2 = abProps3S.copy()
abProps4lt2['G']= 1000.
a = 1e10
abProps4lt2['Gstar'] = (a*R_p**c)**(-1./b)

ab4lt2 = aggBreak(abProps4lt2)
ab4lt2.printSetup()
ab4lt2.solve()
ab4lt2.processData()


lbl3S = "G = {0:3.0f}, G* = {2:3.0f}, Pe = {1:1.3f}, $\\theta$ = {3:.3g}".format(ab3S.G, ab3S.Pe, ab3S.Gstar, ab3S.theta)
lbl4S = "G = {0:3.0f}, G* = {2:3.0f}, Pe = {1:1.3f}, $\\theta$ = {3:.3g}".format(ab4S.G, ab4S.Pe, ab4S.Gstar, ab4S.theta)
lbl4lt1 = "G = {0:3.0f}, G* = {2:3.0f}, Pe = {1:1.3f}, $\\theta$ = {3:.3g}".format(ab4lt1.G, ab4lt1.Pe, ab4lt1.Gstar, ab4lt1.theta)
lbl4lt2 = "G = {0:3.0f}, G* = {2:3.0f}, Pe = {1:1.3f}, $\\theta$ = {3:.3g}".format(ab4lt2.G, ab4lt2.Pe, ab4lt2.Gstar, ab4lt2.theta)


plt.plot(ab3S.t, 100. * ab3S.PA_ts, color='blue', label=lbl3S)
plt.plot(ab4S.t, 100. * ab4S.PA_ts, color='red', label=lbl4S)
plt.plot(ab4lt1.t, 100. * ab4lt1.PA_ts, color='green', label=lbl4lt1)
plt.plot(ab4S.t, 100. * ab4S.PA_ts, color='magenta', label=lbl4lt2)
plt.legend(loc=4)
plt.ylabel('Fraction of Aggregated Platelets, PA (\%)')
plt.xlabel('Time, t (s)')
plt.show()

plt.plot(ab3S.t,  ab3S.v_mean_ts, color='blue', label=lbl3S)
plt.plot(ab4S.t, ab4S.v_mean_ts, color='red', label=lbl4S)
plt.plot(ab4lt1.t,  ab4lt1.v_mean_ts, color='green', label=lbl4lt1)
plt.plot(ab4lt2.t, ab4lt2.v_mean_ts, color='magenta', label=lbl4lt2)
plt.legend(loc=5)
plt.ylabel('Mean Number of Platelets per Aggregate, \\textlangle v \\textrangle')
plt.xlabel('Time, t (s)')
plt.show()


