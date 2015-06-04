#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun 18 5 18:01:01 2014

@author: rudolf

_list    = vector array
_ts      = time series
_list_ts = time series of vector array
_pddf    = pandas.DataFrame()
_pds     = pandas.Series()
"""
import numpy as np
import scipy.integrate as si
from scipy import stats
import pandas as pd
import time

#==============================================================================#

class aggBreak(object):
  def __init__(self, dictionary):
    self.__dict__ = dictionary
    setCMD(self)
    setPBE(self)
    aggBreak.printSetup(self)
    self.isSolved = False
    self.isDataProcessed = False

  def solve(self):
    if self.isSolved:
      print "PBE has already been solved!"
    else:
      print "Solving..."
      tic = time.time()
      if self.isGridUniform:
	      self.C_list_ts = si.odeint(PBE_uniformGrid,
	                                 self.C0,
                                   self.t,
                                   args=(self,),
	                                 printmessg=1,
	                                 rtol=1e-10,
	                                 atol=1e-10)
      else:
	      self.C_list_ts = si.odeint(PBE_nonuniformGrid,
	                                 self.C0,
	                                 self.t,
	                                 args=(self,),
	                                 printmessg=1,
	                                 rtol=1e-10,
	                                 atol=1e-10)
      toc = time.time()
      runTime = toc - tic
      print "\tRun time = {0:3g} s\n".format(runTime)
      self.isSolved = True
      
      #check the quality of the simulation by testing mass conservation
      M1_ts = moment(self.C_list_ts, self.v_list, 1)
      conservation_tests = monomer_conservation(M1_ts)
      print "Initial and final total monomers in the box are {0:.2g} \
            and {1:.2g}, respectively.".format(
                                  conservation_tests[0], conservation_tests[1])
      print "Mass is conserved up to {0:.2f} %".format(
                                                    conservation_tests[2]*100.)
      
      print ""
      #check whether simulation hits the grid ceiling, whose results are unphysical.
      PPD_list_ts = PPD(self.C_list_ts, self.v_list)
      print "Relative platelet population in the last bin (grid ceiling)\
             is {0:.3f} %".format(PPD_list_ts[-1,self.n_bin -1]*100.)
      print ""

  def processData(self):
    if self.isSolved:
      print "Processing time series...\n"
      
      self.t_dimLess_s   = self.t / self.t_s
      self.t_dimLess_d   = self.t / self.t_d
      self.t_dimLess_b   = self.t / self.t_b
      self.G_ts          = self.G * np.ones_like(self.t)
      self.tau_ts        = self.mu * self.G_ts
      self.freeMono_ts   = free_monomers(self.C_list_ts)
      self.M0_ts         = moment(self.C_list_ts, self.v_list, 0)
      self.M1_ts         = moment(self.C_list_ts, self.v_list, 1)
      self.M2_ts         = moment(self.C_list_ts, self.v_list, 2)
      self.PPD_list_ts   = PPD(self.C_list_ts, self.v_list)
      self.monoC_list_ts = monomer_distribution(self.C_list_ts, self.v_list)
      self.R_mean_ts     = meanR(self.C_list_ts, self.v_list, self.R_list)
      self.R_rms_ts      = rmsR(self.C_list_ts, self.v_list, self.R_list)
      self.v_mean_ts     = meanV(self.C_list_ts, self.v_list)
      self.PA_ts         = PA(self.C_list_ts)
      self.dCdt_list_ts  = np.empty_like(self.C_list_ts)
      if self.isGridUniform:
        for i, (Ci, ti) in enumerate(zip(self.C_list_ts, self.t)):
          self.dCdt_list_ts[i] = PBE_uniformGrid(Ci, ti, self)
      else:
        for i, (Ci, ti) in enumerate(zip(self.C_list_ts, self.t)):
          self.dCdt_list_ts[i] = PBE_nonuniformGrid(Ci, ti, self)
          
      self.isDataProcessed = True
    else:
      print "Cannot process data before solving the PBE!\n"

  def printSetup(self):
    print "+System setup:"
    if self.isGridUniform:
      print "- Grid: UNIFORM"
    else:
      print "|- Grid: GEOMETRIC"
    print   "|  |- n_bin = {0:d}, v_max = {1:d}, R_p = {2:3g} m".format(
                                                                     self.n_bin,
                                                                     self.v_max,
                                                                     self.R_p)
    print "|"
    print "|- G = {0:2g} s-1".format(self.G)
    print "|"
    if self.isAggregationOn:
      print "|- Aggregation: ON"
      print '|  |- BROWNIAN: {0}, SORENSENIAN: {1}, SHEAR: {2}'.format(
                                                           self.isBrownAggOn,
                                                           self.isSorensenAggOn,
                                                           self.isShearAggOn)
    else:
      print "Aggregation: OFF!!!"
    print "|  |- t_s = {0:2g} s, t_d = {1:2g} s, Pe = {2:2g}".format(self.t_s,
                                                                       self.t_d,
                                                                       self.Pe)
    print "|"
    if self.isBreakupOn:
     print "|- Breakup: ON"
    else:
     print "|- Breakup: OFF"
    print "|  |- G* = {0:2g} s, b = {1:2g} s-1, c = {2:2g}".format(
                                                                   self.Gstar,
                                                                   self.b,
                                                                   self.c)
    print "|  |- t_b = {0:2g} s".format(self.t_b)
    print "|- theta = {0:2g}".format(self.theta)
    print ""

    
  def saveOutput(self):
    if self.isSolved:
      if not self.isDataProcessed:
        aggBreak.processData(self)
        
      filePrefix =   self.saveFolder \
                   + self.jobName + '-' \
                   + self.jobDate + '-' \
                   + self.jobTime + '-'
             
      #Input data:
      input_pds = pd.Series({'jobName':self.jobName,
                             'jobDate':self.jobDate,
                             'jobTime':self.jobTime,
                             'isGridUniform':self.isGridUniform,
                             'v_max':self.v_max,
                             'n_bin':self.n_bin,
                             'tmin':self.t[0],
                             'tmax':self.t[-1],
                             'tstep':self.t[1] - self.t[0],
                             'D_F':self.D_F,
                             'R_p':self.R_p,
                             'V_p':self.V_p,
                             'C_p':self.C_p,
                             'eta':self.eta,
                             'G':self.G,
                             'isBrownAggOn':self.isBrownAggOn,
                             'isShearAggOn':self.isShearAggOn,
                             'isSorensenAggOn':self.isSorensenAggOn,
                             'isAggregationOn':self.isAggregationOn,
                             'T':self.T,
                             'mu':self.mu,
                             't_d':self.t_d,
                             't_s':self.t_s,
                             'Pe':self.Pe,
                             't_a':self.t_a,
                             'isBreakupOn':self.isBreakupOn,
                             'b':self.b, 'c':self.c, 'Gstar':self.Gstar,
                             'lambd':self.lambd,
                             't_b':self.t_b,
                             'theta':self.theta})
                             
      inputFileName = filePrefix + 'inputParameters.csv'
      print "Saving input in " + inputFileName
      input_pds.to_csv(inputFileName)

      # v_list and R_list as a table:
      vR_list_pddf = pd.DataFrame()
      vR_list_pddf['v'] = self.v_list.astype(int)
      vR_list_pddf['R'] = self.R_list
                                   
      vR_listFileName = filePrefix + 'vR_lists.csv'
      print "Saving cluster sizes in " + vR_listFileName
      vR_list_pddf.to_csv(vR_listFileName)
      
      #Fragment mass distribution:
      fragDist_pddf = pd.DataFrame(self.fragDistr)
      
      fragFileName = filePrefix + 'fragDist.csv'
      print "Saving fragment distribution in " + fragFileName
      fragDist_pddf.to_csv(fragFileName)

      #Output time series (ts):
      C_nameList = []
      for i in xrange(self.v_list.shape[0]):
        C_nameList.append('C_v[{0:d}]'.format(i))

      output_pddf_ts = pd.DataFrame(self.C_list_ts,columns=C_nameList)
      output_pddf_ts['t']           = self.t
      output_pddf_ts['t_dimLess_s'] = self.t_dimLess_s
      output_pddf_ts['t_dimLess_d'] = self.t_dimLess_d
      output_pddf_ts['t_dimLess_b'] = self.t_dimLess_b
      output_pddf_ts['G']           = self.G_ts
      output_pddf_ts['tau']         = self.tau_ts
      output_pddf_ts['freeMono']    = self.freeMono_ts
      output_pddf_ts['M_0']         = self.M0_ts
      output_pddf_ts['M_1']         = self.M1_ts
      output_pddf_ts['M_2']         = self.M2_ts
      output_pddf_ts['R_mean']      = self.R_mean_ts
      output_pddf_ts['R_rms']       = self.R_rms_ts
      output_pddf_ts['v_mean']      = self.v_mean_ts
      output_pddf_ts['PA']          = self.PA_ts

      for i in xrange(self.v_list.shape[0]):
        output_pddf_ts['PPD_v[{0:d}]'.format(i)] = self.PPD_list_ts[:,i]
      for i in xrange(self.v_list.shape[0]):
        output_pddf_ts['dCdt_v[{0:d}]'.format(i)]  = self.dCdt_list_ts[:,i]

      outputFileName = filePrefix + 'output_timeSeries.csv'
      print "Saving output time series in " + outputFileName
      output_pddf_ts.to_csv(outputFileName)


#---------------------system initialization functions--------------------------#

def setCMD(self):
  """ n_bin:  number of bins  i = [0, 1, 2, ..., n]
      chi:    interpolation operator
      R_list: cluster radius of gyration
      v_list: number of platelets in the cluster"""
  
  if self.isGridUniform:
    #list of number of monomers in the cluster
    self.v_list = np.linspace(1.,
                             self.v_max,
                             self.v_max)
    self.n_bin  = int(len(self.v_list))
    self.chi = []
  else:
    #geometric grid
    self.n_bin = int(np.floor(np.log2(self.v_max)) + 1)
    self.v_list = np.logspace(0.,
                             self.n_bin - 1,
                             self.n_bin,
                             base=2.0)
    self.chi = interpolate(self.v_list)
  
  self.v_max = int(self.v_list[-1])         #v_max label gets truncated
  self.R_list = radiusOfGyration(self.v_list, self.D_F, self.R_p)


def setPBE(self):
  self.D_list = []                               # diffusivity of each cluster

  self.k_c = np.zeros([self.n_bin, self.n_bin])  # Smoluchowski collision kernel
  self.k_d = np.zeros([self.n_bin, self.n_bin])  # Diffusion collision kernel
  self.k_s = np.zeros([self.n_bin, self.n_bin])  # Shear collision kernel
  self.k_b = np.zeros(self.n_bin)                # Breakup kernel

  self.t_a = 1e30
  self.t_d = 1e30
  self.t_s = 1e30
  self.t_b = 1e30

  self.isAggregationOn = False
  self.isShearAggOn    = False
  self.isSorensenAggOn = False
  self.isBrownAggOn    = False
  
  #checks the given aggregation physics and sums the collision kernels into k_c
  if 'Brownian' in self.aggPhysics or 'Sorensenian' in self.aggPhysics:
    self.isBrownAggOn    = True
    self.isAggregationOn = True
    self.D_list = EinsteinStokes(self.R_list, self.T, self.mu)
    if 'Sorensenian' in self.aggPhysics:
      self.isSorensenAggOn = True
      self.k_d = Sorensenian_kernel(self.R_list, self.D_list, self.G)
      self.t_d = 1.0 / (  4.0*np.pi \
                        * 2.0 * self.D_list[0] \
                        * 2.0 * self.R_p \
                        * self.C_p) 
    else:
      self.k_d = Brownian_kernel(self.R_list, self.D_list)
      self.t_d = 1.0 / (  (1.0 + 1.05*self.G) \
                        * 4.0*np.pi \
                        * 2.0 * self.D_list[0] \
                        * 2.0 * self.R_p \
                        * self.C_p) 

  if 'shear'in self.aggPhysics:
    self.isAggregationOn = True
    self.isShearAggOn    = True
    self.k_s = shear_kernel(self.R_list, self.G)
    self.t_s = 1.0 / (4./3. * self.G * (2.0 * self.R_p)**3 * self.C_p)
    
  self.k_c = self.k_s + self.k_d
#  else:
#    print   "No valid aggregation physics was given! Aggregation is OFF!" \
#          + " Try: 'shear', 'Brownian', or 'Sorensenian'."

  if self.isBreakupOn:
    self.k_b = breakup_kernel(self.R_list, self.G, self.Gstar, self.b, self.c)
    if self.fragModel is 'binary':
      if self.isGridUniform:
        self.fragDistr = binaryFrag(self.n_bin)
      else:
        self.fragDistr = nonunifBinaryFrag(self.n_bin)
#    elif self.fragModel is 'normalFrag':
#      if self.isGridUniform:
#        self.fragDistr = normalFrag()
#      else:
#        self.fragDistr = nonunifNormalFrag()
#    elif self.fragModel is 'uniform':
#      if self.isGridUniform:
#        self.fragDistr = uniformFrag()
#      else:
#        self.fragDistr = nonunifUniformFrag()
#    elif self.fragModel is 'erosion':
#      if self.isGridUniform:
#        self.fragDistr = erosionFrag()
#      else:
#        self.fragDistr = nonunifErosionFrag()
      self.t_b = 1.0 / ((self.G / self.Gstar)**self.b)
    else:
      print "No valid fragmentation model was given! Turning breakup OFF!"
      self.isBreakupOn = False 
      self.fragDistr = np.zeros((self.n_bin, self.n_bin))
  else:
    self.k_b = breakup_kernel(self.R_list, 0.0, 1.0, 1.0, 1.0)
    self.fragDistr = np.zeros((self.n_bin, self.n_bin))
  
  #generate initial condition (no aggregates, just monomers)
  self.C0    = np.zeros(self.n_bin)
  self.C0[0] = self.C_p
  
  self.Pe    = self.t_s / self.t_d
  self.t_a   = 1.0 / (self.eta * (1.0 / self.t_s + 1.0 / self.t_d))
  self.theta = self.t_a / self.t_b

#-----------------------------other functions----------------------------------#

def EinsteinStokes(R, T, mu):
  """Diffusivity constant using Einsteins-Stokes equation."""
  k_B = 1.386488e-23  # (J K-1) Boltzmann constant  
  return (k_B * T) / (6.0 * np.pi * mu * R)


def radiusOfGyration(v, D_F, R_p):
  """Radius of gyration of a fractal aggregate ball."""
  return v**(1.0 / D_F) * R_p
  

def interpolate(v_list):
  """Interpolation function of Kumar & Ramkrishna (1997), also presented by \
Garrick, Lehtinen & Zachariah (2006)."""
  n_bin  = int(len(v_list))
  chi = np.zeros((n_bin, n_bin, n_bin))
  # Vsum_ij = (v_i + v_j) is a matrix 
  Vsum = v_list[:,np.newaxis] + v_list[np.newaxis,:]
  for k in xrange(1, n_bin): # no aggregation to the last bin
    if k < n_bin -1:
      chi[:,:,k] = np.where((Vsum <= v_list[k+1]) & (Vsum >= v_list[k]),
	                          (v_list[k+1] - Vsum) / (v_list[k+1] - v_list[k]),
	                          chi[:,:,k] )
    chi[:,:,k] = np.where((Vsum <= v_list[k]) & (Vsum >= v_list[k-1]),
		                      (Vsum - v_list[k-1]) / (v_list[k] - v_list[k-1]),
		                      chi[:,:,k])
  return chi
        
#----------------------------kernel functions----------------------------------#


def Brownian_kernel(R_list, D_list):
  k_d = 4.0 * np.pi \
	      * (R_list[:,np.newaxis] + R_list[np.newaxis,:]) \
	      * (D_list[:,np.newaxis] + D_list[np.newaxis,:])
  return removeIndexN(k_d)

def Sorensenian_kernel(R_list, D_list, G):
	return (1.0 + 1.05 * G) * Brownian_kernel(R_list, D_list)
	
def shear_kernel(R_list, G):
  k_s = 4.0 / 3.0 * G * (R_list[:,np.newaxis] + R_list[np.newaxis,:])**3.0
  return removeIndexN(k_s)
	
def breakup_kernel(R_list, G, Gstar, b, c):
  kern_mat = (G / Gstar)**b * (R_list[:] / R_list[0])**c
  kern_mat[0] = 0.0 # monomer doesn't break
  return kern_mat
  
def removeIndexN(k):
  """Remove i = n from Smoluchowski collision kernel. The aggregation with the\
  largest aggregate causes mass loss."""
  n_bin = k.shape[0]
  for i in xrange(n_bin):
    k[i, n_bin-1] = 0.0
    k[n_bin-1, i] = 0.0
  return k

#-----------------------fragment mass distributions----------------------------#

def erosionFrag(n_bin):
  massDistr = np.zeros((n_bin, n_bin))
  for j in xrange(n_bin):
    p = np.zeros(n_bin)
    for i in xrange(np.int((j+1)/2)):
      p[0] = 1.0
    for i in xrange(n_bin):
      massDistr[i][j] = p[i] + p[j-i-1]
  return massDistr
    
#def nonunifErosionFrag(n_bin):
#  massDistr = np.zeros((n_bin, n_bin))
#  for j in xrange(n_bin):
#    p = np.zeros(n_bin)
#    for i in xrange(np.int((j+1)/2)):
#      p[0] = 1.0
#    for i in xrange(n_bin):
#      massDistr[i][j] = 
#  return massDistr
    
def binaryFrag(n_bin):
  massDistr = np.zeros((n_bin, n_bin))
  for j in xrange(n_bin):
    p = np.zeros(n_bin)
    for i in xrange(np.int((j+1)/2)):
      p[np.int((j)/2)] = 1.0
    for i in xrange(n_bin):
      massDistr[i][j] = p[i] + p[j-i-1]
  return massDistr

def nonunifBinaryFrag(n_bin):
  massDistr = np.zeros((n_bin, n_bin))
  for j in xrange(1, n_bin):
    for i in xrange(0, n_bin):
      if i == j - 1:
        massDistr[i,j] = 2.0
  return massDistr
    
def normalFrag(n_bin, lambd = 1):
    massDistr = np.empty((n_bin,n_bin))
    for j in xrange(n_bin):
        p = np.zeros(n_bin)
        
        std = (j+1)/2 / float(lambd)
        mean = (j+1)/2
        N = stats.norm(loc=mean, scale=std)

        for i in xrange(np.int((j+1)/2)):
          p[i] = N.cdf(i+1) - N.cdf(i)

        # to regulate the distribution to the probability p[i<0] = 1 and sum(p[i>0]) = 1      
        sump = p.sum()
        if sump != 0:
          par = 1.0 / sump
        else:
          par = 2.0
 
        for i in xrange(n_bin):
          massDistr[i][j] =  par * (p[i] + p[j-i-1])

    return massDistr

#def nonunifNormalFrag():
#    massDistr = np.zeros((n_bin, n_bin))   
#    nIntP = 500 #number of integration points
#    for j in xrange(1, n_bin):
#        p = np.zeros(n_bin)
#        std = (v_list[j] / 2) / float(lambd)
#        mean = v_list[j] / 2
#        N = stats.norm(loc=mean, scale=std)
#        
#        for i in xrange(j+1):
#          if i > 0: #don't do this if i == 0
#            #chi is the interpolation vector
#            chi = np.linspace(0.,
#                              1.,
#                              nIntP)
#            ps  = 2*N.pdf(np.linspace(v_list[i-1],
#                                    v_list[i],
#                                    nIntP))
#            x = np.linspace(v_list[i-1], v_list[i], nIntP)
#            #Simpson integration of vector chi * ps
#            p[i] += si.simps(chi * ps, x)

#          if i < n_bin-1: #don't do this is i == n_bin
#            chi = np.linspace(1.,
#                              0.,
#                              nIntP)
#            ps  = 2*N.pdf(np.linspace(v_list[i],
#                                    v_list[i+1],
#                                    nIntP))
#            x = np.linspace(v_list[i], v_list[i+1], nIntP)
#            p[i] += si.simps(chi * ps, x)

##          if v_list[i] == v_list[j]/2: #if it splits in half -> two daughters
##            p[i] += 2*(N.cdf(v_list[i] + 1) - N.cdf(v_list[i] - 1)) / 2.0
#        # to regulate the distribution to the probability p[i<0] = and sum(p[i>0]) = 1      
#        sump = p.sum()
#        print sump
#        if sump != 0:
#          par = 2.0 / sump
#        for i in xrange(n_bin):
#          massDistr[i][j] =  par * p[i]
#    return massDistr


def uniformFrag(n_bin):
    massDistr = np.empty((n_bin, n_bin))
    for j in xrange(n_bin):
        p = np.zeros(n_bin)
        for i in xrange(np.int((j+1)/2)):
          p[i] = 1.0 / np.floor((j+1)/2)
        for i in xrange(n_bin):
          massDistr[i][j] = p[i] + p[j-i-1]   
    return massDistr

#def unifDistr(v_j, intPoints): #used in nonunifUniformFrag()
#  p = 2.0 / v_j
#  ps = p * np.ones_like(intPoints)
#  for i in xrange(np.shape(intPoints)[0]):
#    if intPoints[i] > v_j:
#      ps[i] = 0
#  return ps

#def nonunifUniformFrag():
#    massDistr = np.zeros((n_bin, n_bin))
#    nIntP = 500 #number of integration points
#    for j in xrange(1, n_bin):
#        p = np.zeros(n_bin)
#        unif = 2.0 / v_list[j]
#        for i in xrange(j+1):
#          if i > 0: #don't do this if i == 0
#            #chi is the interpolation vector
#            chi = np.linspace(0.,
#                              1.,
#                              nIntP)
##            ps  = unif * np.ones(nIntP)
#            x = np.linspace(v_list[i-1], v_list[i], nIntP)
#            ps = unifDistr(v_list[j], x)
#            #Simpson integration of vector chi * ps
#            p[i] += si.simps(chi * ps, x)

#          if i < n_bin-1: #don't do this is i == n_bin
#            chi = np.linspace(1.,
#                              0.,
#                              nIntP)
##            ps  = unif * np.ones(nIntP)
#            x = np.linspace(v_list[i], v_list[i+1], nIntP)
#            ps = unifDistr(v_list[j], x)            
#            p[i] += si.simps(chi * ps, x)
#          if v_list[i] == v_list[j]/2:
#            p[i] += unif
#        
##        sump = p.sum()
##        if sump != 0:
##          par = 2.0 / sump
#        par = 1
#        for i in xrange(n_bin):
#          massDistr[i][j] = par * p[i]
#    return massDistr

#---------------------------------PBE's----------------------------------------#

def PBE_uniformGrid(y, t, self):
  aggregation = np.zeros(self.n_bin)
  breakage    = np.zeros(self.n_bin)

  # Smoluchowski (1917)
  if self.isAggregationOn:
    creation    = np.zeros(self.n_bin)
    destruction = np.zeros(self.n_bin)

    destruction = -np.dot(np.transpose(self.k_c), y) * y 

    kyn = self.k_c * y[:,np.newaxis] * y[np.newaxis,:]
    for i in xrange(self.n_bin):
      creation[i] = np.sum(kyn[np.arange(i), i - np.arange(i) - 1])
    creation = 0.5 * creation

    aggregation = self.eta * (creation + destruction)
	
  # Pandya & Spielman (1982)
  if self.isBreakupOn:
    creation    = np.zeros(self.n_bin)
    destruction = np.zeros(self.n_bin)

    destruction = -self.k_b * y
    creation = np.dot(self.fragDistr, self.k_b * y)

    breakage = creation + destruction

  return aggregation + breakage


def PBE_nonuniformGrid(y, t, self):
  aggregation = np.zeros(self.n_bin)
  breakage    = np.zeros(self.n_bin)

  # aggregation PBE by Smoluchowski (1917)
  creation    = np.zeros(self.n_bin)
  destruction = np.zeros(self.n_bin)

  if self.isAggregationOn:
    for i in xrange(self.n_bin):
       creation[i] = np.dot(np.dot(self.k_c * self.chi[:, :, i], y), y)
    creation = 0.5 * creation

    destruction = -np.dot(self.k_c, y) * y  

    aggregation = self.eta * (creation + destruction)
	
  # breakage PBE by Pandya & Spielman (1982)
  if self.isBreakupOn:
    creation    = np.zeros(self.n_bin)
    destruction = np.zeros(self.n_bin)

    destruction = -self.k_b * y
    creation    = np.dot(self.fragDistr, self.k_b * y)

    breakage = creation + destruction
	
  return aggregation + breakage


#--------------------------Data process functions------------------------------#

def moment(C_list_ts, v_list, order):
  """Integral moment of the CMD."""
  out = np.sum(C_list_ts * v_list**order, axis=1)
  return out

def monomer_conservation(M1_ts):
  """Evalluates mass concervation. 
     Returns [initial monomer, final monomer, final / initial ratio]."""
  ini_monomerC = M1_ts[0]
  fin_monomerC = M1_ts[-1]
  conservation = fin_monomerC / ini_monomerC

  results = np.array([ini_monomerC, fin_monomerC, conservation])
  return results

def free_monomers(C_list_ts):
  """Concentration of free platelets."""
  return C_list_ts[:,0]

def meanV(C_list_ts, v_list):
  """Mean number of platelets per aggregate."""
  M0_ts = moment(C_list_ts, v_list, 0)
  M1_ts = moment(C_list_ts, v_list, 1)
  return M1_ts / M0_ts
         
def meanR(C_list_ts, v_list, R_list):
  """Mean radius of gyration."""
  numerator   = np.sum(C_list_ts * R_list, axis=1) 
  denominator = moment(C_list_ts, v_list, 0)
  return numerator / denominator

def rmsR(C_list_ts, v_list, R_list):
  """Root-mean-square radius of gyration"""
  numerator = np.sum(C_list_ts * (R_list * v_list)**2, axis=1)
  denominator = moment(C_list_ts, v_list, 2)
  return np.sqrt(numerator / denominator)

def monomer_distribution(C_list_ts, v_list):
  """Calculates the distribution of monomers in cluster."""
  out = C_list_ts * v_list
  return out
   
def PA(C_list_ts):
  """Fraction of aggregated platelets."""
  N_t = free_monomers(C_list_ts)
  return (1 - N_t/N_t[0])
   
def PPD(C_list_ts, v_list):
  """Platelet population distribution"""
  n_bin = v_list.shape[0]
  M1_ts = moment(C_list_ts, v_list, 1)
  return C_list_ts * v_list / np.repeat(M1_ts[:,np.newaxis], n_bin, axis=1)

#def fitting_stat(x, y):
#   slope, intercept, r, prob2, see = stats.linregress(x, y)

#   if len(x) > 2:
#       see=see*np.sqrt(len(x)/(len(x)-2.))

#       mx = x.mean()
#       sx2 = ((x-mx)**2).sum()
#       sd_intercept = see * np.sqrt(1./len(x) + mx*mx/sx2)
#       sd_slope = see * np.sqrt(1./sx2)

#   results=np.zeros(5)

#   results[0]=slope
#   results[1]=intercept
#   results[2]=r
#   if len(x) > 2:
#       results[3]=sd_slope
#       results[4]=sd_intercept

#   return results


