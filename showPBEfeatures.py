#! /usr/bin/env python
# -*- coding: utf-8 -*-

from aggBreakPBE import *

import numpy as np
import matplotlib.pyplot as plt


#-Plot with LaTeX
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#------------------------------Display Kernels---------------------------------#

v_list = np.arange(1,121)
D_F = 2.2
R_p = 1.16e-6
R_list = radiusOfGyration(v_list, D_F, R_p)
G = 1.0

#Brownian collision kernel:
T = 310.0
mu = 1.25e-3
D_list = EinsteinStokes(R_list, T, mu)

k_b = Brownian_kernel(R_list, D_list)

plt.pcolor(k_b)
plt.colorbar()
plt.xlabel('v_j')
plt.ylabel('v_i')
plt.show()

#Shear collision kernel:
k_s = shear_kernel(R_list, G)

plt.pcolor(k_s)
plt.colorbar()
plt.xlabel('v_j')
plt.ylabel('v_i')
plt.show()

#Breakup kernel: 
b = 2.
Gstar = 1000.0

c = 1.
k_b_c1 = breakup_kernel(R_list, G, Gstar, b, c)

c = 2.
k_b_c2 = breakup_kernel(R_list, G, Gstar, b, c)

c = 3.
k_b_c3 = breakup_kernel(R_list, G, Gstar, b, c)

plt.subplot(3,1,1)
plt.plot(v_list, k_b_c1, label='c = 1')
plt.legend(loc=2)
plt.xlabel('v_j')
plt.ylabel('k_b')

plt.subplot(3,1,2)
plt.plot(v_list, k_b_c2, label='c = 2')
plt.legend(loc=2)
plt.xlabel('v_j')
plt.ylabel('k_b')

plt.subplot(3,1,3)
plt.plot(v_list, k_b_c3, label='c = 3')
plt.legend(loc=2)
plt.xlabel('v_j')
plt.ylabel('k_b')
plt.show()

#---------------------Display Fragment Mass Distributions----------------------#

def showMatrix(fragDist, v_list):
  n_bin = fragDist.shape[0]
  fig = plt.figure(1, figsize=(5, 5))
  plt.clf()
  ax = fig.add_subplot(111)
  ax.set_aspect(1)
  ax.matshow(fragDist,cmap=plt.cm.binary)
  for i in xrange(n_bin):
    for j in xrange(n_bin):
      if fragDist[i,j] >= 1.5:
        color = "white"
      else:
        color = "black"
      ax.annotate("{:.2f}".format(fragDist[i][j]),
                                  xy=(j,i),
                                  horizontalalignment="center",
                                  verticalalignment="center",
                                  color=color,
                                  fontsize=10)
  plt.xticks(range(n_bin), [str(int(i)) for i in v_list])
  plt.yticks(range(n_bin), [str(int(i)) for i in v_list])
  plt.xlabel("v_j")
  plt.ylabel("v_i")

n_bin = 10
v_list = np.linspace(1, n_bin, n_bin)
#Binary splitting:
showMatrix(binaryFrag(n_bin), v_list)
plt.show()

#Normal splitting:
showMatrix(normalFrag(n_bin = 10, lambd = 4.), v_list)
plt.show()

#Uniform splitting:
showMatrix(uniformFrag(n_bin = 10), v_list)
plt.show()

#Erosion:
showMatrix(erosionFrag(n_bin = 10), v_list)
plt.show()

#-----------------------Display Interpolation Function-------------------------#

def showMatrix(fragDist, v_list):
  n_bin = fragDist.shape[0]
  fig = plt.figure(1, figsize=(5, 5))
  plt.clf()
  ax = fig.add_subplot(111)
  ax.set_aspect(1)
  ax.matshow(fragDist,cmap=plt.cm.binary)
  for i in xrange(n_bin):
    for j in xrange(n_bin):
      if fragDist[i,j] >= 0.75:
        color = "white"
      else:
        color = "black"
      ax.annotate("{:.2f}".format(fragDist[i][j]),
                                  xy=(j,i),
                                  horizontalalignment="center",
                                  verticalalignment="center",
                                  color=color,
                                  fontsize=10)
  plt.xticks(range(n_bin), [str(int(i)) for i in v_list])
  plt.yticks(range(n_bin), [str(int(i)) for i in v_list])
  plt.xlabel("v_j")
  plt.ylabel("v_i")
  plt.title("\chi_{ij}(v_i, v_j, v_k = 256)")


n_bin = 10
v_list = np.logspace(0., n_bin - 1, n_bin, base=2.0)

chi = interpolate(v_list)

showMatrix(chi[:, :, 8],v_list)
plt.show()



