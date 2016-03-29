# -*- coding:utf-8 -*-
__author__ = 'Tao Rui'
__date__   = '2016-1-10'
import numpy as np
import scipy.linalg as lina
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

n = input('n = ') # 空间剖分数
m = input('m = ') # 时间剖分数
x = np.linspace(0, 1, n+1)
t = np.linspace(0, 1, m+1)
r = n * n / float(m) # 网比
A = np.diagflat(np.ones((1, n - 1), dtype = 'float64') * (1 - 2 * r))\
	+ np.diagflat(np.ones((1, n - 2), dtype = 'float64') * r, -1)\
	+ np.diagflat(np.ones((1, n - 2), dtype = 'float64') * r, 1)
#---------------------精确解求解--------------------
u = np.zeros((n + 1, m + 1), dtype = 'float64')
for i in range(0, n + 1):
	for j in range(0, m + 1):
		u[i,j] = np.exp(x[i] + t[j])
#--------------------初值条件求解-------------------
u1 = np.zeros((n + 1, m + 1), dtype = 'float64')

for i in range(0, n + 1):
	u1[i, 0] = np.exp(x[i])
for j in range(0, m + 1):
	u1[0, j] = np.exp(t[j])
for j in range(0, m + 1):
	u1[n, j] = np.exp(1 + t[j])
#---------------------数值解求解--------------------
for j in range(0, m):
	f = np.zeros((n));
	f[1] = r*u1[0,j];
	f[n - 1] = r*u1[n,j]
	u1[1:n, j+1] = np.dot(A, u1[1:n, j]) + f[1:n]
#--------------------画图------------------------
#fig = plt.figure()
#--------------------误差曲面图--------------------
#error = abs(u-u1)
#error = error.transpose()
#ax = fig.gca(projection='3d')
#X,T = np.meshgrid(x, t)
#surf = ax.plot_surface(T[1:50,:], X[1:50,:], error[1:50,:],rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)
#surf = ax.plot_surface(T, X, u1.transpose(),rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)
#--------------------误差曲线图---------------------
#ax = fig.gca()
E = abs(u[:,m]-u1[:,m])
ax.plot(x,E,'bo-');
plt.xlabel('x')
plt.xticks([0,0.2,0.4,0.6,0.8,1.0])
plt.ylabel(u'误差')
#plt.title(u't=1时误差曲线，n= %d m = %d'%(n,m))
#-------------------------------------------------
#plt.show()
#plt.savefig('1.jpg')