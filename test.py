# -*- coding:utf-8 -*-
# -------------------------------------------------------
# Python Cubic Spline 三次样条插值
# Author: 陶睿 122345615@qq.com
# Date	: 2015-10-31
# version 1.0
# -------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# -------------------------------------------------------
# 函数名：Cubic_Spline(X, Y, t)
# Cubic Spline Interpolation
# 功能：三次样条插值算法(边界条件为端点处二阶微商已知且为0，即s''(x0) = 0, S''(xn) = 0)
# 说明：
#	X为自变量取值向量。共n+1个值。X[0]...X[n]。
#	Y为对应X的函数值向量。共n+1个值。Y[0]...Y[n]。
#	t是一个值或向量。计算t处的插值结果。如果t是向量，返回一个插值结果的向量。t应满足 X[0]<=t<=X[n]。
#	插值多项式Pn(x)的次数为3次。
#	!需要import numpy。 X, Y, Yd为numpy.ndarray类型，t为numpy.float64类型或numpy.ndarray类型。
# 算法：
#	计算出M向量后。再调用分段三次插值的算法(pchip(X, Y, Yd, t))。
# -------------------------------------------------------
def Hermite_Interpolation(X, Y, Yd, t):
	n = X.size - 1
	P = 0  # P:  2n+1次埃尔米特插值多项式 Hermite polynomial of degree 2n+1
	for i in range(0, n + 1):
		L = 1  # L: 拉格朗日插值基函数 Lagrange basis polynomials
		Ld = 0  # Ld: L在t = X[i]处的导数
		for j in range(0, n + 1):
			if (j != i):
				L *= (t - X[j]) / (X[i] - X[j])
				Ld += 1 / (X[i] - X[j])
		h = L**2 * (1 - 2 * (t - X[i]) * Ld)
		H = L**2 * (t - X[i])
		P += Y[i]*h + Yd[i]*H
	return P

def pchip(X, Y, Yd, t):
	import numpy as np
	n = X.size - 1

	def _pchip(t):
		yi = 0
		pos = 0  # t在区间[X[pos] , x[pos+1]]中
		for i in range(0, n):
			if ((X[i] <= t) and (t <= X[i+1])):
				pos = i
		yi = Hermite_Interpolation(X[pos: pos+2], Y[pos: pos+2], Yd[pos: pos+2], t)
		return yi

	if (type(t) == np.float64 or type(t) == float):
		return _pchip(t)
	elif (type(t) == np.ndarray):
		P = np.zeros((t.size))
		for i in range(0, t.size):
			P[i] = _pchip(t[i])
		return P
		
def Cubic_Spline(X, Y, t):
	import numpy as np
	n = X.size - 1
	h = np.zeros(X.size)
	h[0: n] = X[1: n+1] - X[0: n]
	dds0 = 0.  # S''(x0)
	ddsn = 0.  # S''(xn)
	A = np.zeros((n + 1))
	B = np.zeros((n + 1))
	beta = np.zeros((n + 1))
	alpha = np.zeros((n + 1))
	alpha[0] = 1.
	alpha[n] = 0.
	beta[0] = 3./h[0]   * (Y[1] - Y[0])   - h[0]/2.   * dds0
	beta[n] = 3./h[n-1] * (Y[n] - Y[n-1]) - h[n-1]/2. * ddsn
	for i in range(1, n):
		alpha[i] = h[i-1]/(h[i-1] + h[i])
		beta[i] = 3.*( (1-alpha[i])/h[i-1]*(Y[i] - Y[i-1]) \
					  +(  alpha[i])/h[i]  *(Y[i+1] - Y[i]) )
	#追赶法解三对角矩阵
	A[0] = -alpha[0]/2.
	B[0] = beta[0]/2.
	for i in range(1, n + 1):
		if (i < n):
			A[i] = -alpha[i]/(2.+(1-alpha[i])*A[i-1])   
			B[i] = (beta[i] - (1.-alpha[i])*B[i-1]) / (2. + (1.-alpha[i])*A[i-1])
	M = np.zeros((n + 1))
	M[n] = B[n]
	for i in range(n - 1, -1, -1):
		M[i] = A[i]*M[i+1]+B[i]
	#M = np.linalg.solve(A, beta)  #numpy.linalg.solve(A, B)是numpy中求解线性方程组的函数
	return pchip(X, Y, M, t)
# -------------------------------------------------------
# End of Cubic_Spline(X, Y, t)
# -------------------------------------------------------

X = np.array([-3, -2, -1, 0, 1, 2, 3])
Y = np.array([-1, -1, -1, 0, 1, 1, 1])
Yd= np.array([ 0,  0,  0, 1, 0, 0, 0])
testX = np.linspace(-3, 3, 1001)
pc = pchip(X, Y, Yd, testX)
sp = Cubic_Spline(X, Y, testX)
plt.plot(testX, pc, color = 'g')
plt.plot(testX, sp, color = 'r')
plt.plot(X, Y, 'o', color = 'b')
plt.ylim(-1.5, 1.5)
plt.show()