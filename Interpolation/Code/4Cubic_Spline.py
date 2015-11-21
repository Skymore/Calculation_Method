# -*- coding:utf-8 -*-
# -------------------------------------------------------
# Python Cubic Spline 三次样条插值
# Author: 陶睿 122345615@qq.com
# Date	: 2015-10-31
# version 1.2
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


if __name__ == '__main__':

	y = lambda x: 1 / (1 + 25 * x**2)
	a = -1
	b = 1
	testX = np.linspace(a, b, 2001)
	testY = y(testX)
	path = u'/Users/sky/Desktop/计算方法/Interpolation/Figure/'  # 图片存储地址
	
	# 计算err随n的变化
	nMax = 40
	nBest = nMax
	errBest = 9999
	err = np.zeros(50)
	for n in range(nMax, 1, -1):  # n = nMax, nMax-1, ..., 2
		X = np.linspace(a, b, n + 1)
		Y = y(X)
		H = lambda x: Cubic_Spline(X, Y, x)
		integrand = lambda x: abs(H(x) - y(x))
		err[n] = integrate.quad(integrand, a, b, limit = 2001)[0]
		print '#n = ', n, 'err = ', err[n]
		if(err[n] < errBest):
			nBest = n
			errBest = err[n]
	print 'nBest = ', nBest
	print 'errBest = ', errBest
	nList = np.nonzero(err)[0] 
	errList = np.log10(err[err!=0])

	# 图3：对散点插值对比(PCHIP VS Cubic_Spline)
	X = np.array([-3, -2, -1, 0, 1, 2, 3])
	Y = np.array([-1, -1, -1, 0, 1, 1, 1])
	Yd= np.array([ 0,  0,  0, 1, 0, 0, 0])
	testX2 = np.linspace(-3, 3, 2001)
	pc = pchip(X, Y, Yd, testX2)
	sp = Cubic_Spline(X, Y, testX2)
	fig3 = plt.figure(43)
	ax3_1 = fig3.add_subplot(111)
	ax3_1.set_ylim(-1.5, 1.5)
	ax3_1.plot(testX2, pc, color = 'g', label = 'PCHIP')
	ax3_1.plot(testX2, sp, color = 'r', label = 'Cubic Spline')
	ax3_1.plot(X, Y, 'o', color = 'b', label = 'data')
	ax3_1.grid(True)
	ax3_1.legend(loc = 'lower right')
	fig3.savefig(path + u'散点插值对比2.jpg')

	# 图2：err随节点数的变化(PCHIP VS Cubic_Spline)
	fig2 = plt.figure(32)
	plt.title(u'误差面积随节点数的变化', fontsize = 15)
	ax2_1 = fig2.add_subplot(111)
	ax2_1.set_xlabel(u'节点数')
	ax2_1.set_ylabel(u'[-1, 1]区间上误差面积')
	ax2_1.yaxis.set_ticks((-5, -4, -3, -2, -1, 0))
	ax2_1.yaxis.set_ticklabels(('1e-5', '1e-4', '1e-3', '0.01', '0.1', '1'))
	ax2_1.set_ylim(-5, 0)
	ax2_1.set_xlim(0, 40)	
	ax2_1.plot(nList + 1, errList, ls = '-', color = 'green', label = 'Cubic Spline')
	ax2_1.plot(nList + 1, errList, 'o', color = 'green')
	ax2_1.grid(True)
	ax2_1.legend(loc = 'upper right')
	fig2.savefig(path + u'误差对比2.jpg')

	#图1：Cubic_Spline
	n = 5
	X1 = np.linspace(a, b, n + 1)
	X2 = np.linspace(a, b, n + 2)
	X3 = np.linspace(a, b, n + 3)
	Y1 = y(X1)
	Y2 = y(X2)
	Y3 = y(X3)
	testF1 = Cubic_Spline(X1, Y1, testX)
	testF2 = Cubic_Spline(X2, Y2, testX)
	testF3 = Cubic_Spline(X3, Y3, testX)

	fig1 = plt.figure(41)
	plt.title('Cubic_Spline   ' + r'$f(x) = \frac{1}{25x^2 + 1}$', fontsize = 15)
	ax1_1 = fig1.add_subplot(111)
	ax1_1.set_ylim(-1, 1)
	ax1_1.plot(testX, testF1, color = 'b', linestyle = '-', linewidth = 1, label = u'%d个点'%(n+1))
	ax1_1.plot(testX, testF2, color = 'r', linestyle = '-', linewidth = 1, label = u'%d个点'%(n+2))
	ax1_1.plot(testX, testF3, color = 'g', linestyle = '-', linewidth = 1, label = u'%d个点'%(n+3))
	ax1_1.plot(testX, testY, color = 'black', linestyle = '-', linewidth = 1, label = r'$f(x)$')
	ax1_1.plot(X1, Y1, 'o', color = 'blue')
	ax1_1.plot(X2, Y2, 'o', color = 'red')
	ax1_1.plot(X3, Y3, 'o', color = 'green')
	ax1_1.grid(True)
	ax1_1.legend(loc = 'upper right')
	fig1.savefig(path + u'Cubic_Spline_1.jpg')