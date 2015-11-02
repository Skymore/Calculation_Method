# -*- coding:utf-8 -*-
# -------------------------------------------------------
# Python PCHIP 分段三次Hermite插值
# Author: 陶睿 122345615@qq.com
# Date	: 2015-10-30
# version 1.2
# -------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# -------------------------------------------------------
# 函数名：pchip(X, Y, Yd, t)
# Piecewise Cubic Hermite Interpolating Polynomial
# 功能：分段三次Hermite插值算法
# 说明：
#	X为自变量取值向量。共n+1个值。X[0]...X[n]。
#	Y为对应X的函数值向量。共n+1个值。Y[0]...Y[n]。
#	Yd为对应X的导数值向量。共n+1个值。Yd[0]...Yd[n]。
#	t是一个值或向量。计算t处的插值结果。如果t是向量，返回一个插值结果的向量。t应满足 X[0]<=t<=X[n]。
#	插值多项式Pn(x)的次数为3次。
#	!需要import numpy。 X, Y, Yd为numpy.ndarray类型，t为numpy.float64类型或numpy.ndarray类型。
# 算法：
#	对每一段[Xi, Xi+1]调用Hermite_Interpolation算法。
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
	import numpy
	n = X.size - 1

	def _pchip(t):
		yi = 0
		pos = 0  # t在区间[X[pos] , x[pos+1]]中
		for i in range(0, n):
			if ((X[i] <= t) and (t <= X[i+1])):
				pos = i
		yi = Hermite_Interpolation(X[pos: pos+2], Y[pos: pos+2], Yd[pos: pos+2], t)
		return yi

	if (type(t) == numpy.float64 or type(t) == float):
		return _pchip(t)
	elif (type(t) == numpy.ndarray):
		P = np.zeros((t.size))
		for i in range(0, t.size):
			P[i] = _pchip(t[i])
		return P
# -------------------------------------------------------
# End of pchip(X, Y, Yd, t)
# -------------------------------------------------------

if __name__ == '__main__':

	y = lambda x: 1 / (1 + 25 * x**2)
	yd = lambda x: -50*x / (625 * x**4 + 50 * x**2 + 1)
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
		Yd = yd(X)
		H = lambda x: pchip(X, Y, Yd, x)
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
	ax2_1.plot(nList + 1, errList, linestyle = '-', color = 'blue', label = 'PCHIP')
	ax2_1.plot(nList + 1, errList, 'o', color = 'blue')
	ax2_1.grid(True)
	ax2_1.legend(loc = 'upper right')
	fig2.savefig(path + u'误差对比2.jpg')


	#图1：PCHIP
	n = 5
	X1 = np.linspace(a, b, n + 1)
	X2 = np.linspace(a, b, n + 2)
	X3 = np.linspace(a, b, n + 3)
	Y1 = y(X1)
	Y2 = y(X2)
	Y3 = y(X3)
	Yd1 = yd(X1)
	Yd2 = yd(X2)
	Yd3 = yd(X3)
	testF1 = pchip(X1, Y1, Yd1, testX)
	testF2 = pchip(X2, Y2, Yd2, testX)
	testF3 = pchip(X3, Y3, Yd3, testX)

	fig1 = plt.figure(31)
	plt.title('PCHIP   ' + r'$f(x) = \frac{1}{25x^2 + 1}$', fontsize = 15)
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
	fig1.savefig(path + u'PCHIP_1.jpg')

# 误差