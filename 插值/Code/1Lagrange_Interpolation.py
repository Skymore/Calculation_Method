# -*- coding:utf-8 -*-
# -------------------------------------------------------
# Python Lagrange-interpolation 拉格朗日插值
# Author: 陶睿 122345615@qq.com
# Date  : 2015-10-18
# version 1.2
# -------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate  # 求定积分函数

# -------------------------------------------------------
# 函数名：Lagrange_Interpolation(X, Y, Yd, t)
# 功能：拉格朗日插值算法
# 说明：
#	X为自变量取值向量。共n+1个值。X[0]...X[n]。
#	Y为对应X的函数值向量。共n+1个值。Y[0]...Y[n]。
#	t是一个值或向量。计算t处的插值结果。如果t是向量，返回一个插值结果的向量。t应满足 X[0]<=t<=X[n]
#	插值多项式Pn(x)的次数为n次。
#	!需要import numpy。 X, Y, Yd为numpy.ndarray类型，t为numpy.float64类型或numpy.ndarray类型。
# 算法：
#	Li(x) = 连乘j = 0..n, j!=i  (x - X[j]) / (X[i] - X[j])
#	Pn(x) = 累加i = 0..n, Li(x)
# -------------------------------------------------------
def Lagrange_Interpolation(X, Y, t):
	n = X.size - 1
	Pn = 0  # Pn: 拉格朗日插值多项式 Lagrange polynomial
	for i in range(0, n + 1):
		L = 1  # L: 拉格朗日插值基函数 Lagrange basis polynomials
		for j in range(0, n + 1):
			if (j != i):
				L *= (t - X[j]) / (X[i] - X[j])
		Pn += Y[i] * L
	return Pn


if __name__ == '__main__':

	y = lambda x: 1 / (1 + 25 * x**2)
	a = -1
	b = 1
	testX = np.linspace(a, b, 2001)
	testY = y(testX)
	path = u'/Users/sky/Desktop/计算方法/Interpolation/Figure/'  # 图片存储地址
	# 图3：Runge(龙格)现象，并记录误差err
	fig3 = plt.figure(13)
	plt.title(u'Runge(龙格)现象', fontsize = 15)
	ax3_1 = fig3.add_subplot(111)
	ax3_1.set_ylim(-2, 2)

	# 计算err随n的变化
	nMax = 40
	nBest = nMax
	errBest = 9999
	err = np.zeros(50)
	for n in range(nMax, 1, -1):  # n = nMax, nMax-1, ..., 2
		X = np.linspace(a, b, n + 1)
		Y = y(X)
		Pn = lambda x: Lagrange_Interpolation(X, Y, x)
		integrand = lambda x: abs(Pn(x) - y(x))
		err[n] = integrate.quad(integrand, a, b, limit = 2001)[0]
		#图3：Runge
		testF = Pn(testX)
		if (n <= 20):
			ax3_1.plot(testX, testF, color = 'gray', linestyle = '-', linewidth = 1)
		print u"#节点数 =", n + 1, u'次数 =', 2*n + 1, "err =", err[n]
		if(err[n] < errBest):
			nBest = n
			errBest = err[n]
	print 'nBest = ', nBest
	print 'errBest = ', errBest
	nList = np.nonzero(err)[0] 
	errList = np.log10(err[err!=0])

	#图3：Lagrange Runge
	ax3_1.grid(True)
	ax3_1.plot(testX, testY, color = 'red', linestyle = '-', linewidth = 2, label = u'原函数')
	fig3.savefig(path + u'Lagrange_Runge.jpg')

	#图4：err随节点数的变化(Lagrange VS Hermite)
	fig4 = plt.figure(14)
	plt.title(u'插值误差面积随节点数的变化', fontsize = 15)
	ax4_1 = fig4.add_subplot(111)
	ax4_1.set_xlabel(u'节点数')
	ax4_1.set_ylabel(u'[-1, 1]区间上误差面积')
	ax4_1.yaxis.set_ticks((-1, -0.52288, 0, 0.47712, 1, 1.47712, 2, 2.47712))
	ax4_1.yaxis.set_ticklabels(('0.1', '0.3', '1', '3', '10', '30', '100', '300'))
	ax4_1.set_ylim(-1.1, 2.9)
	ax4_1.set_xlim(0, 40)
	ax4_1.plot(nList + 1, errList, linestyle = '-', color = 'blue', label = 'Lagrange')
	ax4_1.plot(nList + 1, errList, 'o', color = 'blue')
	ax4_1.grid(True)
	ax4_1.legend(loc = 'upper right')
	fig4.savefig(path + u'误差对比1_2.jpg')

	#图2：err随多项式次数的变化(Lagrange VS Hermite)
	fig2 = plt.figure(12)
	plt.title(u'插值误差面积随多项式次数的变化', fontsize = 15)
	ax2_1 = fig2.add_subplot(111)
	ax2_1.set_xlabel(u'多项式次数')
	ax2_1.set_ylabel(u'[-1, 1]区间上误差面积')
	ax2_1.yaxis.set_ticks((-1, -0.52288, 0, 0.47712, 1, 1.47712, 2, 2.47712))
	ax2_1.yaxis.set_ticklabels(('0.1', '0.3', '1', '3', '10', '30', '100', '300'))
	ax2_1.set_ylim(-1.1, 2.9)
	ax2_1.set_xlim(0, 40)
	ax2_1.plot(nList, errList, linestyle = '-', color = 'blue', label = 'Lagrange')
	ax2_1.plot(nList, errList, 'o', color = 'blue')
	ax2_1.grid(True)
	ax2_1.legend(loc = 'upper right')
	fig2.savefig(path + u'误差对比1_1.jpg')


	#图1：Lagrange插值
	n = 5
	X1 = np.linspace(a, b, n + 1)
	X2 = np.linspace(a, b, n + 2)
	X3 = np.linspace(a, b, n + 3)
	Y1 = y(X1)
	Y2 = y(X2)
	Y3 = y(X3)
	testF1 = Lagrange_Interpolation(X1, Y1, testX)
	testF2 = Lagrange_Interpolation(X2, Y2, testX)
	testF3 = Lagrange_Interpolation(X3, Y3, testX)

	fig1 = plt.figure(11)
	plt.title('Lagrange Interpolation   ' + r'$f(x) = \frac{1}{25x^2 + 1}$', fontsize = 15)
	ax1_1 = fig1.add_subplot(111)
	ax1_1.set_ylim(-1, 1)
	ax1_1.plot(testX, testF1, color = 'b', linestyle = '-', linewidth = 1, label = u'%d个点'%(n + 1))
	ax1_1.plot(testX, testF2, color = 'r', linestyle = '-', linewidth = 1, label = u'%d个点'%(n + 2))
	ax1_1.plot(testX, testF3, color = 'g', linestyle = '-', linewidth = 1, label = u'%d个点'%(n + 3))
	ax1_1.plot(testX, testY, color = 'black', linestyle = '-', linewidth = 1, label = r'$f(x)$')
	ax1_1.plot(X1, Y1, 'o', color = 'b')
	ax1_1.plot(X2, Y2, 'o', color = 'r')
	ax1_1.plot(X3, Y3, 'o', color = 'g')
	ax1_1.grid(True)
	ax1_1.legend(loc = 'upper right')
	fig1.savefig(path + u'Lagrange_1.jpg')