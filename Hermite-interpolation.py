# -*- coding: utf-8 -*-

#-------------------------------------------------------
# Python Hermite-interpolation 埃尔米特插值
#Author	: Skymore 122345615@qq.com
#Date	: 2015-10-11
#-------------------------------------------------------
#H2n+1(x) = sum(j=0,n) [y[j]*Aj(x) + y'[j]*Bj(x)]
#Li(x) = 连乘j = 0..n,j!=i  (x - X[j]) / (X[i] - X[j])
#Ai(x) = [1 - 2(x-X[i])
#		*sum(j=0,n, j!=i) (1/(X[i]-X[j]))]*Li^2(x)
#Bi(x) = (x - X[i])*Li^2(x)
#-------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

#拉格朗日插值基函数
def L(i, t, X):
	Li = np.float64(1)
	for j in range(0, X.size):
		if(j != i):
			Li *= (t - X[j]) / (X[i] - X[j])
	return Li

#Li(t)在 t = X[i]的导数
def Ld(i, t, X):
	summ = np.float64(0)
	for j in range(0, X.size):
		if(j != i):
			summ += 1/(X[i] - X[j])
	return summ

def A(i, t, X):
	Ai = (1 - 2 * (t - X[i]) * Ld(i, t, X)) * (L(i, t, X)**2)
	return Ai

def B(i, t, X):
	Bi = (t - X[i]) * (L(i, t, X)**2)
	return Bi

def H(t, X, Y, Yd):
	ans = 0
	for i in range(0, X.size):
		ans += Y[i] * A(i, t, X) + Yd[i] * B(i, t, X)
	return ans

y = lambda x: 1 / (1 + 25 * x**2)
yd = lambda x: -50*x / (625 * x**4 + 50 * x**2 + 1)

if __name__ == '__main__':
	#X Y为生成插值函数的数据，testX testY为测试插值函数图像的数据
	#n为插值次数,X[i](i = 0...n)    X[0] = a, X[n] = b 等间距把[a,b]分为n份

	a = -10
	b = 10
	nBest = 2
	errBest = 9999
	for n in range(2, 21, 2):
		X = np.linspace(a, b, n + 1)
		Y = y(X)
		Yd = yd(X)

		#Hn为进行插值后得到的函数。
		Hn = lambda t: H(t, X, Y, Yd)

		integrand = lambda x: abs(Hn(x) - y(x))
		err, error = integrate.quad(integrand, a, b, limit = 101)
		print "n = ", n, "err = ", err
		if(err < errBest):
			nBest = n
			errBest = err
	print "nBest = ", nBest
	print "errBest = ", errBest

	n = nBest

	X = np.linspace(a, b, n + 1)
	Y = y(X)
	Yd = yd(X)
	Hn = lambda t: H(t, X, Y, Yd)
	testX = np.linspace(a,b,201)
	testY = y(testX)
	testYd = yd(testX)
	testF = H(testX,X,Y,Yd)

	Fig = plt.figure(1)
	#subplot(221)把绘图区域等分为2行*2列共4个区域,
	#然后在区域1(左上区域)中创建一个轴对象.
	#subplot(222)在区域2(左上区域)中创建一个轴对象.
	#subplot(212)把绘图区域等分为2行*1列共2个区域,
	#pl.subplot(212)在区域2(下区域)创建一个轴对象.
	ax1 = plt.subplot(311)
	ax2 = plt.subplot(312)
	ax3 = plt.subplot(313)
	
	#左上画f(x)及其取点和Hn(x)图像
	plt.sca(ax1)
	plt.title(u"埃尔米特插值 , n = %d, err = %.2f" % (nBest, errBest))
	plt.plot(X, Y, 'ro')
	plt.plot(testX, testY, color = "r", 
			linestyle = "-", label = "f(x)")
	plt.plot(testX, testF, color = "b", 
			linestyle = "-", label = "Hn(x)")
	plt.legend(loc='upper left')

	#右上画f'(x)及其取点
	plt.sca(ax2)
	plt.plot(X, Yd, 'ro')
	plt.plot(testX, testYd, color = "r", 
			linestyle = "-", label = "f'(x)")
	plt.plot(testX, 0*testY, color = "black", linestyle = "--")
	plt.legend(loc='upper left')

	#下方画Hn(x) - f(x)，即误差
	plt.sca(ax3)
	plt.plot(testX, testF - testY, color = "g", 
			linestyle = "-", label = "Pn(x) - f(x)")
	plt.plot(testX, 0*testY, color = "black", linestyle = "--")
	plt.legend(loc='upper left')

	#plt.show()
	Fig.savefig("Hermite-interpolation.pdf")