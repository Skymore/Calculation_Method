# -*- coding: utf-8 -*-

#-------------------------------------------------------
# Python Hermite-interpolation 埃尔米特插值
#Author	: Skymore 122345615@qq.com
#Date	: 2015-10-11
#-------------------------------------------------------
#H2n+1(x) = sum(j=0,n) [y[j]*Aj(x) + y'[j]*Bj(x)]
#Li(x) = 连乘j = 0..n,j!=i  (x - X[j]) / (X[i] - X[j])
#Ai(x) = [1 - 2(x-X[i])*Ld(i,x)*Li^2(x)
#Bi(x) = (x - X[i])*Li^2(x)

#设H(x) = sum(i=0...n) (Y[i]*Ai(x) + Yd[i]*Bi(x))
#Ai(X[i])  = delta(i,j) j=0..n
#Ai'(X[j]) = 0			j=0..n
#设Ai(x) = (a + b(x-X[i]))*Li^2(x)
#Ai


#Bi(X[j])  = 0 j=0...n
#Bi'(X[i]) = delta(i,j) j=0..n
#-------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

#拉格朗日插值基函数
def L(i, t, X):
	Li = 1
	for j in range(0, X.size):
		if(j != i):
			Li *= (t - X[j]) / (X[i] - X[j])
	return Li

#Li(t)在 t = X[i]的导数
Ld = lambda i, X: sum([ 1/(X[i] - X[j])
	for j in range(0, X.size) if (j != i) ])

A = lambda i, t, X:(1 - 2 * (t - X[i]) * Ld(i, X)) * (L(i, t, X)**2)
B = lambda i, t, X:(t - X[i]) * (L(i, t, X)**2)

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
	a = -3
	b = 3
	nBest = 2
	errBest = 9999
	for n in range(1, 21, 2):
		X = np.linspace(a, b, n + 1)
		Y = y(X)
		Yd = yd(X)

		#Hn为对X[i],Y[i]进行插值后得到的函数。
		Hn = lambda t: H(t, X, Y, Yd)

		integrand = lambda x: abs(Hn(x) - y(x))
		err, error = integrate.quad(integrand, a, b, limit = 1001)
		
		print "n = ", n, "err = ", err
		if(err < errBest):
			nBest = n
			errBest = err
	print "nBest = ", nBest
	print "errBest = ", errBest

	n = 10
	X = np.linspace(a, b, n + 1)
	Y = y(X)
	Yd = yd(X)
	Hn = lambda t: H(t, X, Y, Yd)

	#计算误差面积
	integrand = lambda x: abs(Hn(x) - y(x))
	err, error = integrate.quad(integrand, a, b, limit = 1001)
	
	#与原来函数进行比较
	testX = np.linspace(a,b,201)
	testY = y(testX)
	testYd = yd(testX)
	testF = H(testX,X,Y,Yd)

	#画图
	Fig = plt.figure(1)
	#subplot(221)把绘图区域等分为2行*2列共4个区域,
	#然后在区域1(左上区域)中创建一个轴对象.
	#subplot(222)在区域2(左上区域)中创建一个轴对象.
	#subplot(212)把绘图区域等分为2行*1列共2个区域,
	#pl.subplot(212)在区域2(下区域)创建一个轴对象.
	ax1 = plt.subplot(211)
	#ax2 = plt.subplot(312)
	ax3 = plt.subplot(212)
	
	#左上画f(x)及其取点和Hn(x)图像
	plt.sca(ax1)
	plt.title(u"Hermite, n = %d, err = %.3f" % (n, err))
	plt.ylim(-4, 4)
	plt.plot(X, Y, 'ro')
	plt.plot(testX, testY, color = "r", 
			linestyle = "-", label = "f(x)")
	plt.plot(testX, testF, color = "b", 
			linestyle = "-", label = "Hn(x)")
	plt.legend(loc='upper right')

	#右上画f'(x)及其取点
	#plt.sca(ax2)
	#plt.plot(X, Yd, 'ro')
	#plt.plot(testX, testYd, color = "r", 
	#		linestyle = "-", label = "f'(x)")
	#plt.plot(testX, 0*testY, color = "black", linestyle = "--")
	#plt.legend(loc='upper right')

	#下方画Hn(x) - f(x)，即误差
	plt.sca(ax3)
	plt.ylim(-4, 4)
	plt.plot(testX, testF - testY, color = "g", 
			linestyle = "-", label = "Pn(x) - f(x)")
	plt.plot(testX, 0*testY, color = "black", linestyle = "--")
	plt.legend(loc='upper right')

	plt.show()
	Fig.savefig("Hermite-interpolation.pdf")