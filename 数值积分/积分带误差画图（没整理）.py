# -*- coding:utf-8 -*-
# -------------------------------------------------------
# Python integrate
# Author: 陶睿 122345615@qq.com
# Date	: 2015-12-21
# version 1.0
# -------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
path = u'/Users/sky/Desktop/计算方法/'
# -------------------------------------------------------
# 函数名:Newton_Cotes_integrate(f, a, b, n)
# 功能:Newton-Cotes公式
# 说明:
# Parameters:
#	(a,b)积分区间。
#	f:function object 被积函数
#	n:等分数量 X[0]...X[n] 需满足n<=7.
#		n=8的时候由于龙格现象会出现数值不稳定
# Returns:
#	res:一个浮点数,Y在整个区间(a,b)上的积分
# -------------------------------------------------------
def Newton_Cotes_integrate(f, a, b, n):
	import numpy as np
	h = (b-a)/(1.*n)
	X = np.zeros(n + 1)
	for i in range(0, n+1):
		X[i] = a + i*h
	C = np.zeros((8, 8)) #C[n][i] 即书中C^{(n)}_{i}
	C = [[]]
	C.append([1./2,  1./2])
	C.append([1./6,  2./3,  1./6])
	C.append([1./8,  3./8,  3./8,  1./8])
	C.append([7./90,  16./45,  2./15,  16./45,  7./90])
	C.append([19./288,  25./96,  25./144,  25./144,  25./96,  19./288])
	C.append([41./840,  9./35,  9./280,  34./105,  9./280,  9./35,  41./840])
	C.append([751./17280,  3577./17280,  1323./17280,  2989./17280,  \
			2989./17280,  1323./17280,  3577./17280,  751./17280])
	res = sum((b-a)*C[n][i]*f(X[i]) for i in range(0, n+1))
	return res

# -------------------------------------------------------
# 函数名:integrate_cumtrapz(Y, X)
# 功能:复化梯形公式求累积积分
# 说明:
# Parameters:
#	Y:array_like 为对应X的函数值向量。共n+1个值。Y[0]...Y[n]。
#	X:array_like 为自变量取值向量。共n+1个值。X[0]...X[n]。等距节点或者非等距节点均可
# Returns:
#	res:ndarray The result of cumulative integration of y.
#		y的累积积分。res[0](=0)...res[n]
# -------------------------------------------------------
def integrate_cumtrapz(Y, X):
	import numpy as np
	n = Y.size - 1 #n = 19
	h = np.zeros(n + 1);
	h[0:n] = X[1:n+1] - X[0:n] # h[i] = X[i+1] - X[i]
	res = np.zeros(n + 1)
	for i in range(1, n+1):
		res[i] = (Y[i-1]+Y[i])*h[i-1]/2.0 + res[i-1]
	return res
def c_integrate(f, a, b, n):
	import numpy as np
	X = np.linspace(a, b, n + 1)
	Y = f(X)
	return integrate_cumtrapz(Y,X)[-1]
# -------------------------------------------------------
# 函数名:integrate_parabola(Y, a, b, n)
# 功能:复化抛物线公式求积分
# 说明:
# Parameters:
#	(a,b)积分区间。
#	n:把(a,b)区间等分为2n份。断点为X[0]...X[2n]
#	如n等于3，则区间分为6份。X[0]...X[7]
#	Y:array_like 为对应每个X的函数值向量。共2n+1个值。Y[0]...Y[2n]。
# Returns:
#	res:一个浮点数,Y在整个区间(a,b)上的积分
# -------------------------------------------------------
def integrate_parabola(Y, a, b, n):
	import numpy as np
	h = (b - a)/(2. * n)
	sum1 = sum(Y[2*k-1] for k in range(1, n+1))
	sum2 = sum(Y[2*k] for k in range(1, n))
	res = h/3.*(Y[0]+Y[2*n]+4.*sum1+2.*sum2)
	return res
def p_integrate(f, a, b, n):
	import numpy as np
	X = np.linspace(a,b,2*n+1)
	Y = f(X)
	return integrate_parabola(Y,a,b,n)

# -------------------------------------------------------
# 函数名:half_integrate_trapezoid(f, a, b, eps)
# 功能:梯形公式的逐次分半法求积分
# 说明:
# Parameters:
#	(a,b)积分区间。
#	f:function object 被积函数
#	eps最大允许误差
# Returns:
#	res:一个浮点数,Y在整个区间(a,b)上的积分
# -------------------------------------------------------
def half_integrate_trapezoid(f, a, b, eps):
	import numpy as np
	n = 1;
	h = (b-a)/2.
	Tn = h*(f(a)+f(b))
	T2n = 0. #Tn为书中T0 T2n为书中T
	while(True):
		F = sum(f(a + (2*i-1)*h) for i in range(1, n+1))
		T2n = Tn/2. + h*F
		if (np.fabs(T2n - Tn) < (3*eps)):
			res = T2n
			break
		n *= 2
		h /= 2.
		Tn = T2n
	return res

# -------------------------------------------------------
# 函数名:half_integrate_parabola(f, a, b, eps)
# 功能:抛物线公式的逐次分半法求积分
# 说明:
# Parameters:
#	(a,b)积分区间。
#	f:function object 被积函数
#	eps最大允许误差
# Returns:
#	res:一个浮点数,Y在整个区间(a,b)上的积分
# -------------------------------------------------------
def half_integrate_parabola(f, a, b, eps):
	import numpy as np
	F1 = f(a) + f(b)
	F2 = f((b + a) / 2.)
	Sn = (b - a)/6. * (F1 + 4*F2)
	n = 2
	h = (b - a)/4.
	while(True):
		F3 = sum(f(a + (2*i-1)*h) for i in range(1, n+1))
		S2n = h/3.*(F1 + 2*F2 + 4*F3)
		if (np.fabs(S2n - Sn) < (3*eps)):
			res = S2n
			return n
			break
		n *= 2
		h /= 2.
		F2 = F2 + F3
		Sn = S2n
	return res

# -------------------------------------------------------
# 函数名:romberg_integrate(f, a, b, eps)
# 功能:Romberg求积
# 说明:
# Parameters:
#	(a,b)积分区间。
#	f:function object 被积函数
#	eps最大允许误差
#	out 是否输出运算过程 1为输出 0为不输出
# Returns:
#	res:一个浮点数,Y在整个区间(a,b)上的积分
# Tmat[i][j] = T_{i}^{j}
# -------------------------------------------------------
def romberg_integrate(f, a, b, mmax, out):
	n = 1
	h = (b-a)/2.
	Tmat = [[]]
	Tmat[0].append(h * (f(a) + f(b)))
	k = 1
	while(True):
		Tmat.append([])
		F = sum(f(a + (2*i-1)*h) for i in range(1, n+1))
		Tmat[0].append(Tmat[0][k-1]/2. + h*F)
		for m in range(1, k+1):
			m4 = 4**m
			app = (m4 * Tmat[m-1][k-m+1] - Tmat[m-1][k-m]) / (m4 - 1)
			Tmat[m].append(app)
		if (out):
			print '----------',k,'----------'
			for i in range(0, k+1):
				print(Tmat[i])
			print 'err: ', np.fabs(Tmat[m][0] - Tmat[m-1][0]) 
			print '-----------------------'
			print ' '
		#if (np.fabs(Tmat[m][0] - Tmat[m-1][0]) < eps):
		if (m == mmax):
			res = Tmat[m][0]
			break
		n *= 2
		h /= 2.
		k += 1
	return res

if __name__ == '__main__':
	def testf1(x):
		import numpy as np
		if (x==0):
			return 1
		else:
			return np.sin(x)/x
	def testf2(x):
		import numpy as np
		return np.e**(-x**2)*np.sqrt(1.+x**2)
	ans1 = integrate.quad(testf2,0,1)[0]
	err1 = np.zeros((101))
	err2 = np.zeros((101))
	elog1 = np.zeros((101))
	elog2 = np.zeros((101))
	for i in range(1,101):
		err1[i] = abs(c_integrate(testf2,0,1,i) - ans1)
		err2[i] = abs(p_integrate(testf2,0,1,i) - ans1)
	elog1[1:101] = np.log10(err1[1:101])
	elog2[1:101] = np.log10(err2[1:101])
	n = np.linspace(1,100,100)

#	err3 = np.zeros((7))
#	elog3 = np.zeros((7))
#	for i in range(1,7):
#		err3[i] = abs(romberg_integrate(testf2,0,1,i,0) - ans1)
#	elog3[1:7] = np.log10(err3[1:7])
#	plt.plot(n[0:6],elog3[1:7],'b*')
#	plt.plot(n[0:6],elog3[1:7],'b')
#	plt.ylim([elog3[-1],0])
#	plt.title(u'romberg求积法误差随外推次数的变化')
#	plt.ylabel(u'误差(取10的对数)')
#	plt.xlabel(u'外推次数')
#	plt.savefig(path + u'fig5.jpg')
	epsl = np.linspace(-1,-6,6)

	N1 = np.zeros(7)
	N2 = np.zeros(7)
	for i in epsl:
		print i
		N1[-i] = half_integrate_trapezoid(testf2,0,1,10**i)
		N2[-i] = half_integrate_parabola(testf2,0,1,10**i)
	plt.plot(-epsl,N1[1:7],'r*',label = u'梯形公式]')
	plt.plot(-epsl,N2[1:7],'b*',label = u'抛物线公式')
	plt.plot(-epsl,N1[1:7],'r')
	plt.plot(-epsl,N2[1:7],'b')
	plt.legend(loc ='upper right')
	plt.title(u'逐次分半法区间等分数随误差的变化')
	plt.ylabel(u'区间等分数')
	plt.xlabel(u'误差(取10的负对数)')
	plt.savefig(path + u'fig4.jpg')
	#plt.plot(n,elog1[1:101],'r',label = u'复化梯形公式')
	#plt.plot(n[1:100:2],elog2[1:51],'b',label = u'复化抛物线公式')
	#plt.legend(loc ='upper right')
	#plt.ylim([min(elog1[-1],elog2[-1]),0])
	#plt.title(u'复化公式误差随n的变化')
	#plt.ylabel(u'误差(取10的对数)')
	#plt.xlabel(u'区间等分数')
	#plt.savefig(path + u'fig3.jpg')