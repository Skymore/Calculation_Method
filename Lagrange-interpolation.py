# -*- coding: utf-8 -*-

# Python Lagrange-interpolation 拉格朗日n次插值
#Author	: Skymore 122345615@qq.com
#Date	: 2015-10-10

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

#对数据（Xi, Yi）进行插值，次数为X.size
def P(x, X, Y):
	def L(i, x):
		Li = 1
		for j in range(0, X.size):
			if(j != i):
				Li *= (x - X[j]) / (X[i] - X[j])
		return Li
	ans = 0
	for i in range(0, X.size):
		ans += Y[i]*L(i, x)
	return ans

y = lambda x: 1 / (1 + 25 * x ** 2)
#y = lambda x: 3 * x**3 + 2.3 * x **2 + 1.5

if __name__ == '__main__':
	#X Y为生成插值函数的数据，testX testY为测试插值函数图像的数据
	#n为插值次数,X[i](i = 0...n)    X[0] = a, X[n] = b 等间距把[a,b]分为n份

	a = -1
	b = 1

	nBest = 2
	errBest = 9999
	for n in range(2, 31, 2):
		X = np.linspace(a, b, n + 1)
		Y = y(X)

		#Pn为进行n(n = X.size)次插值后得到的函数。
		Pn = lambda x: P(x, X, Y)

		integrand = lambda x: abs(Pn(x) - y(x))

		err ,error = integrate.quad(integrand, a, b, limit = 101)

		print "n = ", n, "err = ", err
		if(err < errBest):
			nBest = n
			errBest = err
	print "nBest = ",nBest

	n = 10
	X = np.linspace(a, b, n + 1)
	Y = y(X)
	Pn = lambda x: P(x, X, Y)
	testX = np.linspace(a,b,200001)
	testY = y(testX)
	testF = Pn(testX)
	#for i in range(0, testX.size):
	#	if (abs(testY[i] - testF[i])> 0.01):
	#		print i
	#		print "ERROR"

	#up把插值所得的函数图像向上提升up个数值，方便比较
	up = 0
	figure = plt.figure()
	plt.plot(X, Y, 'ro')
	plt.plot(testX, testY, color = "green")
	plt.plot(testX, testF + up, color = "blue")
	plt.show()