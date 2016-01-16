# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

x1 = np.array([1, 2, 3, 4, 5, 6, 7, 8])
x2 = np.array([2, 3, 4, 5, 7, 8, 10, 11, 14, 15, 16, 18, 19])
y1 = np.array([15.3, 20.5, 27.4, 36.6, 49.1, 65.6, 87.8, 117.6])
y2 = np.array([106.42, 108.20, 109.58, 109.50, 110.00, 109.93, 110.49, 110.59, 110.60, 110.90, 110.76, 111.00, 111.20])
x = x1
y = y1

print('经验公式:')
print('1、y = a+b*x')
print('2、y = a+b*x^2')
print('3、y = a*e^(bx)')
print('4、y = a*x^b')
print('5、y = a+b/x')

def fit(x, y, choose):
	if (choose == 1):
		u = y
		v = x
	elif (choose == 2):
		u = y
		v = x*x
	elif (choose == 3):
		u = np.log(y)
		v = x
	elif (choose == 4):
		u = np.log(y)
		v = np.log(x)
	elif (choose == 5):
		u = y
		v = 1.0/x

	#u = A + Bv
	n = v.shape[0]
	v1 = np.ones((n))
	v2 = v
	V = np.column_stack((v1, v2))
	VT = np.transpose(V)
	#VT*V*alpha = VT*u
	alpha = np.linalg.solve(np.dot(VT, V), np.dot(VT, u))
	A = alpha[0]
	B = alpha[1]

	u = A + B*v
	testX = np.linspace(x.min(), x.max(), 100)
	if (choose == 1):
		yfit = u
		a = A
		b = B
	elif (choose == 2):
		yfit = u
		a = A
		b = B
	elif (choose == 3):
		yfit = np.e**u
		a = np.e**A
		b = B
	elif (choose == 4):
		yfit = np.e**u
		a = np.e**A
		b = B
	elif (choose == 5):
		yfit = u
		a = A
		b = B
	Q = sum((yfit-y)**2)
	return (yfit, Q)

#----------main------------
plt.plot(x, y,    'o', color = 'r', label = u'data')
minChoose = 1
minQ = 9999999
for choose in range(1, 6):
	label = ['', 
			 r'$y^* = a + bx$  ', 
			 r'$y^* = a + bx^2$', 
			 r'$y^* = ae^{bx}$    ', 
			 r'$y^* = ax^b$     ', 
			 r'$y^* = a + \frac{b}{x}$    ']
	color = ['', 
			 'blue', 
			 'green', 
			 'grey', 
			 'yellow', 
			 'pink']
	ystar = fit(x, y, choose)
	if (ystar[1]<minQ):
		minQ = ystar[1]
		minChoose = choose
	plt.plot(x, ystar[0], '-', lw = 1.5, color = color[choose], label = label[choose]+' Q = %.2f' %ystar[1])
plt.plot(x, y,    'o', color = 'r')
plt.legend(loc = 'upper left', fontsize = '10')
plt.show()

label = ['', 
		 r'$y^* = a + bx$  ', 
		 r'$y^* = a + bx^2$', 
		 r'$y^* = ae^{bx}$    ', 
		 r'$y^* = ax^b$     ', 
		 r'$y^* = a + \frac{b}{x}$    ']
color = ['', 
		 'blue', 
		 'green', 
		 'grey', 
		 'yellow', 
		 'pink']
choose = minChoose
ystar = fit(x, y, choose)
plt.plot(x, y,    'o', color = 'r', label = u'data')
plt.plot(x, ystar[0], '-', lw = 1.5, color = color[choose], label = label[choose]+' Q = %.2f' %ystar[1])
plt.legend(loc = 'upper left', fontsize = '8')
plt.show()