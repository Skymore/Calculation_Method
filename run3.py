# -*- coding: utf-8 -*-

#------------------------------------------
#author	:陶睿
#email	:122345615@qq.com
#date	:2015-10-20
#第一个坐标为经度（或序号）x（1-159） 第二个坐标为纬度（或序号）y（1-185）
#------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

#xMax是经度序号最大值，yMAX是纬度序号最大值，tMAX是时间序号最大值
xMax = 159
yMax = 185
zMax = 31 * 24

path = ['/Users/sky/Desktop/SRDP/data/','/Users/sky/Desktop/SRDP/data2/']
figPath = '/Users/sky/Desktop/SRDP/fig/'
tidal = ['K1  ', 'M2  ']
years = [2008, 1985]

runNum = 1
tidalNum = 1

for runNum in range(0, 2): # 0, 1
	pha = np.loadtxt(path[runNum] + tidal[tidalNum] + '_pha.dat')
	amp = np.loadtxt(path[runNum] + tidal[tidalNum] + '_amp.dat')
	pos = np.loadtxt(path[runNum] + 'waterPoint.dat', skiprows = 2)

	fInput2 = open(path[runNum] + 'POM_GRID')
	lon = np.zeros((xMax + 1))
	lat = np.zeros((yMax + 1))
	coast = np.zeros((xMax + 1, yMax + 1)) # 999代表陆地
	lines = fInput2.readlines()
	# 读取经度
	for x in range(1, xMax + 1): 	#x = 1...xMax
		l = []
		for num in lines[(x - 1) * yMax].split():
			l.append(float(num))
		lon[x] = l[7]
	# 读取纬度
	for y in range(1, yMax + 1): 	#y = 1...yMax
		l = []
		for num in lines[y - 1].split():
			l.append(float(num))
		lat[y] = l[8]
	# 读取水深
	for x in range(1, xMax + 1):
		for y in range(1, yMax + 1):
			l = []
			for num in lines[(x - 1)*yMax + y - 1].split():
				l.append(float(num))
			if l[5] == 999:
				coast[x][y] = l[5]
	fInput2.close()

	print "第%d次运行"%runNum + tidal[tidalNum][0: 2] + '分潮:'
	print 'amp.min =', amp.min()
	print 'amp.max =', amp.max()
	print 'pha.min =', pha.min()
	print 'pha.max =', pha.max()
	eps = 0.01
	wetNum = pos.shape[0]
	amp2 = np.zeros((xMax + 1, yMax + 1))
	pha2 = np.zeros((xMax + 1, yMax + 1))
	for i in range(0, wetNum):
		x = pos[i][1]
		y = pos[i][2]
		amp2[x][y] = amp[i]
		pha2[x][y] = pha[i]

	#------------------------------------------
	#绘图
	#------------------------------------------

	fig = plt.figure(10*runNum + tidalNum)
	plt.title(u'%d年'%years[runNum] + tidal[tidalNum][0: 2] + u'分潮')
	ax = fig.add_subplot(111)
	Y, X = np.meshgrid(lat[10:180+1], lon[10:147+1])

	# 绘制边界
	lels = [-100, 998, 999]
	CS = ax.contourf(X, Y, coast[10:147+1, 10:180+1], lels, 
					 colors=('#AAAAAA', '#555555'))
	#CS = ax.contour(X, Y, coast[10:147+1, 10:180+1], lw = 3, colors='#555555')

	# 绘制等振幅填色图
	lels1 = np.linspace(1.2, 1.475, 11) # M2
	amp2 = np.ma.masked_outside(amp2, 1.2, 1.5) # M2
#	lels1 = np.linspace(0.32, 0.345, 1001) # K1	
#	amp2 = np.ma.masked_outside(amp2, 0.1, 0.5) # K1
	plt.contourf(X, Y, amp2[10:147+1, 10:180+1], 
				 lels1, extend = 'both')
	plt.colorbar(ticks = [1.2,1.25,1.3,1.35,1.4,1.45]) # M2
	plt.contour(X, Y, amp2[10:147+1, 10:180+1], 
				 lels1, linewidths = 2, linestyles = '--', extend = 'both',colors = 'blue')
#	plt.colorbar(ticks = [0.32, 0.325, 0.33, 0.335, 0.34, 0.345]) # K1

	# 绘制等迟角线
	lels2 = np.linspace(133,139,7) # M2
	pha2 = np.ma.masked_outside(pha2, 128, 139) # M2
#	lels2 = np.linspace(350.5,352,4) # K1
#	pha2 = np.ma.masked_outside(pha2, 300, 356) # K1
	CS2 = ax.contour(X, Y, pha2[10:147+1, 10:180+1], 
		lels2, linewidths = 2, colors = 'black')
	plt.clabel(CS2, 
			   inline=1, 
			   fmt = '%d', # M2
#			   fmt = '%.1f', # K1
			   fontsize=14)

	# 绘制坐标轴
	ticks = np.linspace(120.1,120.4,4)
	ax.xaxis.set_ticks(ticks)
	ax.xaxis.set_ticklabels(ticks)
	fig.savefig(figPath + u'%d年'%years[runNum] + 
				tidal[tidalNum][0: 2] + u'分潮4.jpg')
	plt.show()