# -*- coding: utf-8 -*-

#------------------------------------------
#author	:陶睿
#email	:122345615@qq.com
#date	:2015-10-20
#第一个坐标为经度x（1-xMax） 第二个坐标为纬度y（1-yMax）
#------------------------------------------

import numpy as np

#xMax是经度序号最大值，yMAX是纬度序号最大值，tMAX是时间序号最大值
xMax = 159
yMax = 185
tMax = 31 * 24

path = '/Users/sky/Desktop/SRDP/data/'

#------------------------------------------
#读入潮位数据，h为潮位数据
#------------------------------------------
month = 5 
dayMin = 1
dayMax = 31
hourMin = 0
hourMax = 23
timeCount = 0

h = np.zeros((xMax + 1, yMax + 1, tMax + 1))
h = h-1000
for day in range(dayMin, dayMax + 1):			#day = 1...31
	for hour in range(hourMin, hourMax + 1):	#hour = 0...23
		timeCount += 1							#时间序号
		#print("%02d%02d%02d"%(monthYY, day, hour))
		fInput1 = open("/Users/sky/Desktop/SRDP/POM_output/hydro/hydro2011%02d%02d%02d.dat"%(month,day,hour))
		for line in fInput1:
		    l = []
		    for num in line.split():
		        l.append(float(num))
		    x = int(l[0])
		    y = int(l[1])
		    h[x,y,timeCount] = l[8]
		fInput1.close()

#------------------------------------------
#读入经纬度数据
#lon[1:xMax + 1]共xMax个经度对应x:1-xMax
#lat[1:yMax + 1]共yMax个纬度对应y:1-yMax
#------------------------------------------
lon = np.zeros((xMax + 1))
lat = np.zeros((yMax + 1))
fInput2 = open("/Users/sky/Desktop/SRDP/POM_output/POM_GRID")
lines = fInput2.readlines()
for x in range(1, xMax + 1): 	#x = 1...xMax
	l = []
	for num in lines[(x - 1) * yMax].split():
		l.append(float(num))
	lon[x] = l[7]
for y in range(1,yMax + 1): 	#y = 1...yMax
	l = []
	for num in lines[y - 1].split():
		l.append(float(num))
	lat[y] = l[8]	
fInput2.close()

#------------------------------------------
#求至少tMax-5个时间节点都为水点的点
#标记水点的wetFlag[x,y]= 1
#wetNum 水点的总个数
#------------------------------------------
wetFlag = np.zeros((xMax + 1, yMax + 1))
wetNum = 0
for x in range(1, xMax + 1):
	for y in range(1, yMax + 1):
		noZero = np.nonzero(h[x,y]+1000)
		length = len(h[x,y][noZero])
		if (length >= tMax - 5):
			wetFlag[x,y] = 1
			wetNum += 1
print("wetNum = %d"%wetNum)

#------------------------------------------
#输出水点的经纬度及其序号到文件(供matlab使用)
#------------------------------------------
f = open(path + "waterPoint.dat", "w")
f.write("共有 %d 个水点\n"%wetNum)
f.write("数据格式:wetNum x y lon lat\n")
i = 0
for x in range(1, xMax + 1):
	for y in range(1, yMax + 1):
		if (wetFlag[x,y] == 1):
			i += 1
			f.write(str(i))
			f.write(' ')
			f.write(str(x))
			f.write(' ')
			f.write(str(y))
			f.write(' ')
			f.write(str(lon[x]))
			f.write(' ')
			f.write(str(lat[y]))
			f.write(' ')
			f.write('\n')

#------------------------------------------
#输出水点的潮位时间序列(供matlab使用)
#------------------------------------------
f = open(path + "waterLevel.dat", "w")
f.write("共有 %d 个水点\n"%wetNum)
f.write("数据格式:wetNum time h\n")
i = 0
for x in range(1, xMax + 1):
	for y in range(1, yMax + 1):
		if (wetFlag[x, y] == 1):
			i += 1
			for timeCount in range(1, tMax + 1):
				f.write(str(i))
				f.write(' ')
				f.write(str(timeCount))
				f.write(' ')
				f.write(str(h[x,y,timeCount]))
				f.write('\n')
f.close()
#------------------------------------------
#画潮汐变化图
#------------------------------------------
levels = np.linspace(-4.36,3.85,200)
H = np.ma.masked_equal(h,-1000)
H = H[1:,1:,:]

fig = plt.figure()   
ax = plt.axes(xlim=(0, 160), ylim=(0, 186))   
line, = ax.contourf(H[:,:,1],levels)  

def init():   
    line.set_data([],[],[])   
    return line,

def animate(i):     
    line.set_data(H[:,:,i])   
    return line,
anim = animation.funcanimation(fig, animate, init_func=init,   
                               frames=200, interval=20, blit=true)