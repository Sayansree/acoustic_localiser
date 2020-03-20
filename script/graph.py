import matplotlib.pyplot as plt
import numpy as np




x1 = []
y1 = []
z1=[]
for i in open("tdata1.txt",'r'):
    x1.append(float(i[:-1].split()[0]))
    y1.append(float(i[:-1].split()[1]))
plt.figure("test")
plt.subplot(421)
plt.plot(y1,x1)

x1 = []
y1 = []
z1=[]
for i in open("tdata2.txt",'r'):
    x1.append(float(i[:-1].split()[0]))
    y1.append(float(i[:-1].split()[1]))

plt.subplot(423)
plt.plot(y1,x1)

x1 = []
y1 = []
for i in open("fdata1.txt",'r'):
    x1.append(float(i[:-1].split()[0]))
    y1.append(float(i[:-1].split()[1]))
plt.subplot(422)
plt.plot(y1,x1)
x1 = []
y1 = []
for i in open("fdata2.txt",'r'):
    x1.append(float(i[:-1].split()[0]))
    y1.append(float(i[:-1].split()[1]))
plt.subplot(424)
plt.plot(y1,x1)
x1 = []
y1 = []
for i in open("fcorr.txt",'r'):
    x1.append(float(i[:-1].split()[0]))
    y1.append(float(i[:-1].split()[1]))
plt.subplot(426)
plt.plot(y1,x1)
x1 = []
y1 = []
for i in open("tcorr2.txt",'r'):
    x1.append(float(i[:-1].split()[0]))
    y1.append(float(i[:-1].split()[1]))
plt.subplot(425)
plt.plot(y1,x1)
x1 = []
y1 = []

plt.show()
