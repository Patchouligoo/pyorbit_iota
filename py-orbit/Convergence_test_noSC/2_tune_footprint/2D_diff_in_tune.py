import matplotlib.pylab as plt
plt.style.use("IceCube")
from math import *
import numpy as np




f = open("tune_diff.txt", 'r')
diff = []
while 1:
    temp = f.readline().split()

    if len(temp) < 1:
        break

    diff.append(float(temp[0]))



diff_2D = [[] for i in range(100)]
curr = 0
for i in range(100):
    for j in range(100):
        diff_2D[99 - j].append(diff[curr])
        curr += 1


f = plt.figure()
plt.imshow(diff_2D, extent=(0, 5, 0, 5))
colorbar1 = plt.colorbar()
colorbar1.set_label('normal scale')
plt.xlabel('x in sigma_x')
plt.ylabel('y in sigma_y')


for i in range(len(diff_2D)):
    for j in range(len(diff_2D[i])):
        if diff_2D[i][j] == 0:
            diff_2D[i][j] = np.nan
        else:
            diff_2D[i][j] = log10(diff_2D[i][j])

f2 = plt.figure()
plt.imshow(diff_2D, extent=(0, 5, 0, 5))
colorbar2 = plt.colorbar()
colorbar2.set_label('log 10 sclae')
plt.xlabel('x in sigma_x')
plt.ylabel('y in sigma_y')
plt.show()
