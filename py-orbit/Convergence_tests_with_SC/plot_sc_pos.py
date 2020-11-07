import matplotlib.pylab as plt
plt.style.use("IceCube")


f = open("./intensity_11/n_parts/SC_KICKSslice_01.txt", 'r')
pos_list = []


while 1:
    temp = f.readline().split()

    if len(temp) == 0:
        break


    pos_list.append(float(temp[0]))


pos_list_2 = []
pos_list_num = []

i = 1
pos_list_num.append(1)
pos_list_2.append(pos_list[0])
while i < len(pos_list):
    if pos_list[i] == pos_list_2[-1]:
        pos_list_num[-1] += 1
        i+=1
    else:
        pos_list_2.append(pos_list[i])
        pos_list_num.append(1)
        i+=1


f = plt.figure()

for i in range(len(pos_list_2))[:-1]:
    plt.plot([pos_list_2[i], pos_list_2[i]], [0, pos_list_num[i]])

plt.xlabel("position on lattice")
plt.show()

print(len(pos_list_2))
