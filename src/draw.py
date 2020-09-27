import numpy as np
import matplotlib.pyplot as plt
import sys
from sys import argv
# Enter "python3 draw.py 1", then sys.argv[0] would be "draw.py"

dir = ""
if sys.argv[1]=="1":
    dir += "../log/Brio_Wu"
elif sys.argv[1]=="2":
    dir += "../log/RJ2a"
elif sys.argv[1]=="3":
    dir += "../log/Falle"
elif sys.argv[1]=="4":
    dir += "../log/RJ4d"
elif sys.argv[1]=="5":
    dir += "../log/Komissarov"

data_file = dir + "/d.txt"
data_save = dir + "/d.jpeg"

arr_y = []
count = 0
with open(data_file) as my_file:
    for line in my_file:
        arr_y.append(line)
        count = count + 1
arr_x = np.linspace(0,1,count)
plt.figure(1)
plt.plot(arr_x, arr_y)
plt.xlabel(r"$x$")
plt.ylabel(r"$\rho$")
plt.axis([0, 1, None, None])
plt.rcParams.update({"font.size": 16})
plt.savefig(data_save)
#plt.show()

data_file = dir + "/ux.txt"
data_save = dir + "/ux.jpeg"
arr_y = []
count = 0
with open(data_file) as my_file:
    for line in my_file:
        arr_y.append(line)
        count = count + 1
arr_x = np.linspace(0,1,count)
plt.figure(2)
plt.plot(arr_x, arr_y)
plt.xlabel(r"$x$")
plt.ylabel(r"$u_x$")
plt.axis([0, 1, None, None])
plt.rcParams.update({"font.size": 16})
plt.savefig(data_save)
#plt.show()

data_file = dir + "/uy.txt"
data_save = dir + "/uy.jpeg"
arr_y = []
count = 0
with open(data_file) as my_file:
    for line in my_file:
        arr_y.append(line)
        count = count + 1
arr_x = np.linspace(0,1,count)
plt.figure(3)
plt.plot(arr_x, arr_y)
plt.xlabel(r"$x$")
plt.ylabel(r"$u_y$")
plt.axis([0, 1, None, None])
plt.rcParams.update({"font.size": 16})
plt.savefig(data_save)
#plt.show()

data_file = dir + "/uz.txt"
data_save = dir + "/uz.jpeg"
arr_y = []
count = 0
with open(data_file) as my_file:
    for line in my_file:
        arr_y.append(line)
        count = count + 1
arr_x = np.linspace(0,1,count)
plt.figure(4)
plt.plot(arr_x, arr_y)
plt.xlabel(r"$x$")
plt.ylabel(r"$u_z$")
plt.axis([0, 1, None, None])
plt.rcParams.update({"font.size": 16})
plt.savefig(data_save)
#plt.show()

data_file = dir + "/P.txt"
data_save = dir + "/P.jpeg"
arr_y = []
count = 0
with open(data_file) as my_file:
    for line in my_file:
        arr_y.append(line)
        count = count + 1
arr_x = np.linspace(0,1,count)
plt.figure(5)
plt.plot(arr_x, arr_y)
plt.xlabel(r"$x$")
plt.ylabel(r"$P$")
plt.axis([0, 1, None, None])
plt.rcParams.update({"font.size": 16})
plt.savefig(data_save)
#plt.show()

data_file = dir + "/Bx.txt"
data_save = dir + "/Bx.jpeg"
arr_y = []
count = 0
with open(data_file) as my_file:
    for line in my_file:
        arr_y.append(line)
        count = count + 1
arr_x = np.linspace(0,1,count)
plt.figure(6)
plt.plot(arr_x, arr_y)
plt.xlabel(r"$x$")
plt.ylabel(r"$B_x$")
plt.axis([0, 1, None, None])
plt.rcParams.update({"font.size": 16})
plt.savefig(data_save)
#plt.show()

data_file = dir + "/By.txt"
data_save = dir + "/By.jpeg"
arr_y = []
count = 0
with open(data_file) as my_file:
    for line in my_file:
        arr_y.append(line)
        count = count + 1
arr_x = np.linspace(0,1,count)
plt.figure(7)
plt.plot(arr_x, arr_y)
plt.xlabel(r"$x$")
plt.ylabel(r"$B_y$")
plt.axis([0, 1, None, None])
plt.rcParams.update({"font.size": 16})
plt.savefig(data_save)
#plt.show()

data_file = dir + "/Bz.txt"
data_save = dir + "/Bz.jpeg"
arr_y = []
count = 0
with open(data_file) as my_file:
    for line in my_file:
        arr_y.append(line)
        count = count + 1
arr_x = np.linspace(0,1,count)
plt.figure(8)
plt.plot(arr_x, arr_y)
plt.xlabel(r"$x$")
plt.ylabel(r"$B_z$")
plt.axis([0, 1, None, None])
plt.rcParams.update({"font.size": 16})
plt.savefig(data_save)
#plt.show()