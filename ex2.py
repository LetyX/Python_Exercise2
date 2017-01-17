import os
import re
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math


def get_data(name):  # returns r, theta and the corresponding energies as a list of tuples
    list = []
    dir = os.listdir(name)
    for file in dir:
        filename = str(file)
        m = re.match("\w\d\w\.r(\d\.\d+)theta(\d+)\.\d\.out$", filename)
        r = round(float(m.group(1)), 2)
        theta = int(m.group(2))
        f = open(os.path.join(name, file), "r")
        for line in f:
            if str("SCF Done:") in line:
                l = line.split()
                e = round(float(l[4]), 10)
        f.close()
        list.append((r, theta, e))

    data = np.asarray(list)
    return data


def zmatrix(data, r, theta):  # arranges the z values in a 2D array suitable for a 3D plot
    dim1 = r.size
    dim2 = theta.size
    e = np.zeros((dim1, dim2))

    for k in range(0, int(data.size / 3)):
        rindex = float((data[k, 0] - min(r)) / 0.05)  # index of corresponding r
        i = int(round(rindex, 0))  # convert index to an integer
        j = int(data[k, 1] - min(theta))  # index of corresponding theta
        e[i, j] = data[k, 2]

    return e


def plot_energy(x, y, z, molecule):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x, y = np.meshgrid(x, y)
    ax.plot_surface(x, y, np.transpose(z), rstride=2, cstride=2, linewidth=0.5, cmap='Spectral',
                    antialiased=True)

    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(y), np.max(y))
    ax.set_zlim(np.min(z), np.max(z))

    ax.set_xlabel('r /${\AA}$')
    ax.set_ylabel('Theta /deg')
    ax.set_zlabel('Energy /Hartree')

    ax.zaxis.set_major_locator(LinearLocator(6))  # 6 evenly-spaced tick marks on the z axis
    ax.yaxis.set_major_locator(LinearLocator(6))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_title("Potential Energy Surface for " + molecule)  # title

    ax.view_init(elev=30, azim=285)  # elevation & angle
    ax.dist = 9.2  # distance from the plot

    plt.savefig("myplot.png")


def freq(data, r, theta, e):
    # Determine equilibrium energy and geometry
    eeq = np.min(e)
    for k in range(0, int(data.size / 3)):
        if data[k, 2] == eeq:
            req = data[k, 0]
            teq = data[k, 1]

    rindex = float((req - np.min(r))/0.05)
    a = int(round(rindex, 0))  # index in r corresponding to req
    b = int(teq - np.min(theta))  # index in theta corresponding to teq

    # Calculate force constant k_r
    x1 = np.zeros(6)
    e1 = np.zeros(6)

    for i in range(0, 6):  # Use 6 points around the lowest energy point
        if i < 3:
            x1[i] = 0.5 * (r[a - 3 + i] - req)**2
            e1[i] = e[a - 3 + i, b]
        else:
            x1[i] = 0.5 * (r[a + i - 2] - req)**2
            e1[i] = e[a + i - 2, b]

    p = np.polyfit(x1, e1, 1)
    kr = p[0]

    # Calculate force constant k_theta
    x2 = np.zeros(8)
    e2 = np.zeros(8)

    for i in range(0, 8):  # Use 8 points around the lowest energy point
        if i < 4:
            x2[i] = 0.5 * ((theta[b - 4 + i] - teq) * math.pi/180)**2  # convert theta from deg to rad
            e2[i] = e[a, b - 4 + i]
        else:
            x2[i] = 0.5 * ((theta[b + i - 3] - teq) * math.pi/180)**2
            e2[i] = e[a, b + i - 3]

    p = np.polyfit(x2, e2, 1)
    ktheta = p[0]

    # Calculate frequencies
    nu1 = math.sqrt(kr/2)/(2 * math.pi)
    nu2 = math.sqrt(ktheta / (0.5 * req**2))/(2 * math.pi)

    return nu1, nu2


directory = str(input("Enter a directory pathname: "))
name = directory[0:3]  # get the molecular formula

alldata = get_data(directory)
r = np.arange(min(alldata[:, 0]), max(alldata[:, 0]) + 0.05, 0.05)
theta = np.arange(min(alldata[:, 1]), max(alldata[:, 1]) + 1, 1)
energ = zmatrix(alldata, r, theta)

plot_energy(r, theta, energ, name)

frequencies = freq(alldata, r, theta, energ)
# Data in frequencies has units of (1/angstrom)*sqrt(Hartree/mu)

print("The vibrational frequencies in cm-1 are:")
print("Symmetric stretch", int(frequencies[0] * 17091.7), sep='   ')
print("Bend", int(frequencies[1] * 17091.7), sep='   ')
