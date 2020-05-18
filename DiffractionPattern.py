import pprint
import csv
import numpy as np
import matplotlib.pyplot as plt
import statistics 
import math
from numpy import linalg as LA

cif = open("C:/Users/Charles Cardot/Dropbox/Charles/YBCO_200h_1/YBCO_200hr_1.cif", "r")

pi = math.pi

wavelength = 1.5406*10**(-10)

a1 = np.array([3.8124*10**(-10),0,0])
a2 = np.array([0,3.8807*10**(-10),0])
a3 = np.array([0,0,11.6303*10**(-10)])

denom = np.dot(a1, np.cross(a2,a3))

b1 = (2*pi*np.cross(a2,a3))/denom
b2 = (2*pi*np.cross(a3,a1))/denom
b3 = (2*pi*np.cross(a1,a2))/denom

#print(LA.norm(b1))
#print(LA.norm(b2))
#print(LA.norm(b3))




theta = []
G = []
Gmag = []
AllValues = []
h=0
while(h<8):
    k=0
    while(k<8):
        l=0
        while(l<8):
            if(h==0 and k==0 and l==0):
                null = 0
                #to catch 1/0 --> Infinity edge case
            else:
                d=(2*pi)/(math.sqrt((h*LA.norm(b1))**2 + (k*LA.norm(b2))**2+(l*LA.norm(b3))**2))
                num = wavelength/(2*d)
                if(num <= 1):
                    temp = math.asin(num)
                    #print(h,k,l)
                    #print(temp*180/pi)
                    if(2.5*pi/180 <= temp <= 60*pi/180):
                        temp = temp*180/pi 
                        theta.append(temp)
                        G.append([h,k,l])
                        Gmag.append(math.sqrt((h*LA.norm(b1))**2 + (k*LA.norm(b2))**2+(l*LA.norm(b3))**2))
                        AllValues.append([temp, [h,k,l], math.sqrt((h*LA.norm(b1))**2 + (k*LA.norm(b2))**2+(l*LA.norm(b3))**2)])
            l=l+1
        k=k+1
    h=h+1



atoms = []
arr = []
for x in cif:
    arr.append(x)

#Building array with all atoms in the unit cell as well as fractional coordinates, named "atoms"
x=78
while(x<(78+8)):
    line = arr[x].split()
    name = line[0]
    #print(line)
    line[4] = line[4].replace('(', '')
    line[4] = line[4].replace(')', '')
    line[5] = line[5].replace('(', '')
    line[5] = line[5].replace(')', '')
    line[6] = line[6].replace('(', '')
    line[6] = line[6].replace(')', '')
    cord = [line[4], line[5], line[6]]
    atoms.append([name, cord])
    x=x+1


symmetries = []
#Building array with all symmetries in an array called symmetries
x=52
while(x<(52+8)):
    line = arr[x].split()
    line[0] = line[0].replace("'", '')
    line[0] = line[0].split(",")
    symmetries.append(line[0])
    x=x+1


#atom is the atom I want to find all the symmetries of
#arr is the most updated version of the "atoms" array
#symmetries is the array of opperations I will perform to atom to find all the new "atoms"
def BuildStructure(atom, arr, sym):
    #print(atom)
    for i in sym:
        
        xpos = atom[1][0]
        ypos = atom[1][1]
        zpos = atom[1][2]

        if i[0] == 'x':
            i[0] == 1
        elif i[0] == '-x':
            i[0] == -1
            #print("reached -x")
            xpos = str(1-float(xpos))
        if i[1] == 'y':
            i[1] == 1
        elif i[1] == '-y':
            i[1] == -1
            #print("reached -y")
            ypos = str(1-float(ypos))
        if i[2] == 'z':
            i[2] == 1
        elif i[2] == '-z':
            i[2] == -1
            #print("reached -z")
            zpos = str(1-float(zpos))
        cord = [xpos, ypos, zpos]

        #print(cord)
        skip = False
        if(float(xpos) == 1 or float(ypos) == 1 or float(zpos) == 1):
            skip = True
        for j in arr:
            if(j[1] == cord):
                #Just skips these because the position is already occupied by another atom
                skip = True


        if(skip == False):
            arr.append([atom[0], cord])
        
    return arr
            
holder = atoms
atoms = BuildStructure(holder[0], atoms, symmetries)
atoms = BuildStructure(holder[1], atoms, symmetries)
atoms = BuildStructure(holder[2], atoms, symmetries)
atoms = BuildStructure(holder[3], atoms, symmetries)
atoms = BuildStructure(holder[4], atoms, symmetries)
atoms = BuildStructure(holder[5], atoms, symmetries)
atoms = BuildStructure(holder[6], atoms, symmetries)
atoms = BuildStructure(holder[7], atoms, symmetries)

#Gives atomic number and approximate radius of each atom
atomdata = [['Ba1', [56, 149*10**(-12)]], ['Y1', [39, 104*10**(-12)]], ['Cu1', [29, 87*10**(-12)]], ['Cu2', [29, 87*10**(-12)]], ['O1', [8, 126*10**(-12)]], ['O2', [8, 126*10**(-12)]], ['O3', [8, 126*10**(-12)]], ['O4', [8, 126*10**(-12)]]]

#takes in the atom name and the magnitude of G
def f(atomname, mag):
    for i in atomdata:
        if i[0] == atomname:
            x = float(mag)*i[1][1]
            num = 3*i[1][0]*(np.sin(x)-x*np.cos(x))/(x**(3))
            return num



#Finds the structure factor for every allowable hkl value
#g has hkl info
#mag has magnitude of G vector info
#arr is array of atoms
def AbsSquareStructFactor(g, mag, arr):
    F = 0
    f_temp = 0
    for atom in arr:
        #print(atom)
        #print(mag)
        f_temp = f(atom[0], mag)
        num = f_temp*math.e ** (2*pi * 1j)*(float(g[0])*float(atom[1][0])+float(g[1])*float(atom[1][1])+float(g[2])*float(atom[1][2]))
        F = F + num
    return abs(F)*abs(F)


x = []
y = []
for i in AllValues:
    #making the x values 2 theta values
    x.append(round(2*float(i[0]), 2))
    F = AbsSquareStructFactor(i[1], i[2], atoms)
    #making the y values the structure factor squared
    y.append(round(float(F)))


holder = []
temp = [x,y]
for i in range(len(temp[0])):
    holder.append([temp[0][i], temp[1][i]])
holder = np.array(holder)
holder = holder[np.argsort(holder[:,0])]

x = []
y = []

for i in range(len(holder)):
    x.append(holder[i][0])
    y.append(holder[i][1])

#adds in zero values for y at every 0.02 x whenever there is no peak
z = 5
num = x[0]
i = 0
while(z <= 120):
    while(z < num):
        x.insert(i, z)
        y.insert(i, 0)
        z = z + 0.02
        i = i + 1
        #print(z)
    if(z > 118):
        z = 121
    else:
        z = z + 0.02
        num = x[i+1]
        i = i + 1

ymax = 0
for i in y:
    if i > ymax:
        ymax = i

for i in range(len(y)):
    y[i] = y[i]/(ymax)


title = "X-Ray Diffraction Theoretical Plot"
name = "Xray_YBCO_Theory"

fig = plt.figure(figsize = (6,4))
plt.scatter(x, y, label = "Data")
plt.plot(x, y)
plt.xlabel("2 Theta")
plt.ylabel("Intensity")
plt.legend(loc='best')
plt.title(title)
fig.savefig(name)

plt.show()
