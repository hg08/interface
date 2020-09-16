#Filename: 0_recenter_traj.py 
# Purpose: 
# 1. To recenter the atoms along the normal direction 

# To run this function: python3 Filename

import periodictable
import math 
import numpy as np

#=========
#CONSTANTS
#=========
atob=1.88971616463 #Factor to convert from Ang to bohr

#===============
#Basic questions
#===============
name=input("What is the name of the position file? (input)\n")
file_i = open(name,"r")

line=input("What are the crystal parameters? (a b c in Ang)\n")
line_split=line.split(' ')
a=float(line_split[0])
b=float(line_split[1])
c=float(line_split[2])

whish_size=float(input("What is aproximatively the size of a grid division? (in Ang)\n"))

n_axis=int(input("What is the normal axis? (0->x, 1->y, 2->z)\n"))

nb_steps=int(input("What is the total number of steps?\n"))

name_split=name.split('.')

name3 = 'recenter_{}.xyz'.format(name_split[0]) 
file_o3 = open(name3, 'w') #the recenter trajectory file as output
#END OF QUESTIONS

#============================
#Reading the main information
nb_atoms=int(file_i.readline()) #Number of atoms (the first line is the number of atoms)

# Some parameters
nb_divx=round(a/whish_size)  # number of division along x
nb_divy=round(b/whish_size)  #                          y
nb_divz=round(c/whish_size)  #                          z
divx=a/nb_divx  # length of each grid along x
divy=b/nb_divy  #                           y
divz=c/nb_divz  #                           z

#========================================
#nb_atoms lines about the atomic position
#atomic_number   random_float    x   y  z
#========================================
#Skip the commentary line of the xyz file
file_i.readline()

#1)Read all the atomic positions (and symbol)
#and calculate the position of the mass center
lsymb = [None]*nb_atoms # use a list of length 733 to store symb of atoms
lname = [None]*nb_atoms 
lx = [None]*nb_atoms
ly = [None]*nb_atoms
lz = [None]*nb_atoms
sum_mass=0
center_pos=0

# Go to start of file
file_i.seek(0)

for s in range(nb_steps):
    
    #skip the first line
    file_i.readline()
    # read the second line
    comment = file_i.readline()
    
    # write into file_o3
    file_o3.write("{0:5d}\n".format(nb_atoms))
    file_o3.write("{}".format(comment))
    
    # Read the information of the atoms    
    for i in range(nb_atoms):
        symb,x,y,z=file_i.readline().split()

        #Record the values inside a list with good format
        lsymb[i]=periodictable.elements.symbol(symb) #lsymb is an object
        lname[i]=symb
        lx[i]=float(x)
        ly[i]=float(y)
        lz[i]=float(z)
    
        #Calculate the mass center along the normal axis
        if n_axis==0:
            center_pos += lx[i]*lsymb[i].mass
        elif n_axis==1:
            center_pos += ly[i]*lsymb[i].mass
        else: 
            center_pos += lz[i]*lsymb[i].mass
        sum_mass +=lsymb[i].mass 

    #For each step s, we calculate the center_pos
    center_pos = center_pos / sum_mass

    #2)Recenter (along the normal axis only) the atoms according to the mass center
    if n_axis==0:
        for i in range(nb_atoms):
            lx[i]=lx[i]-center_pos +a/2
    elif n_axis ==1:
        for i in range(nb_atoms):
            ly[i]=ly[i]-center_pos +b/2
    else:
        for i in range(nb_atoms):
            lz[i]=lz[i]-center_pos +c/2 

    #3) Recenter all the atoms inside the cell (along all the directions)
    # and write the data into the cube file
    for i in range(nb_atoms):
        lx[i]=lx[i]%a; ly[i]=ly[i]%b; lz[i]=lz[i]%c
        #file_o3.write(" {0:5s} {1:12.6f} {2:12.6f} {3:12.6f}\n".format(lname[i],lx[i]*atob,ly[i]*atob,lz[i]*atob))
        # Do not change the units.
        file_o3.write(" {0:5s} {1:12.6f} {2:12.6f} {3:12.6f}\n".format(lname[i],lx[i],ly[i],lz[i]))
