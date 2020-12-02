#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import sys
from math import acos, degrees, cos, pi, exp
from IPython.core.interactiveshell import InteractiveShell
from time import time
import pandas as pd 

InteractiveShell.ast_node_interactivity = "all"
np.set_printoptions(threshold=sys.maxsize)


# In[2]:


#get the position from 1 to 10000
natoms=20
#index start from 1
st=1
ed=10000
pos=np.zeros((ed-st+1,natoms,3))
f=open("XDATCAR")
tmp = f.readlines()
f.close()
lattice = np.array([9.047747, 9.047747, 12.940000]).reshape(1,1,3)
for step in range(st, ed+1):
    for atom in range(natoms):
        ln = 8 + (step-1)*(natoms+1) + atom
        pos[step-st, atom] = np.array([float(x) for x in tmp[ln].strip().split()])

#pos *= lattice
#np.save("pos.npy", pos)


# In[4]:


import matplotlib.pyplot as plt


# In[5]:


x = np.arange(-1.3*pi, 2*pi, 0.1)
y = (1+np.cos(x-2))**1


# In[6]:


plt.plot(x,y)# 


# In[7]:


x = np.arange(-4,8, 0.1)
y = np.exp(-0.15*(x-2)**2)


# In[8]:


plt.plot(x,y)# 


# In[9]:


#Minimum image convention

def mic_shortest_dis_xyz(mol_pos):
    x_pos = mol_pos[:,0]
    y_pos = mol_pos[:,1]
    z_pos = mol_pos[:,2]
    
    x_dis = np.zeros((20,20))
    y_dis = np.zeros((20,20))
    z_dis = np.zeros((20,20))
    
    #calculate shorest distance along x axis 
    for i in range(0,20):
        short_dis = []
        for x in x_pos:
            if abs(x_pos[i] - x) > 0.5:
                x_del = min(abs(x_pos[i] - x), abs(x_pos[i] - x+1), abs(x_pos[i] - x-1))
                short_dis.append(x_del)
            else:
                x_del = abs(x_pos[i] - x)
                short_dis.append(x_del)

        x_dis[i] = short_dis

    #calculate shorest distance along y axis of all pairs of atoms
    #To make it clear, 3 dimentions are separate 
    for i in range(0,20):
        short_dis = []
        for y in y_pos:
            if abs(y_pos[i] - y) > 0.5:
                y_del = min(abs(y_pos[i] - y), abs(y_pos[i] - y+1), abs(y_pos[i] - y-1))
                short_dis.append(y_del)
            else:
                y_del = abs(y_pos[i] - y)
                short_dis.append(y_del)

        y_dis[i] = short_dis

    #calculate shorest distance along z axis
    for i in range(0,20):
        short_dis = []
        for z in z_pos:
            if abs(z_pos[i] - z) > 0.5:
                z_del = min(abs(z_pos[i] - z), abs(z_pos[i] - z+1), abs(z_pos[i] - z-1))
                short_dis.append(z_del)
            else:
                z_del = abs(z_pos[i] - z)
                short_dis.append(z_del)

        z_dis[i] = short_dis
        
    return(x_dis, y_dis, z_dis)



# In[11]:


def mic_shortest_dis(mic_shortest_dis_xyz):
    lattice_x = lattice[0][0][0]
    lattice_y = lattice[0][0][1]
    lattice_z = lattice[0][0][2]        
    
    coord = mic_shortest_dis_xyz
    
    dis = np.zeros((20,20))
    
    for m in range (0, 20):  
        mic_shortest_dis = []
        for i in range (0, 20):
            x = coord[0][m][i]*lattice_x
            y = coord[1][m][i]*lattice_y
            z = coord[2][m][i]*lattice_z
            mic_shortest_dis.append(np.sqrt(x*x + y*y + z*z))
            
        dis[m] = mic_shortest_dis
    
    return(dis)
    
    
    
    
    
    
    
    


# In[12]:


lattice_x = lattice[0][0][0]
lattice_y = lattice[0][0][1]
lattice_z = lattice[0][0][2]

np.sqrt((0.5*lattice_x)**2 + (0.5*lattice_y)**2 + (0.5*lattice_z)**2)


# In[13]:


def min_max_dis (dis_data):
    #min & max Cs dis
    dis_Cs_Cs = []
    dis_Cs_Pb = []
    dis_Cs_I = []
    for i in range (0,4):
        for m in range (i+1, 4):
            dis_Cs_Cs.append(dis_data[i][m])
        for n in range (4, 8):  
            dis_Cs_Pb.append(dis_data[i][n])
        for k in range (8, 20):
            dis_Cs_I.append(dis_data[i][k])


    min_dis_Cs_Cs = min(dis_Cs_Cs)
    max_dis_Cs_Cs = max(dis_Cs_Cs)
    min_dis_Cs_Pb = min(dis_Cs_Pb)
    max_dis_Cs_Pb = max(dis_Cs_Pb)
    min_dis_Cs_I = min(dis_Cs_I)   
    max_dis_Cs_I = max(dis_Cs_I)

    #min & max Pb dis
    dis_Pb_Cs = dis_Cs_Pb
    dis_Pb_Pb = []
    dis_Pb_I = []


    for i in range(4,8):
        for m in range (i+1, 8):
            dis_Pb_Pb.append(dis_data[i][m])
        for n in range (8,20):
            dis_Pb_I.append(dis_data[i][n])

    min_dis_Pb_Cs = min_dis_Cs_Pb
    max_dis_Pb_Cs = max_dis_Cs_Pb
    min_dis_Pb_Pb = min(dis_Pb_Pb )
    max_dis_Pb_Pb = max(dis_Pb_Pb )
    min_dis_Pb_I = min(dis_Pb_I)   
    max_dis_Pb_I = max(dis_Pb_I)

    #min & max I dis
    dis_I_Cs = dis_Cs_I
    dis_I_Pb = dis_Pb_I
    dis_I_I = []


    for i in range (8,20):
        for m in range (i+1, 20):
            dis_I_I.append(dis_data[i][m])

    min_dis_I_Cs = min_dis_Cs_I
    max_dis_I_Cs = max_dis_Cs_I
    min_dis_I_Pb = min_dis_Pb_I
    max_dis_I_Pb = max_dis_Pb_I
    min_dis_I_I = min(dis_I_I)   
    max_dis_I_I = max(dis_I_I)

    min_dis_data = np.zeros((3,3))
    max_dis_data = np.zeros((3,3))
    min_dis_data[0][0], min_dis_data[0][1], min_dis_data[0][2] = min_dis_Cs_Cs, min_dis_Cs_Pb, min_dis_Cs_I
    min_dis_data[1][0], min_dis_data[1][1], min_dis_data[1][2] = min_dis_Pb_Cs, min_dis_Pb_Pb, min_dis_Pb_I
    min_dis_data[2][0], min_dis_data[2][1], min_dis_data[2][2] = min_dis_I_Cs, min_dis_I_Pb, min_dis_I_I
    
    max_dis_data[0][0], max_dis_data[0][1], max_dis_data[0][2] = max_dis_Cs_Cs, max_dis_Cs_Pb, max_dis_Cs_I
    max_dis_data[1][0], max_dis_data[1][1], max_dis_data[1][2] = max_dis_Pb_Cs, max_dis_Pb_Pb, max_dis_Pb_I
    max_dis_data[2][0], max_dis_data[2][1], max_dis_data[2][2] = max_dis_I_Cs, max_dis_I_Pb, max_dis_I_I

    return (min_dis_data, max_dis_data)


# In[14]:


def dis_delta (min_max_dis):
    min_dis = np.array(min_max_dis[0])
    max_dis = np.array(min_max_dis[1])    
    dis_delta = max_dis - min_dis    
    return(dis_delta)
    


# In[15]:


def Rs_decision(pos):
    step1 = mic_shortest_dis_xyz(pos)
    step2 = mic_shortest_dis(step1)
    step3 = min_max_dis(step2)    
    return(step3)


# In[ ]:


dis_min_data = []
dis_max_data = []

for i in range(0,10000):
    dis_min_data.append(Rs_decision(pos[i])[0])
    dis_max_data.append(Rs_decision(pos[i])[1])
    
dis_min_data = np.array(dis_min_data)
dis_max_data = np.array(dis_max_data)

dis_minavg_data = np.average(dis_min_data, axis = 0)
dis_maxavg_data = np.average(dis_max_data, axis = 0)
    
dis_avg_delta = dis_maxavg_data - dis_minavg_data
    
    
    
    
    


# In[ ]:


dis_minavg_data


# In[ ]:


dis_maxavg_data


# In[ ]:


dis_avg_delta


# In[16]:


def angle_BAC (mic_shortest_dis):
    
    dis_data = mic_shortest_dis    
    angle_data = np.zeros((20,190))
       
    for n in range(0,20):
        angle_data_each_atom = []
        for i in range(0,20):
            for m in range(i+1,20):
                A = dis_data[n][i]
                B = dis_data[n][m]
                C = dis_data[i][m]

                if A==0. or B==0. or C==0.:
                    angle = 0
                else:
                    angle = acos((A * A + B * B - C * C)/(2.0 * A * B))         
                angle_data_each_atom.append(angle)
        angle_data[n] = angle_data_each_atom

    return (angle_data)
    
    
    
    
    
    
    


# In[17]:


def angle_preprocess(angle_data):
    ind_sec_atom = np.arange(19,0,-1,dtype = int)
    
    angle_fir_atom = []
    
    for n in range(0,20):
        angle_sec_atom = []
        i = 0
        m = 0
        while i < 190:
            angle_sec_atom.append(angle_data[n][i: i+ind_sec_atom[m]])
            i = i+ ind_sec_atom[m]
            m = m +1
        angle_fir_atom.append(angle_sec_atom)
    return (angle_fir_atom)
    
    
    
    


# In[18]:


def ang_Cs(angle_Cs_data):
    # angle Cs-Cs-Cs, Cs-Cs-Pb, Cs-Cs-I
    ang_CsCsCs = []
    ang_CsCsPb = []
    ang_CsCsI = []

    for i in range(0,4):
        for m in range (0, 3-i):
            if angle_Cs_data[i][m] == 0:
                pass
            else:
                ang_CsCsCs.append(angle_Cs_data[i][m])

        for n in range(3-i, 7-i):
            if angle_Cs_data[i][n] == 0:
                pass
            else:
                ang_CsCsPb.append(angle_Cs_data[i][n])

        for k in range(7-i, 19-i):
            if angle_Cs_data[i][k] == 0:
                pass
            else:
                ang_CsCsI.append(angle_Cs_data[i][k])

    #angle Pb-Cs-I, Pb-Cs-Pb

    ang_PbCsPb = []
    ang_PbCsI = []
    q = 0

    for i in range(4, 8):
        for m in range(0, 3-q):
            if angle_Cs_data[i][m] == 0:
                pass          
            else:
                ang_PbCsPb.append(angle_Cs_data[i][m])

        for n in range(3-q, 15-q):
            if angle_data[i][n] == 0:
                pass          
            else:
                ang_PbCsI.append(angle_Cs_data[i][n])        
        q += 1

    #angle I-Cs-I
    q = 0
    ang_ICSI = []
    for i in range (8, 19):
        for m in range(0, 11-q):
            if angle_Cs_data[i][m] == 0:
                pass
            else:
                ang_ICSI.append(angle_Cs_data[i][m])
        q +=1

    return(ang_CsCsCs, ang_CsCsPb, ang_CsCsI, ang_PbCsPb, ang_PbCsI, ang_ICSI)


# In[19]:


def ang_Pb(angle_Pb_data):
    # angle Cs-Pb-Cs, Cs-Pb-Pb, Cs-Pb-I
    ang_CsPbCs = []
    ang_CsPbPb = []
    ang_CsPbI = []

    for i in range(0,4):
        for m in range (0, 3-i):
            if angle_Pb_data[i][m] == 0:
                pass
            else:
                ang_CsPbCs.append(angle_Pb_data[i][m])

        for n in range(3-i, 7-i):
            if angle_Pb_data[i][n] == 0:
                pass
            else:
                ang_CsPbPb.append(angle_Pb_data[i][n])

        for k in range(7-i, 19-i):
            if angle_Pb_data[i][k] == 0:
                pass
            else:
                ang_CsPbI.append(angle_Pb_data[i][k])

    #angle Pb-Pb-Pb, Pb-Pb-I

    ang_PbPbPb = []
    ang_PbPbI = []
    q = 0

    for i in range(4, 8):
        for m in range(0, 3-q):
            if angle_Pb_data[i][m] == 0:
                pass          
            else:
                ang_PbPbPb.append(angle_Pb_data[i][m])

        for n in range(3-q, 15-q):
            if angle_Pb_data[i][n] == 0:
                pass          
            else:
                ang_PbPbI.append(angle_Pb_data[i][n])        
        q += 1

    #angle I-Pb-I
    q = 0
    ang_IPbI = []
    for i in range (8, 19):
        for m in range(0, 11-q):
            if angle_Pb_data[i][m] == 0:
                pass
            else:
                ang_IPbI.append(angle_Pb_data[i][m])
        q +=1

    return(ang_CsPbCs, ang_CsPbPb, ang_CsPbI, ang_PbPbPb, ang_PbPbI, ang_IPbI)



# In[20]:


def ang_I(angle_I_data):
    # angle Cs-I-Cs, Cs-I-Pb, Cs-I-I
    ang_CsICs = []
    ang_CsIPb = []
    ang_CsII = []

    for i in range(0,4):
        for m in range (0, 3-i):
            if angle_I_data[i][m] == 0:
                pass
            else:
                ang_CsICs.append(angle_I_data[i][m])

        for n in range(3-i, 7-i):
            if angle_I_data[i][n] == 0:
                pass
            else:
                ang_CsIPb.append(angle_I_data[i][n])

        for k in range(7-i, 19-i):
            if angle_I_data[i][k] == 0:
                pass
            else:
                ang_CsII.append(angle_I_data[i][k])

    #angle Pb-I-Pb, Pb-I-I

    ang_PbIPb = []
    ang_PbII = []
    q = 0

    for i in range(4, 8):
        for m in range(0, 3-q):
            if angle_I_data[i][m] == 0:
                pass          
            else:
                ang_PbIPb.append(angle_I_data[i][m])

        for n in range(3-q, 15-q):
            if angle_I_data[i][n] == 0:
                pass          
            else:
                ang_PbII.append(angle_I_data[i][n])        
        q += 1

    #angle I-I-I
    q = 0
    ang_III = []
    for i in range (8, 19):
        for m in range(0, 11-q):
            if angle_I_data[i][m] == 0:
                pass
            else:
                ang_III.append(angle_I_data[i][m])
        q +=1

    return(ang_CsICs, ang_CsIPb, ang_CsII, ang_PbIPb, ang_PbII, ang_III)




# In[21]:


def thetas_decision(pos):
    step1 = mic_shortest_dis_xyz(pos)
    step2 = mic_shortest_dis(step1)
    step3 = angle_BAC(step2)
    step4 = angle_preprocess(step3)
    
    ang_CsCsCs = []
    ang_CsCsPb = []
    ang_CsCsI = []
    ang_PbCsPb = []
    ang_PbCsI = []
    ang_ICSI = []
    
    for i in range (0,4):
        Cs = ang_Cs(step4[i])
        ang_CsCsCs.append(Cs[0])
        ang_CsCsPb.append(Cs[1])
        ang_CsCsI.append(Cs[2])
        ang_PbCsPb.append(Cs[3])
        ang_PbCsI.append(Cs[4])
        ang_ICSI.append(Cs[5])
        
    ang_CsPbCs = []
    ang_CsPbPb = []
    ang_CsPbI = []
    ang_PbPbPb = []
    ang_PbPbI = []   
    ang_IPbI = []
    
    for i in range (4,8):
        Pb = ang_Pb(step4[i])
        ang_CsPbCs.append(Pb[0])
        ang_CsPbPb.append(Pb[1])
        ang_CsPbI.append(Pb[2])
        ang_PbPbPb.append(Pb[3])
        ang_PbPbI.append(Pb[4])
        ang_IPbI.append(Pb[5])
        
    ang_CsICs = []    
    ang_CsIPb = []   
    ang_CsII = []    
    ang_PbIPb = []    
    ang_PbII = []    
    ang_III = []
        
    for i in range (8, 20):
        I = ang_I(step4[i])
        ang_CsICs.append(I[0])
        ang_CsIPb.append(I[1])
        ang_CsII.append(I[2])
        ang_PbIPb.append(I[3])
        ang_PbII.append(I[4])
        ang_III.append(I[5])
    
    max_ang_CsCsCs = max(np.array(ang_CsCsCs).ravel())
    max_ang_CsCsPb = max(np.array(ang_CsCsPb).ravel())
    max_ang_CsCsI = max(np.array(ang_CsCsI).ravel())
    max_ang_PbCsPb = max(np.array(ang_PbCsPb).ravel())    
    max_ang_PbCsI = max(np.array(ang_PbCsI).ravel())
    max_ang_ICSI = max(np.array(ang_ICSI).ravel())
    
    max_ang_CsPbCs = max(np.array(ang_CsPbCs).ravel())
    max_ang_CsPbPb = max(np.array(ang_CsPbPb).ravel())    
    max_ang_CsPbI = max(np.array(ang_CsPbI).ravel())
    max_ang_PbPbPb = max(np.array(ang_PbPbPb).ravel())    
    max_ang_PbPbI = max(np.array(ang_PbPbI).ravel())
    max_ang_IPbI = max(np.array(ang_IPbI).ravel())
    
    max_ang_CsICs = max(np.array(ang_CsICs).ravel())
    max_ang_CsIPb = max(np.array(ang_CsIPb).ravel())    
    max_ang_CsII = max(np.array(ang_CsII).ravel())
    max_ang_PbIPb = max(np.array(ang_PbIPb).ravel()) 
    max_ang_PbII = max(np.array(ang_PbII).ravel())
    max_ang_III = max(np.array(ang_III).ravel())      
    
    min_ang_CsCsCs = min(np.array(ang_CsCsCs).ravel())
    min_ang_CsCsPb = min(np.array(ang_CsCsPb).ravel())
    min_ang_CsCsI = min(np.array(ang_CsCsI).ravel())
    min_ang_PbCsPb = min(np.array(ang_PbCsPb).ravel())    
    min_ang_PbCsI = min(np.array(ang_PbCsI).ravel())
    min_ang_ICSI = min(np.array(ang_ICSI).ravel())
    
    min_ang_CsPbCs = min(np.array(ang_CsPbCs).ravel())
    min_ang_CsPbPb = min(np.array(ang_CsPbPb).ravel())    
    min_ang_CsPbI = min(np.array(ang_CsPbI).ravel())
    min_ang_PbPbPb = min(np.array(ang_PbPbPb).ravel())    
    min_ang_PbPbI = min(np.array(ang_PbPbI).ravel())
    min_ang_IPbI = min(np.array(ang_IPbI).ravel())
    
    min_ang_CsICs = min(np.array(ang_CsICs).ravel())
    min_ang_CsIPb = min(np.array(ang_CsIPb).ravel())    
    min_ang_CsII = min(np.array(ang_CsII).ravel())
    min_ang_PbIPb = min(np.array(ang_PbIPb).ravel()) 
    min_ang_PbII = min(np.array(ang_PbII).ravel())
    min_ang_III = min(np.array(ang_III).ravel())
    
    max_ang_data = np.zeros((1,18))
    min_ang_data = np.zeros((1,18))
    
    max_ang_data[0][0], max_ang_data[0][1], max_ang_data[0][2],max_ang_data[0][3], max_ang_data[0][4], max_ang_data[0][5] = max_ang_CsCsCs, max_ang_CsCsPb, max_ang_CsCsI, max_ang_PbCsPb, max_ang_PbCsI, max_ang_ICSI
    max_ang_data[0][6], max_ang_data[0][7], max_ang_data[0][8],max_ang_data[0][9], max_ang_data[0][10], max_ang_data[0][11] = max_ang_CsPbCs, max_ang_CsPbPb, max_ang_CsPbI, max_ang_PbPbPb, max_ang_PbPbI, max_ang_IPbI
    max_ang_data[0][12], max_ang_data[0][13], max_ang_data[0][14],max_ang_data[0][15], max_ang_data[0][16], max_ang_data[0][17] = max_ang_CsICs, max_ang_CsIPb, max_ang_CsII, max_ang_PbIPb, max_ang_PbII, max_ang_III   
        
    min_ang_data[0][0], min_ang_data[0][1], min_ang_data[0][2],min_ang_data[0][3], min_ang_data[0][4], min_ang_data[0][5] = min_ang_CsCsCs, min_ang_CsCsPb, min_ang_CsCsI, min_ang_PbCsPb, min_ang_PbCsI, min_ang_ICSI
    min_ang_data[0][6], min_ang_data[0][7], min_ang_data[0][8],min_ang_data[0][9], min_ang_data[0][10], min_ang_data[0][11] = min_ang_CsPbCs, min_ang_CsPbPb, min_ang_CsPbI, min_ang_PbPbPb, min_ang_PbPbI, min_ang_IPbI
    min_ang_data[0][12], min_ang_data[0][13], min_ang_data[0][14],min_ang_data[0][15], min_ang_data[0][16], min_ang_data[0][17] = min_ang_CsICs, min_ang_CsIPb, min_ang_CsII, min_ang_PbIPb, min_ang_PbII, min_ang_III           
        
    return(max_ang_data, min_ang_data)    
        
        
        
        
        
        

        


# In[ ]:


t1 = time()
max_angle = []
min_angle = []

for i in range (0, 10000):
    max_angle.append(np.array(thetas_decision(pos[i])[0]))
    min_angle.append(np.array(thetas_decision(pos[i])[1]))
    
time()-t1    
    
    
    


# In[ ]:


max_angle = np.array(max_angle)
min_angle = np.array(min_angle)
max_avg_angle = np.average(max_angle, axis = 0)
min_avg_angle = np.average(min_angle, axis = 0)
max_min_delta = max_avg_angle - min_avg_angle


# In[ ]:


max_avg_angle 


# In[ ]:


min_avg_angle


# In[ ]:


max_min_delta


# In[ ]:





# In[22]:


def cutoff (dis):
    if dis > Rcut:
        f_cut = 0
    if dis == 0:
        f_cut = 0
    else:
        f_cut = 0.5*cos(pi * dis/Rcut) + 0.5
    
    return f_cut
    
    


# In[23]:


#Cs-I, Cs-Pb, Cs-Cs
Rs_Cs = [3.7514617, 4.96424887, 5.82601472]

#Pb-I, Pb-Cs, Pb-Pb
Rs_Pb = [2.98938695, 4.96424887, 6.17388088]

#I-Pb, I-Cs, I-I
Rs_I = [2.98938695, 3.7514617, 4.10984614]


# In[24]:


#Cs-Cs-Pb, Cs-Cs-Cs, Cs-Cs-I, Pb-Cs-I, Pb-Cs-Pb, I-Cs-I
thetas_Cs = [1.0765401, 1.59292253, 1.66209905, 1.7899382, 2.01433819, 3.01143415]

#Cs-Pb-Pb, Pb-Pb-Pb, Pb-Pb-I, Cs-Pb-Cs, Cs-Pb-I, I-Pb-I
thetas_Pb = [1.03929047, 1.56521586, 1.63728606, 1.9439263, 2.30627383, 3.0252016]

#Cs-I-Pb, Cs-I-I, Pb-I-I, I-I-I, Pb-I-Pb, Cs-I-Cs 
thetas_I = [1.82739854, 2.31958872, 2.57124806, 2.90291669, 2.91610732, 3.00980311]


# In[25]:


Rcut =  9.098997905703929
zeta = 1 #enough to conver ang_del_max = 2.7
#thetas = 0
eta = 0.15 #to conver dis_del_max = 4.7
#Rs = [0.5,1,2,3,4,5,6,7,8,9.1]


# In[27]:


def features_Cs_cal (dis_data, angle_data):
       
    coeff = 2**(1-zeta)
    Rs_data = Rs_Cs
    thetas_data = thetas_Cs
    feature_over_thetas = []
    
    for thetas in thetas_data:
        for Rs in Rs_data:
            feature_over_Rs = []
            for n in range(0, 4): # atom index starting from first Cs
                feature_Cs_over_all_atom = []
                for i in range (0,20): #20 distances in a list, starting from Cs-Cs
                    for m in range(i+1,20):
                        r_ij = dis_data[n][i]
                        r_ik = dis_data[n][m]
                        fc_ij = cutoff(r_ij)
                        fc_ik = cutoff(r_ik)

                        #radial term
                        rad_ftu = exp(-eta * ((r_ij+r_ik-2*Rs)/2)**2)

                        angle_ijk = angle_data[n][i][m-i-1]

                        #angular term
                        if angle_ijk == 0.:
                            ang_ftu = 0
                        else:
                            ang_ftu = (1+cos(angle_ijk - thetas))**zeta

                        feature_Cs = ang_ftu * rad_ftu * fc_ij * fc_ik
                        feature_Cs_over_all_atom.append(feature_Cs)

                Sum = coeff * np.sum(feature_Cs_over_all_atom) #keep it here, all in one atom
                feature_over_Rs.append(Sum)
            feature_over_thetas.append(feature_over_Rs)

    features = np.array(feature_over_thetas)
    
    return features.ravel()
    













# In[28]:


def features_Pb_cal (dis_data, angle_data):
       
    coeff = 2**(1-zeta)
    Rs_data = Rs_Pb
    thetas_data = thetas_Pb
    feature_over_thetas = []
    
    for thetas in thetas_data:
        for Rs in Rs_data:
            feature_over_Rs = []
            for n in range(4, 8): # atom index starting from first Pb
                feature_Pb_over_all_atom = []
                for i in range (0,20): #20 distances in a list, starting from Pb-Cs
                    for m in range(i+1,20):
                        r_ij = dis_data[n][i]
                        r_ik = dis_data[n][m]
                        fc_ij = cutoff(r_ij)
                        fc_ik = cutoff(r_ik)

                        #radial term
                        rad_ftu = exp(-eta * ((r_ij+r_ik-2*Rs)/2)**2)

                        angle_ijk = angle_data[n][i][m-i-1]

                        #angular term
                        if angle_ijk == 0.:
                            ang_ftu = 0
                        else:
                            ang_ftu = (1+cos(angle_ijk - thetas))**zeta

                        feature_Pb = ang_ftu * rad_ftu * fc_ij * fc_ik
                        feature_Pb_over_all_atom.append(feature_Pb)

                Sum = coeff * np.sum(feature_Pb_over_all_atom) #keep it here, all in one atom
                feature_over_Rs.append(Sum)
            feature_over_thetas.append(feature_over_Rs)

    features = np.array(feature_over_thetas)
    
    return features.ravel()


# In[29]:


def features_I_cal (dis_data, angle_data):
       
    coeff = 2**(1-zeta)
    Rs_data = Rs_I
    thetas_data = thetas_I
    feature_over_thetas = []
    
    for thetas in thetas_data:
        for Rs in Rs_data:
            feature_over_Rs = []
            for n in range(8, 20): # atom index starting from first I
                feature_I_over_all_atom = []
                for i in range (0,20): #20 distances in a list, starting from I-Cs
                    for m in range(i+1,20):
                        r_ij = dis_data[n][i]
                        r_ik = dis_data[n][m]
                        fc_ij = cutoff(r_ij)
                        fc_ik = cutoff(r_ik)

                        #radial term
                        rad_ftu = exp(-eta * ((r_ij+r_ik-2*Rs)/2)**2)

                        angle_ijk = angle_data[n][i][m-i-1]

                        #angular term
                        if angle_ijk == 0.:
                            ang_ftu = 0
                        else:
                            ang_ftu = (1+cos(angle_ijk - thetas))**zeta

                        feature_I = ang_ftu * rad_ftu * fc_ij * fc_ik
                        feature_I_over_all_atom.append(feature_I)

                Sum = coeff * np.sum(feature_I_over_all_atom) #keep it here, all in one atom
                feature_over_Rs.append(Sum)
            feature_over_thetas.append(feature_over_Rs)

    features = np.array(feature_over_thetas)
    
    return features.ravel()


# In[30]:


def features_concat (dis_data, angle_data):
    Cs_feat = np.array(features_Cs_cal (dis_data, angle_data))
    Pb_feat = np.array(features_Pb_cal (dis_data, angle_data))
    I_feat = np.array(features_I_cal (dis_data, angle_data))
    
    feat = np.concatenate((Cs_feat, Pb_feat, I_feat))
    return feat
    


# In[31]:


def mod_SF_feature (pos):
    
    step1 = mic_shortest_dis_xyz(pos)
    step2 = mic_shortest_dis(step1)
    
    step3 = angle_BAC(step2)
    step4 = angle_preprocess(step3)
    
    step5 = features_concat(step2, step4)
    
    return step5
    
    


# In[32]:


fingerprint = []

t1 = time()
for i in range(0,10000):
    fingerprint.append(mod_SF_feature(pos[i]))
    
t2 = time()
t_del = t2-t1
t_del


# In[33]:


len(fingerprint)


# In[34]:


fingerprint


# In[35]:


np.save("features_v2.npy", fingerprint)


# In[36]:


pd.DataFrame(fingerprint).to_csv("features_v2.csv")


# In[ ]:





# In[ ]:




