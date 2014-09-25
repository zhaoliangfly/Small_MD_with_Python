#!/usr/bin/python
from Function import *

# We have only 2 atoms #
dt=0.0001
N=10000
m=[14.005,14.005]

Q=[0.00,0.00]

Sig_Eps_1=[0.0,0.00]
Sig_Eps_2=[0.0,0.00]
Sig_Eps=[Sig_Eps_1,Sig_Eps_2]

p_1=[-0.40,0.00,0.00]
p_2=[+0.40,0.00,0.00]
p=[p_1,p_2]

v_1=[0.00,0.00,0.00]
v_2=[0.00,0.00,0.00]
v=[v_1,v_2]

a_1=[0.00,0.00,0.00]
a_2=[0.00,0.00,0.00]
a=[a_1,a_2]

step=0

while step < N :
	step = step +1
# To calculate the a(t) #	
	for i in [0,1]:
		for j in [0,1]:
			if j!=i :
				for k in [0,1,2]:
					a[i][k]+= Ele_Force(Q[i],Q[j],p[i],p[j])*Direction_Vector(p[i],p[j])[k]/m[i]+Vdw_Force(Sig_Eps[i],Sig_Eps[j],p[i],p[j])*Direction_Vector(p[i],p[j])[k]/m[i]+Bond_Force(15000,0.142,p[i],p[j])*Direction_Vector(p[i],p[j])[k]/m[i]
# To calculate the r(t+dt) #
		for l in [0,1,2]:
			p[i][l]+=dt*v[i][l]+0.5*dt*dt*a[i][l]
# Reset the a(t) #
	for i in [0,1]:
		for k in [0,1,2]:
			a[i][k]=0.00;

# To calculate the a(t+dt) #			
	for i in [0,1]:
		for j in [0,1]:
			if j!=i :
				for k in [0,1,2]:
					a[i][k]+= Ele_Force(Q[i],Q[j],p[i],p[j])*Direction_Vector(p[i],p[j])[k]/m[i]+Vdw_Force(Sig_Eps[i],Sig_Eps[j],p[i],p[j])*Direction_Vector(p[i],p[j])[k]/m[i]+Bond_Force(15000,0.5,p[i],p[j])*Direction_Vector(p[i],p[j])[k]/m[i]

        b=[[0.00,0.00,0.00],[0.00,0.00,0.00]]
	for i in [0,1]:
		for k in [0,1,2]:
			b[i][k]=a[i][k]

# To calculate the v(t+dt) #
	for i in [0,1]:
		for k in [0,1,2]:
			v[i][k]+=0.5*dt*(a[i][k]+b[i][k])

#	print step,p[0][0],p[1][0],v[0][0],v[1][0]
	print step,Bond_Energy(15000,0.5,p[1],p[0]),Kinetic_Energy(m[0],v[0]),Kinetic_Energy(m[1],v[1])



