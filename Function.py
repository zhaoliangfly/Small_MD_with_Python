#!/usr/bin/python
from math import *

# P_1 and P_2 are the list including coordinates of two atoms #
# Q_1 and Q_2 are the charges of two atoms                    #
# Sig_Eps_1 and Sig_Eps_2 are the force parameters            #

def Distance(P_1, P_2):
	distance=sqrt(pow((P_1[0]-P_2[0]),2)+pow((P_1[1]-P_2[1]),2)+pow((P_1[2]-P_2[2]),2))
	return distance

def Direction_Vector(P_1, P_2):
	direction_vector=[0.00,0.00,0.00]
	direction_vector[0]=(P_2[0]-P_1[0])/Distance(P_1, P_2)
	direction_vector[1]=(P_2[1]-P_1[1])/Distance(P_1, P_2)
	direction_vector[2]=(P_2[2]-P_1[2])/Distance(P_1, P_2)
	return direction_vector

def Ele_Force(Q_1, Q_2, P_1, P_2):
	ele_force=138.935485*Q_1*Q_2/pow(Distance(P_1, P_2),2)
	return ele_force

def Ele_Energy(Q_1, Q_2, P_1, P_2):
	ele_energy=138.935485*Q_1*Q_2/Distance(P_1, P_2)
	return ele_energy

def Vdw_Force(Sig_Eps_1, Sig_Eps_2, P_1, P_2):
	Sig_1_2=0.5*(Sig_Eps_1[0]+Sig_Eps_2[0])
	Eps_1_2=sqrt(Sig_Eps_1[1]+Sig_Eps_2[1])
	vdw_force=4.0*Eps_1_2*(pow(Sig_1_2,12)*(-12)/pow(Distance(P_1, P_2),13)+pow(Sig_1_2,6)*6/pow(Distance(P_1, P_2),7))
	return vdw_force

def Vdw_Energy(Sig_Eps_1, Sig_Eps_2, P_1, P_2):
	Sig_1_2=0.5*(Sig_Eps_1[0]+Sig_Eps_2[0])
	Eps_1_2=sqrt(Sig_Eps_1[1]+Sig_Eps_2[1])
	vdw_energy=4.0*Eps_1_2*(pow(Sig_1_2/Distance(P_1, P_2),12)-pow(Sig_1_2/Distance(P_1, P_2),6))
	return vdw_energy

def Bond_Force(kb, a0, P_1, P_2):
	bond_force=kb*(Distance(P_1, P_2)-a0)
	return bond_force

def Bond_Energy(kb, a0, P_1, P_2):
	bond_energy=0.5*kb*pow(Distance(P_1, P_2)-a0,2)
	return bond_energy

def Kinetic_Energy(m,v):
	kinetic_energy=0.5*m*pow(v[0],2)+0.5*m*pow(v[1],2)+0.5*m*pow(v[2],2)
	return kinetic_energy
# The codes below decide the algorithm we choose           #
