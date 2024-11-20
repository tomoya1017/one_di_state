import numpy as np
import scipy.linalg as linalg
import argparse
import copy
import MPS
import matplotlib 
import matplotlib.pyplot as plt
import ED
from basic import basic_MPS


class model_cal(basic_MPS):
    def __init__(self,N,m,chi_max,seed,Jz,Jxy,hx,D):
        super().__init__(N,m,chi_max,seed)

        self.Jz = Jz
        self.Jxy = Jxy
        self.hx = hx
        self.D = D

    def make_model(self):
        print("Model parameters: Jz, Jxy, hx = "+repr(self.Jz)+", "+repr(self.Jxy)+", "+repr(self.hx))
        print("Parameters: N, m, chi_max = "+repr(self.N)+", "+repr(self.m)+ ", "+repr(self.chi_max))

        
        eig_val,eig_vec = ED.Calc_GS(self.m,self.Jz,self.Jxy,self.hx,self.D,self.N,k=1)


        print("Ground state energy per bond = "+repr(eig_val[0]/(self.N-1)))

        return eig_val, eig_vec
    
    def calc_exact_energy(self,Tn_ex,lam_ex):
        Env_left=[]
        Env_right=[]
        for i in range(self.N):
            Env_left.append(np.identity((lam_ex[i].shape[0])))
            Env_right.append(np.dot(np.dot(np.diag(lam_ex[i+1]),np.identity((lam_ex[i+1].shape[0]))),np.diag(lam_ex[i+1])))

        # Tn_ex = copy.deepcopy(Tn)
        # lam_ex = copy.deepcopy(lam)
        E_exact = MPS.Calc_Energy(Env_left,Env_right,Tn_ex,lam_ex,self.Jz,self.Jxy,self.hx,self.D)
        print("Energy of Exact MPS = "+repr(E_exact))

        return  E_exact
    
    def calc_truncated_energy(self,Tn,lam,):
        ## Calculate Energy for truncated MPS
        Env_left=[]
        Env_right=[]
        for i in range(self.N):
            Env_left.append(np.identity((lam[i].shape[0])))
            Env_right.append(np.dot(np.dot(np.diag(lam[i+1]),np.identity((lam[i+1].shape[0]))),np.diag(lam[i+1])))

        print("Truncation: chi_max = "+repr(self.chi_max))
        E_truncated = MPS.Calc_Energy(Env_left,Env_right,Tn,lam,self.Jz,self.Jxy,self.hx,self.D)
        print("Energy of MPS with truncation = "+repr(E_truncated))

        return E_truncated
        

    def difference(self,E_exact,E_truncated):
        print("Energy difference: E_truncated - E_exact =" + repr(E_truncated - E_exact))

    def compare(self,min_chi_max, max_chi_max,d_chi_max,Tn_ex,lam_ex,E_exact):
        chi_max_list = np.arange(min_chi_max, max_chi_max+1, d_chi_max, dtype=int)
        chi_list = np.ones((self.N+1,),dtype=int)
        vec_ex = MPS.remake_vec(Tn_ex,lam_ex)

        distances=[]
        energies=[]
        for chi_max in chi_max_list:
            for i in range(self.N-1):
                chi = min(chi_max,lam_ex[i+1].shape[0])
                chi_list[i+1] = chi
            lam = [np.ones((1,))]
            Tn = []
            for i in range(N):
                lam.append(lam_ex[i+1][:chi_list[i+1]])
                Tn.append(Tn_ex[i][:,:chi_list[i],:chi_list[i+1]])
            vec_ap = MPS.remake_vec(Tn,lam)
            distances.append(linalg.norm(vec_ex - vec_ap))
            
            ## Calculate Energy for truncated MPS
            Env_left=[]
            Env_right=[]
            for i in range(self.N):
                Env_left.append(np.identity((lam[i].shape[0])))
                Env_right.append(np.dot(np.dot(np.diag(lam[i+1]),np.identity((lam[i+1].shape[0]))),np.diag(lam[i+1])))
                
            energies.append(MPS.Calc_Energy(Env_left,Env_right,Tn,lam,self.Jz,self.Jxy,self.hx,self.D))

        ## plot distances

        plt.title("Distances for "+repr(self.N)+" sites spin chain")
        plt.plot(chi_max_list,distances,"o")
        plt.xlabel("chi_max")
        plt.ylabel("Distance")
        plt.yscale("log")
        plt.show()

        plt.title("Energy for "+repr(self.N)+" sites spin chain")
        plt.plot(chi_max_list,(energies-E_exact)/np.abs(E_exact),"o")
        plt.xlabel("chi_max")
        plt.ylabel("(Energy - $E_{ex}$)/$|E_{ex}|$")
        plt.yscale("log")
        plt.show()

    
## Set parameters
N = 16 ## set "system size" N 
m = 2 ## vector size m: total dimension is m^N
chi_max = 2 ## maximum bond dimension at truncation
seed = None ## The seed for random numnber generator. 
Jz = -1.0 ## SzSz interaction
Jxy = 0.0 ## SxSx and SySy interactions
hx = 0.4 ## extarnal transverse magnetic field
D=0.0

model = model_cal(N,m,chi_max,seed,Jz,Jxy,hx,D)
eig_val, eig_vec = model.make_model()
Tn_ex, lam_ex = model.make_MPS_left(eig_vec)
E_exact = model.calc_exact_energy(Tn_ex,lam_ex)
model.plot_schmidt_coefficient_at_half(lam_ex)
Tn_ex1, lam_ex1,Tn, lam = model.truncation(Tn_ex, lam_ex)
E_truncated = model.calc_truncated_energy(Tn,lam)
model.difference(E_exact,E_truncated)
model.distance(Tn_ex,lam_ex,Tn,lam)
model.compare(1,20,1,Tn_ex,lam_ex,E_exact)


