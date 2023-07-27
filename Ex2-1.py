import numpy as np
import scipy.linalg as linalg
import argparse
import copy
import MPS
import matplotlib 
import matplotlib.pyplot as plt


## Set parameters
N = 16 ## set "system size" N 
m = 2 ## vector size m: total dimension is m^N
chi_max = 20 ## maximum bond dimension at truncation
seed = None ## The seed for random numnber generator. 

class model_cal():
    def __init__(self,N,m,chi_max,seed):
        self.N = N
        self.m = m
        self.chi_max = chi_max
        self.seed = seed


    def calc_innerproduct(self,Tn1,lam1,Tn2,lam2):


        chi1 = Tn1[0].shape[2]
        chi2 = Tn2[0].shape[2]

        vec = np.tensordot(Tn1[0],Tn2[0].conj(),axes=(0,0)).reshape(chi1,chi2)

        for i in range(1,len(Tn1)):
            vec = np.tensordot(np.tensordot(np.tensordot(np.tensordot(vec,np.diag(lam1[i]),(0,0)),np.diag(lam2[i]),(0,0)),Tn1[i],(0,1)),Tn2[i].conj(),([0,1],[1,0]))

        return vec.reshape(1)[0]


    def remake_vec(self,Tn,lam):

        chi = Tn[0].shape[2]
        m = Tn[0].shape[0]
        vec = np.reshape(Tn[0],(m,chi))

        for i in range(1,len(Tn)):
            vec = np.tensordot(np.tensordot(vec,np.diag(lam[i]),(i,0)),Tn[i],(i,1))
        return vec.flatten()


    ## Main calculation
    def make_vector(self):
        if seed != None:
            np.random.seed(seed)
                
            
        print("Parameters: N, m, chi_max = "+repr(self.N)+", "+repr(self.m)+ ", "+repr(self.chi_max))
        print("Random seed: = "+repr(self.seed))

        ## create random vecgtor
        eig_vec = ((np.random.rand(self.m**self.N)-0.5) + 1.0j * (np.random.rand(self.m**self.N)-0.5)).reshape(self.m**self.N)
        ## normalization
        norm = np.tensordot(eig_vec,eig_vec.conj(),axes=(0,0))
        eig_vec /= np.sqrt(np.abs(norm))

        return eig_vec


    ## Make exact MPS (from "left")
    def make_MPS_left(self,eig_vec):
        Tn_ex = []
        lam_ex = [np.ones((1,))]
        lam_inv = 1.0/lam_ex[0]
        R_mat = eig_vec[:].reshape(self.m,self.m**(self.N-1))

        chi_l=1
        for i in range(N-1):
            U,s,VT = linalg.svd(R_mat,full_matrices=False)
            chi_r = s.size

            Tn_ex.append(np.tensordot(np.diag(lam_inv),U.reshape(chi_l,self.m,chi_r),(1,0)).transpose(1,0,2))
            lam_ex.append(s)
            lam_inv = 1.0/s
            R_mat = np.dot(np.diag(s),VT).reshape(chi_r*self.m,self.m**(self.N-i-2))
            chi_l = chi_r
        Tn_ex.append(VT.reshape(self.m,self.m,1).transpose(1,0,2))
        lam_ex.append(np.ones((1,)))

        return Tn_ex, lam_ex

    ## Truncation to chi_max
    def truncation(self,Tn_ex,lam_ex):
        # Tn_ex = copy.deepcopy(Tn)
        # lam_ex = copy.deepcopy(lam)

        #Tn_ex = Tn
        #lam_ex = lam
        # Tn_exとlam_exをTruncationしてTn_ex = Tn, lam_ex = lamとする
        # これを実行するとTn=Tn_exとなってしまうため、リターンでもう一度Tn_exとlam_exを返す。
        Tn_ex = copy.deepcopy(Tn_ex)
        lam_ex = copy.deepcopy(lam_ex)
        for i in range(self.N-1):
            chi = min(self.chi_max,lam_ex[i+1].shape[0])
            lam_ex[i+1]=lam_ex[i+1][:chi]
            Tn_ex[i]=Tn_ex[i][:,:,:chi]
            Tn_ex[i+1]=Tn_ex[i+1][:,:chi,:]

        print("Truncation: chi_max = "+repr(chi_max))
        Tn = Tn_ex
        lam = lam_ex
        
        return Tn_ex, lam_ex, Tn, lam

    ## plot Schmidt coefficient at N/2
    ## Red line indicates the position of chi_max
    def plot_schmidt_coefficient_at_half(self,lam_ex):
        plt.title("Schmidt coefficients for "+"(N, m) = ("+repr(self.N)+", "+repr(self.m)+") random vector")
        plt.plot(np.arange(len(lam_ex[self.N//2]))+1,lam_ex[self.N//2]**2,"o",label="Schmidt coefficients")
        plt.axvline([chi_max],0,1,  c="red", linestyle='dashed', label="chi_max") ## position of chi_max
        plt.xlabel("index")
        plt.xscale("log")
        plt.ylabel("Schmidt coefficients")
        plt.yscale("log")
        plt.legend()
        plt.show()

    ## Distance between Exact MPS and truncated MPS
    def distance(self,Tn_ex, lam_ex, Tn, lam):
        vec_ex = MPS.remake_vec(Tn_ex,lam_ex)
        vec_ap = MPS.remake_vec(Tn,lam)
        distance = linalg.norm(vec_ex - vec_ap)
        print("Distance between exact and truncated MPS = "+repr(linalg.norm(vec_ex - vec_ap)))

        return distance
    
    # chiの違いでどのくらいdistanceが変わるのかを図示
    def compare_chi(self,min_chi_max,max_chi_max,d_chi_max,Tn_ex,lam_ex):
        chi_max_list = np.arange(min_chi_max, max_chi_max+1, d_chi_max, dtype=int)
        chi_list = np.ones((self.N+1,),dtype=int)
        vec_ex = MPS.remake_vec(Tn_ex,lam_ex)

        distances=[]
        for chi_max in chi_max_list:
            for i in range(self.N-1):
                chi = min(chi_max,lam_ex[i+1].shape[0])
                chi_list[i+1] = chi
            lam = [np.ones((1,))]
            Tn = []
            for i in range(self.N):
                lam.append(lam_ex[i+1][:chi_list[i+1]])
                Tn.append(Tn_ex[i][:,:chi_list[i],:chi_list[i+1]])
            vec_ap = MPS.remake_vec(Tn,lam)
            distances.append(linalg.norm(vec_ex - vec_ap))       
        ## plot distances

        plt.title("Distances for "+"(N, m) = ("+repr(N)+", "+repr(m)+") random vector")
        plt.plot(chi_max_list,distances,"o")
        plt.xlabel("chi_max")
        plt.ylabel("Distance")
        #plt.yscale("log")
        plt.show()

    
calc = model_cal(N,m,chi_max,seed)
eig_vec = calc.make_vector()
Tn_ex, lam_ex = calc.make_MPS_left(eig_vec)
Tn_ex1, lam_ex1,Tn, lam = calc.truncation(Tn_ex, lam_ex)
calc.plot_schmidt_coefficient_at_half(lam_ex)
distacne = calc.distance(Tn_ex, lam_ex, Tn, lam)
calc.compare_chi(10,260,10,Tn_ex,lam_ex)