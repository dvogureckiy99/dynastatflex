import sympy as sm
import inspect
from IPython.display import display, Math
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import time
import bisect
import os 

# https://www.geeksforgeeks.org/inheritance-in-python-inner-class/

"""
Idea:
1. add saving files with parameters with different names _1,_2,_3 e.c. with checking
 if this number file_name already exist.
2. when loading: open in sequence files until file with appropriate FEM and Ldivide parameters will be found.
    2.1 discover all files in root folder and finding with right name
    2.2 checking all files with different names _1,_2,_3 e.c. until file 
    with appropriate FEM and Ldivide parameters will be found
"""

class Flex_beam(object):
    def __init__(self,L=0.06,E=2.95*1e9,h=0.01,w=0.001):
        """
        Creating Flex_beam Class instance
        # Parameters
        ----------
        L: float, optional
            beam length in [m]
        """
        # converting all in [mm] and applied multiplier
        self.mult = 1
        self.L = L * 1e3 * self.mult
        self.E = E * 1e-6 /self.mult**2
        self.h = h * 1e3 * self.mult
        self.w = w * 1e3 * self.mult
        self.rho = 1270 * 1e-9 / self.mult**3

    def Create_Simulation(self,disp=False):
        """
        Creating Child Class Simulation instance and passing parameters to it 
        # Parameters
        ---------- 
        disp: bool, optional
            display data
        """
        self.Simulating = self.Simulating(self.L,self.E,self.h,self.w,self.rho,self.mult,disp) # creating Child class instance and passing parameters to it

    def FEM(self,Ne=10,disp=False):
        """
        Splitting a beam into finite elements
        # Parameters
        ---------- 
        Ne: int
            element's number,
        disp: bool, optional
            display data
        """
        self.Ne = Ne
        self.dl = self.L/Ne
        self.Ldl = np.linspace(0,self.L,self.Ne+1)
        try:
            self.Simulating.set_FEM_data(Ne,self.dl,self.Ldl) # inheritance of parameters by a child class
        except:
            raise ValueError("First call Create_Simulation method!") from None
        if disp:
            display(Math("\mathcal{N}_{0,\dots,Ne}=\text{"+np.str_(self.Ldl)+" [mm]}"))   
            
    class Simulating(object):
        def __init__(self,L,E,h,w,rho,mult,disp): # inheritance of parameters by a child class
            self.L = L
            self.E = E
            self.h = h
            self.w = w
            self.I = w**3*h/12
            self.rho = rho
            self.A = w*h
            self.mult = mult
            C_val = np.array([[-6,15,-10,0,0,1],[-3,8,-6,0,1,0],[-0.5,1.5,-1.5,0.5,0,0],
                  [6,-15,10,0,0,0],[-3,7,-4,0,0,0],[0.5,-1,0.5,0,0,0]])
            self.p = []
            self.dp = []
            self.ddp = []
            self.dddp = []
            self.ddddp = []
            for i in range(6):
                self.p.append(np.poly1d(C_val[i])) # for psi
                self.dp.append(np.polyder(self.p[-1], m=1)) # for d^1psi/dlambda^1
                self.ddp.append(np.polyder(self.p[-1], m=2)) # for d^2psi/dlambda^2
                self.dddp.append(np.polyder(self.p[-1], m=3)) # for d^3psi/dlambda^3
                self.ddddp.append(np.polyder(self.p[-1], m=4)) # for d^4psi/dlambda^4
            if disp:
                a = np.ones((1,6))[0]
                # a = np.array([0,1,12,1,1,13])
                x = np.arange(0,1+1e-3,1e-3)
                plt.subplots(1,5,figsize = (20,8))
                plt.subplot(231)
                for (obj,i) in zip(self.p,np.arange(6)):
                    y = np.polyval(obj,x)*a[i]
                    plt.plot(x,y)
                plt.grid(True)
                plt.title("$\psi_i$ - basis functions for 1 FE")
                plt.subplot(232)
                for (obj,i) in zip(self.dp,np.arange(6)):
                    y = np.polyval(obj,x)*a[i]
                    plt.plot(x,y)
                plt.grid(True)
                plt.title("$\\frac{\partial \psi_i(\lambda)}{\partial \lambda}$ - basis functions for 1 FE")
                plt.subplot(233)
                for (obj,i) in zip(self.ddp,np.arange(6)):
                    y = np.polyval(obj,x)*a[i]
                    plt.plot(x,y)
                plt.grid(True)
                plt.title("$\\frac{\partial^2 \psi_i(\lambda)}{\partial \lambda^2}$ - basis functions for 1 FE")
                plt.subplot(234)
                y = np.zeros((1,len(x)))[0]
                for obj,i in zip(self.p,np.arange(6)):
                    y = y + np.polyval(obj,x)*a[i]
                plt.plot(x,y)
                plt.grid(True)
                plt.title("$\sum_i\psi_i$ - sum of basis functions for 1 FE")
                plt.subplot(235)
                y = np.zeros((1,len(x)))[0]
                for obj,i in zip(self.dp,np.arange(6)):
                    y = y + np.polyval(obj,x)*a[i]
                plt.plot(x,y)
                plt.grid(True)
                plt.title("$\sum_i \\frac{\partial \psi_i(\lambda)}{\partial \lambda}$ - sum of basis functions for 1 FE")
                plt.subplot(236)
                y = np.zeros((1,len(x)))[0]
                for obj,i in zip(self.ddp,np.arange(6)):
                    y = y + np.polyval(obj,x)*a[i]
                plt.plot(x,y)
                plt.grid(True)
                plt.title("$\sum_i \\frac{\partial^2 \psi_i(\lambda)}{\partial \lambda^2}$ - sum of basis functions for 1 FE")
                plt.show()

        def __bmatrix(self,a):
            """
            Returns a LaTeX bmatrix
            # Parameters
            -----------    
            :a: numpy array
            :returns: LaTeX bmatrix as a string
            """
            if len(a.shape) > 2:
                raise ValueError('bmatrix can at most display two dimensions')
            a = np.round(a,3)
            lines = np.array2string(a,formatter={'float':lambda x: "%.4f" % x}).replace('\n  ', ' ').replace('[', '').replace(']', '').splitlines()
            rv = [r'\begin{bmatrix}']
            rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
            rv +=  [r'\end{bmatrix}']
            return '\n'.join(rv)

        def __delta1(self,l):
            if l>=0:
                return 1
            else:
                return 0

        # def __get_x_approx(self,l):
        #     return self.int_cumsum_x[self.__search_index(self.l_all_true,l)]
        # def __get_y_approx(self,l):
        #     return self.int_cumsum_y[self.__search_index(self.l_all_true,l)]
        # def __get_dx_approx(self,l):
        #     return self.cos_phi_appr[self.__search_index(self.l_all_true,l)]
        # def __get_dy_approx(self,l):
        #     return self.sin_phi_appr[self.__search_index(self.l_all_true,l)]
        # def __get_ddx_approx(self,l):
        #     return self.ddx_sum[self.__search_index(self.l_all_true,l)]
        # def __get_ddy_approx(self,l):
        #     return self.ddy_sum[self.__search_index(self.l_all_true,l)]

        def __fun_static_optim(self,a_diff,F_perp_ext,dl_F_perp_ext):
            # preparing a vector for each FE, cause a_diff contain only unique values 
            a = np.zeros((1,6*self.Ne))[0]
            a[0] = 0
            a[1:6] = a_diff[:5]
            for i in range(self.Ne-1):
                if i==self.Ne-2:
                    a[6*(i+1):6*(i+1)+3] = a[6*(i+1)-3:6*(i+1)]
                    a[6*(i+1)+3:6*(i+2)] = np.concatenate([a_diff[5+3*i:],np.array([0,0])])
                else:
                    a[6*(i+1):6*(i+1)+3] = a[6*(i+1)-3:6*(i+1)]
                    a[6*(i+1)+3:6*(i+2)] = a_diff[5+3*i:8+3*i]
            
            dphi_appr = np.matmul(self.dpsi,a) # [N,1]
            
            
            # self.ddphi_appr = np.matmul()
            # self.dddphi_appr = np.array([])
            # self.ddddphi_appr = np.array([])
            # self.cos_phi_appr = np.array([])
            # self.sin_phi_appr = np.array([])
            # self.cos_phi_appr = np.append(self.cos_phi_appr,np.cos(self.phi_appr[-1]))
            # self.sin_phi_appr = np.append(self.sin_phi_appr,np.sin(self.phi_appr[-1]))

            # for x,y constraints
            # for x 
            # sin_phi_appr_ddphi_appr = np.multiply(self.sin_phi_appr,self.ddphi_appr)
            # sin_phi_appr_ddphi_appr_psi_L = sin_phi_appr_ddphi_appr[-1]*self.psi_sum[-1]
            # sin_phi_appr_ddphi_appr_psi_0 = sin_phi_appr_ddphi_appr[0]*self.psi_sum[0]
            # M3_x = -np.sum(np.multiply(sin_phi_appr_ddphi_appr,self.dpsi_sum)) + sin_phi_appr_ddphi_appr_psi_L - sin_phi_appr_ddphi_appr_psi_0
            # for y
            # cos_phi_appr_ddphi_appr = np.multiply(self.cos_phi_appr,self.ddphi_appr)
            # cos_phi_appr_ddphi_appr_psi_L = cos_phi_appr_ddphi_appr[-1]*self.psi_sum[-1]
            # cos_phi_appr_ddphi_appr_psi_0 = cos_phi_appr_ddphi_appr[0]*self.psi_sum[0]
            # M3_y = np.sum(np.multiply(cos_phi_appr_ddphi_appr,self.dpsi_sum)) - cos_phi_appr_ddphi_appr_psi_L + cos_phi_appr_ddphi_appr_psi_0

            M1 = -np.sum(np.multiply(self.ddddphi_appr,self.psi_sum))+\
                np.sum(np.multiply(np.multiply(self.dphi_appr**2,self.ddphi_appr),self.psi_sum))
            
            # M2 =    self.psi_sum_F_ext_in*self.__delta1((self.L-dl_F_perp_ext)/2)-\
            #         self.psi_sum_F_ext_out*self.__delta1((self.L+dl_F_perp_ext)/2)
            self.iteration_num += 1   
            # cost=np.sum((M1*self.c1+M2)**2)
            cost = np.array([M1*self.EI+self.int_psi_dF_ext,M3_x,M3_y])
            # cost=np.sum((M1*self.c1+self.int_psi_sum_F_ext)**2+()**2)
            print("iter = {}\ncost={}".format(self.iteration_num,np.sum(cost)))
            return cost
            

        def static_preparing(self,disp=True):
            try:
                self.Ne
                self.dl
                self.Ldl
            except:
                raise ValueError("No FEM formulation was made!") from None
            try:
                self.step
                self.l_all_true
                self.N
            except:
                raise ValueError("Call Ldivide first!") from None
            
            flag_preparing_already_done = 0
            if os.path.isfile('psi_matrix_and_vectors.npz'):
                flag_preparing_already_done = 1

            if flag_preparing_already_done:
                print("Found numpy zip archive with preparing data: psi vectors, F,M matrix. Checking if we can use it!")
                with np.load('psi_matrix_and_vectors.npz') as npzfile: # for closign after using it
                    self.c1 = npzfile['c1']
                    self.EI = npzfile['EI']
                    self.psi = npzfile['psi']
                    self.dpsi = npzfile['dpsi']
                    self.ddpsi = npzfile['ddpsi']
                    self.dddpsi = npzfile['dddpsi']
                    self.ddddpsi = npzfile['ddddpsi']
                    self.F = npzfile['F']
                    self.M = npzfile['M']
                    N = npzfile['N']
                    Ne = npzfile['Ne']
                    dl = npzfile['dl']
                    step = npzfile['step']
                del npzfile
            
            if (not flag_preparing_already_done) or (not N==self.N) or (not Ne==self.Ne) or (not dl==self.dl) or (not step==self.step):
                if flag_preparing_already_done:
                    print("Checking finished. We cannot use this data as FEM or/and Ldivide parameters mismatch. Creating new one:")
                self.c1 = self.E*self.I/(self.rho*self.A)
                self.EI = self.E*self.I
                start_time = time.time_ns()
                time.sleep(0.000001) # sleep 1 us
                # preparing for fast computation next
                self.psi = np.zeros((self.N,6*self.Ne))
                self.dpsi = np.zeros((self.N,6*self.Ne))
                self.ddpsi = np.zeros((self.N,6*self.Ne))
                self.dddpsi = np.zeros((self.N,6*self.Ne))
                self.ddddpsi = np.zeros((self.N,6*self.Ne))
                for (l,i) in zip(self.l_all_true,range(self.N)):  
                    self.psi[i] = self.__get_psi(l)
                    self.dpsi[i] = self.__get_dpsi(l)
                    self.ddpsi[i] =self.__get_ddpsi(l)
                    self.dddpsi[i] =self.__get_dddpsi(l)
                    self.ddddpsi[i] =self.__get_ddddpsi(l)
                
                time_psi_calc = time.time_ns()-start_time-1*1e3
                print("Psi calculation time: %s s" % (round(time_psi_calc*1e-9,3)))

                self.F = np.zeros((6*self.Ne,6*self.Ne))
                for j in range(6*self.Ne):
                    for i in range(6*self.Ne):
                        self.F[j][i]= np.sum( np.multiply( self.ddpsi[:,i],self.ddpsi[:,j] )*self.step ) +\
                            self.dddpsi[-1,i]*self.psi[-1,j]-self.dddpsi[0,i]*self.psi[0,j]-\
                            self.ddpsi[-1,i]*self.dpsi[-1,j]+self.ddpsi[0,i]*self.dpsi[0,j]
                
                self.M = np.zeros((6*self.Ne,6*self.Ne))
                for j in range(6*self.Ne):
                    for i in range(6*self.Ne):
                        self.M[j][i]= np.sum( np.multiply( self.psi[:,i],self.psi[:,j] )*self.step )

                time_end = time.time_ns()-start_time-1*1e3
                if (time_end-time_psi_calc)==0:
                    print("Psi matrix and vectors calculation time is less then 1 ns")
                else:
                    print("Psi matrix and vectors calculation time: %s s" % (round((time_end-time_psi_calc)*1e-9,3)))
                print("Preparing time: %s s" % (round(time_end*1e-9,3)))
            
                np.savez('psi_matrix_and_vectors.npz',psi=self.psi,dpsi=self.dpsi,\
                     ddpsi=self.ddpsi,dddpsi=self.dddpsi,ddddpsi=self.ddddpsi,F=self.F,\
                        c1=self.c1,EI=self.EI,M=self.M,\
                        N=self.N,Ne=self.Ne,step=self.step,dl=self.dl)
            else:
                if flag_preparing_already_done:
                    print("Checking finished. Using loaded data")

            if disp:
                display(Math("\\bm{F}="+self.__bmatrix(self.F)))
                display(Math("\\bm{M}="+self.__bmatrix(self.M)))

        def static(self,F_perp_ext=1,dl_F_perp_ext=1e-2,a0=[1,2]):
            # # self.psi_sum = np.sum(self.psi,axis=1)
            # # self.dpsi_sum = np.sum(self.dpsi,axis=1)
            # # self.dF_ext_const = F_perp_ext/dl_F_perp_ext/2
            # # self.dF_ext = np.zeros((1,self.N))[0]
            # # for i in range(self.N):
            # #     self.dF_ext[i] = self.dF_ext_const*self.__delta1(self.l_all_true[i]-((self.L-dl_F_perp_ext)/2))-\
            # #         (self.dF_ext_const*2)*self.__delta1(self.l_all_true[i]-((self.L)/2))+\
            # #         (self.dF_ext_const)*self.__delta1(self.l_all_true[i]-((self.L+dl_F_perp_ext)/2))
                
            # # plt.plot(self.l_all_true,self.dF_ext)
            # # plt.xlim((self.L-dl_F_perp_ext)/2,(self.L+dl_F_perp_ext)/2)
            # # plt.show()
            # self.int_psi_dF_ext = np.sum(np.multiply(self.psi_sum,self.dF_ext))
            # # display(Math("\\bm{F}_1="+self.__bmatrix(F1)))

            # self.iteration_num = 0
            # if np.shape(a0)[0]<3:
            #     a0 = np.ones((1,6+3*(self.Ne-1)-1-2))[0]
            #     # a0 = np.array([1.39805362,0.52238212,0.98447721,-0.24425985,0.25443791,1.00086218,-0.24874274,0.25506558,0.99999917,
            #     #                -0.24857425,0.25505484,0.99999997,-0.24857546,0.25505483,1.00000002,-0.24857545,0.25505484,1.00000026,
            #     #                -0.24857441,0.25505484,0.99994755,-0.24858544,0.25505376,0.99399755,-0.24949161,0.25491721,1.20695856,
            #     #                -0.1501284,0.26941578,2.049172,1.17788224,0.71017331])
            # """
            # bound_min = np.zeros((1,len(a0)))[0]
            # bound_max = np.zeros((1,len(a0)))[0]
            # bound_min[0]=-5
            # bound_max[0]=5
            # bound_min[1]=-10
            # bound_max[1]=10
            # for i in range(len(a0)-2-1):
            #     if i%3 == 0:
            #         bound_min[i+2]=-np.pi
            #         bound_max[i+2]=np.pi
            #     if i%3 == 1:
            #         bound_min[i+2]=-5
            #         bound_max[i+2]=5
            #     if i%3 == 2:
            #         bound_min[i+2]=-10
            #         bound_max[i+2]=10
            # bound_min[-1]=-np.pi
            # bound_max[-1]=np.pi
            # """
            
            # start_time = time.time()
            # # res = sp.optimize.minimize(self.__fun_static_optim, a0,args=(F_perp_ext,dl_F_perp_ext),method='Nelder-Mead')
            # res = sp.optimize.least_squares(self.__fun_static_optim,a0,args=(F_perp_ext,dl_F_perp_ext),
            #                                 ftol=1e-8,gtol=1e-8,xtol=1e-20,max_nfev=1e6,method='trf')
            # end_time = time.time()-start_time
            # print("status: %s"%(res.message))
            # print("status: %s"%(res.status))
            # print("evaluation time:%s s" % (round(end_time,0)))
            # print("time on 1 iter:%s ms" % (round(1e3*end_time/self.iteration_num,0)))  
            # print("iteration number:%s" % (self.iteration_num))

            # self.a_diff = np.ones((1,6+3*(self.Ne-1)-1-2))[0]
            # for i in range(6+3*(self.Ne-1)-1-2):
            #     self.a_diff[i] = res.x[i]

            # self.a_approx = np.zeros((1,6*self.Ne))[0]
            # self.a_approx[0] = 0
            # self.a_approx[1:6] = self.a_diff[:5]
            # for i in range(self.Ne-1):
            #     if i==self.Ne-2:
            #         self.a_approx[6*(i+1):6*(i+1)+3] = self.a_approx[6*(i+1)-3:6*(i+1)]
            #         self.a_approx[6*(i+1)+3:6*(i+2)] = np.concatenate([ self.a_diff[5+3*i:],np.array([0,0])])
            #     else:
            #         self.a_approx[6*(i+1):6*(i+1)+3] = self.a_approx[6*(i+1)-3:6*(i+1)]
            #         self.a_approx[6*(i+1)+3:6*(i+2)] = self.a_diff[5+3*i:8+3*i]

            # print("cost = {}".format(res.cost))  
              
            """
            'L-BFGS-B' work long
            'Nelder-Mead' work somehow
            least_squares - worked excellent
            """

        def __search_index(self,v,x):
            return bisect.bisect(v, x) - 1 

        def set_FEM_data(self,Ne,dl,Ldl):
            self.Ne = Ne
            self.dl = dl
            self.Ldl = Ldl
            
        def Ldivide(self,step_mult=1,disp=False):
            """
            Discretize beam length on piecies with some step
            # Parameters
            ----------
            phi: sympy.Function
                Sympy Function Instance with the l parameter,
            disp: bool, optional
                Display data
            """
            self.step = self.L*step_mult*1e-3
            self.l_all_true = np.arange(0,self.L+self.step/2,self.step)
            self.N = len(self.l_all_true)
            if disp:
                display(Math("\\text{number of steps in simulation=}"+np.str_(self.N)))   

        def create_a(self,disp=False):
            """
            Creating a vector with phi,dphi,ddphi values based on provided phi function FEM formulation
            # Parameters
            ----------
            disp: bool, optional
                Display data
            """
            try:
                self.Ne
                self.dl
                self.Ldl
            except:
                raise ValueError("No FEM formulation was made!") from None
            try:
                self.fun_phi
                self.fun_dphi
                self.fun_ddphi
            except:
                raise ValueError("Call set_phi first. You should provide test phi function!") from None  
            self.a = np.array([])
            for m in range(6*self.Ne):
                self.a = np.append(self.a,self.__get_a(m)) 
            if disp:
                print("a = {}".format(self.a))       

        def set_a_diff(self,a_diff):
            self.a_diff = a_diff

        def set_a_approx(self,a_approx):
            self.a_approx = a_approx

        def get_a_diff(self):
            try:
                self.a_diff
            except:
                raise ValueError("Optimization wasn't. Don't have an approximation!") from None
            else:
                return self.a_diff

        def get_a_approx(self):
            try:
                self.a_approx
            except:
                raise ValueError("Optimization wasn't. Don't have an approximation!") from None
            else:
                return self.a_approx

        def __get_a(self,m):
            # there m from 0 to 6*Ne-1
            if m%6 == 0:
                return self.fun_phi(self.Ldl[(m)//6])
            if m%6 == 1:
                return self.fun_dphi(self.Ldl[(m)//6])
            if m%6 == 2:
                return self.fun_ddphi(self.Ldl[(m)//6])
            if m%6 == 3:
                return self.fun_phi(self.Ldl[(m)//6+1])
            if m%6 == 4:
                return self.fun_dphi(self.Ldl[(m)//6+1])
            if m%6 == 5:
                return self.fun_ddphi(self.Ldl[(m)//6+1])

        def __psi_choser(self,e,l):
            # e from 0 to Ne-1
            if l<self.Ldl[e+1] and l>=self.Ldl[e]:
                return 1
            else:
                return 0

        def __get_psi(self,l): # psi
            ret = np.array([])
            if l==self.L:
                ret = np.zeros((1,self.Ne*6-6))[0]
                for i in range(6):
                    ret=np.append(ret,np.polyval(self.p[(i)],(l-self.Ldl[self.Ne-1])/self.dl))
                return ret
            else:
                for e in range(self.Ne):
                    for i in range(6):
                        ret = np.append(ret,np.polyval(self.p[(i)],(l-self.Ldl[e])/self.dl))
                    ret[-6:]=ret[-6:]*self.__psi_choser(e,l)
                return ret  
        def __get_dpsi(self,l): # d^1psi/dlambda^1
            if l==self.L:
                return np.zeros((1,self.Ne*6))[0]
            else:
                ret = np.array([])
                for e in range(self.Ne):
                    for i in range(6):
                        ret = np.append(ret,np.polyval(self.dp[(i)],(l-self.Ldl[e])/self.dl))
                    ret[-6:]=ret[-6:]*self.__psi_choser(e,l)
                return ret
        def __get_ddpsi(self,l): # d^2psi/dlambda^2
            if l==self.L:
                return np.zeros((1,self.Ne*6))[0]
            else:
                ret = np.array([])
                for e in range(self.Ne):
                    for i in range(6):
                        ret = np.append(ret,np.polyval(self.ddp[(i)],(l-self.Ldl[e])/self.dl))
                    ret[-6:]=ret[-6:]*self.__psi_choser(e,l)
                return ret
        def __get_dddpsi(self,l): # d^3psi/dlambda^3
            if l==self.L:
                return np.zeros((1,self.Ne*6))[0]
            else:
                ret = np.array([])
                for e in range(self.Ne):
                    for i in range(6):
                        ret = np.append(ret,np.polyval(self.dddp[(i)],(l-self.Ldl[e])/self.dl))
                    ret[-6:]=ret[-6:]*self.__psi_choser(e,l)
                return ret
        def __get_ddddpsi(self,l): # d^4psi/dlambda^4
            if l==self.L:
                return np.zeros((1,self.Ne*6))[0]
            else:
                ret = np.array([])
                for e in range(self.Ne):
                    for i in range(6):
                        ret = np.append(ret,np.polyval(self.ddddp[(i)],(l-self.Ldl[e])/self.dl))
                    ret[-6:]=ret[-6:]*self.__psi_choser(e,l)
                return ret

        def show_one_element_approx(self,e=1):
            """
            Show one element phi's approximation together with base psi_i functions
            # Parameters
            ----------
            e: int, optional
                element number from 1 to Ne
            """
            try:
                self.a
            except:
                raise ValueError("Call create_a first!") from None
            for l in self.l_all_true:
                if not l:
                    psi1 = np.array([])
                    psi2 = np.array([])
                    psi3 = np.array([])
                    psi4 = np.array([])
                    psi5 = np.array([])
                    psi6 = np.array([])
                    psi = np.array([])
                if (l <= self.Ldl[e]) and (l >= self.Ldl[e-1]):
                    psi1 = np.append(psi1,self.a[0+6*(e-1)]*self.__get_psi(0,e-1,l))
                    psi2 = np.append(psi2,self.a[1+6*(e-1)]*self.__get_psi(1,e-1,l))
                    psi3 = np.append(psi3,self.a[2+6*(e-1)]*self.__get_psi(2,e-1,l))
                    psi4 = np.append(psi4,self.a[3+6*(e-1)]*self.__get_psi(3,e-1,l))
                    psi5 = np.append(psi5,self.a[4+6*(e-1)]*self.__get_psi(4,e-1,l))
                    psi6 = np.append(psi6,self.a[5+6*(e-1)]*self.__get_psi(5,e-1,l))
                    psi = np.append(psi,psi1[-1]+psi2[-1]+psi3[-1]+psi4[-1]+psi5[-1]+psi6[-1])
                else:
                    psi1 = np.append(psi1,0)
                    psi2 = np.append(psi2,0)
                    psi3 = np.append(psi3,0)
                    psi4 = np.append(psi4,0)
                    psi5 = np.append(psi5,0)
                    psi6 = np.append(psi6,0)
                    psi = np.append(psi,0)
            plt.subplots(figsize = (20,8))
            plt.ylabel("$\\varphi_{"+str(e)+"}(l,t=0)$ [deg]",fontsize=15)
            plt.xlabel("$l$")
            labels = ['$\psi_1$','$\psi_2$','$\psi_3$','$\psi_4$','$\psi_5$','$\psi_6$','$\psi$']
            colours = ['b','k','c','g','m','y','r']
            plt.title(str(e)+" element")
            plt.plot(self.l_all_true,np.rad2deg(psi1),label=labels[0],color=colours[0])
            plt.plot(self.l_all_true,np.rad2deg(psi2),label=labels[1],color=colours[1])
            plt.plot(self.l_all_true,np.rad2deg(psi3),label=labels[2],color=colours[2])
            plt.plot(self.l_all_true,np.rad2deg(psi4),label=labels[3],color=colours[3])
            plt.plot(self.l_all_true,np.rad2deg(psi5),label=labels[4],color=colours[4])
            plt.plot(self.l_all_true,np.rad2deg(psi6),label=labels[5],color=colours[5])
            plt.plot(self.l_all_true,np.rad2deg(psi),label=labels[6],color=colours[6])
            plt.legend(fontsize="15",loc='upper right',ncol=3)
            plt.grid(True)
            plt.show()

        def set_phi(self,phi: sm.core.expr.Expr,disp=False):
            """
            Creating phi,dphi,ddphi evaluation functions for given phi  test sympy.Function.
            # Parameters
            ---------- 
            phi: sympy.Function
                Sympy Function Instance with l parameter [in rad],
            disp: bool, optional
                Display data
            """
            try:
                self.step
                self.l_all_true
                self.N
            except:
                raise ValueError("Call Ldivide first!") from None
            l = sm.symbols("l")
            self.fun_phi = sm.lambdify(l, phi, modules='numpy')
            self.phi_true = self.fun_phi(self.l_all_true)
            dphi = sm.diff(phi,l)
            self.fun_dphi = sm.lambdify(l, dphi, modules='numpy')
            self.dphi_true = self.fun_dphi(self.l_all_true)
            ddphi = sm.diff(dphi,l)
            self.fun_ddphi = sm.lambdify(l, ddphi, modules='numpy')
            self.ddphi_true = self.fun_ddphi(self.l_all_true)

            self.x_phi_true = np.array([0])
            self.y_phi_true = np.array([0])
            for i in range(len(self.l_all_true)-1):
                self.x_phi_true = np.append(self.x_phi_true,
                    np.cos(self.fun_phi(self.l_all_true[i+1]))*self.step + self.x_phi_true[-1] )
                self.y_phi_true = np.append(self.y_phi_true,
                    np.sin(self.fun_phi(self.l_all_true[i+1]))*self.step + self.y_phi_true[-1])
                        
            if disp:
                display(Math("\\varphi=\Large"+sm.latex(phi)))    
                display(Math("\\frac{\\partial \\varphi}{\\partial l}=\Large"+sm.latex(dphi)))    
                display(Math("\\frac{\\partial^2 \\varphi}{\\partial l^2}=\Large"+sm.latex(ddphi))) 

                plt.subplots(1,5,figsize = (20,6))
                plt.subplot(141)
                plt.plot(self.x_phi_true,self.y_phi_true)
                plt.axis('equal')
                plt.title("given beam shape (x,y)",fontsize=15)
                plt.xlabel("$x$ [mm]",fontsize=15)
                plt.ylabel("$y$ [mm]",fontsize=15)
                plt.grid(True)

                plt.subplot(142)
                plt.plot(self.l_all_true,np.rad2deg(self.phi_true))
                plt.title("given beam shape (phi)",fontsize=15)
                plt.xlabel("$l$ [mm]",fontsize=15)
                plt.ylabel("$\\varphi(l,t=0)$ [deg]",fontsize=15)
                plt.grid(True)

                plt.subplot(143)
                plt.plot(self.l_all_true,np.rad2deg(self.dphi_true))
                plt.title("given beam shape (dphi)",fontsize=15)
                plt.xlabel("$l$ [mm]",fontsize=15)
                plt.ylabel("$\\frac{\partial\\varphi(l,t=0)}{\partial l}$ [deg/mm]",fontsize=15)
                plt.grid(True)

                plt.subplot(144)
                plt.plot(self.l_all_true,np.rad2deg(self.ddphi_true))
                plt.title("given beam shape (ddphi)",fontsize=15)
                plt.xlabel("$l$ [mm]",fontsize=15)
                plt.ylabel("$\\frac{\partial^2\\varphi(l,t=0)}{\partial l^2}$ [deg/mm^2]",fontsize=15)
                plt.grid(True)
                plt.show()

        def phi_approx_preparing(self):
            # preparing for fast computation next
            try:
                self.Ne
                self.dl
                self.Ldl
            except:
                raise ValueError("No FEM formulation was made!") from None
            try:
                self.step
                self.l_all_true
                self.N
            except:
                raise ValueError("Call Ldivide first!") from None
            
            flag_preparing_already_done = 0
            if os.path.isfile('psi_vectors.npz'):
                flag_preparing_already_done = 1

            if flag_preparing_already_done:
                print("Found numpy zip archive with preparing data: psi vectors. Checking if we can use it!")
                with np.load('psi_vectors.npz') as npzfile: # for closign after using it
                    self.psi = npzfile['psi']
                    self.dpsi = npzfile['dpsi']
                    self.ddpsi = npzfile['ddpsi']
                    N = npzfile['N']
                    Ne = npzfile['Ne']
                    dl = npzfile['dl']
                    step = npzfile['step']
                del npzfile
            
            if (not flag_preparing_already_done) or (not N==self.N) or (not Ne==self.Ne) or (not dl==self.dl) or (not step==self.step):
                if flag_preparing_already_done:
                    print("Checking finished. We cannot use this data as FEM or/and Ldivide parameters mismatch. Creating new one:")
                self.c1 = self.E*self.I/(self.rho*self.A)
                self.EI = self.E*self.I
                start_time = time.time_ns()
                time.sleep(0.000001) # sleep 1 us
                self.psi = np.zeros((self.N,6*self.Ne))
                self.dpsi = np.zeros((self.N,6*self.Ne))
                self.ddpsi = np.zeros((self.N,6*self.Ne))
                for (l,i) in zip(self.l_all_true,range(self.N)):  
                    self.psi[i] = self.__get_psi(l)
                    self.dpsi[i] = self.__get_dpsi(l)
                    self.ddpsi[i] =self.__get_ddpsi(l)
                time_end = time.time_ns()-start_time-1*1e3
                print("Preparing time: %s s" % (round(time_end*1e-9,3)))

                np.savez('psi_vectors.npz',psi=self.psi,dpsi=self.dpsi,\
                     ddpsi=self.ddpsi,N=self.N,Ne=self.Ne,step=self.step,dl=self.dl)
            else:
                if flag_preparing_already_done:
                    print("Checking finished. Using loaded data")
                

        def phi_approx(self,disp_time=True,der_num=0):
            """
            Phi function approximation. Opportunity to approximate dphi and ddphi functions.
            # Parameters
            ---------- 
            disp_time: bool, optional
                Display evaluation time of for cycle and one iteration 
            der_num: int, optional
                Number of phi derivative need to compute (from 0 to 2)
            # Note: 
            ----------
                - approximation work only when L>>1. It's very important
                - for convergence to only phi without dphi and ddphi need 6 FE for anough quality and 20 FE for excellent quality
                - for convergence to phi,dphi,ddphi need 20 FE for excellent quality
            """
            try:
                self.a_approx
            except:
                print("Optimization wasn't. Don't have an approximation!We will use a created by create_a fun!")
                flag_a_approx_is = 0
                try:
                    self.a
                except:
                    raise ValueError("Call create_a first!") from None
            else:
                print("Found an approximation. Will use it!")
                self.a = self.a_approx
                flag_a_approx_is = 1
                
            if der_num == 2:
                #evaluation
                start_time = time.time_ns()
                time.sleep(0.000001) # sleep 1 us
                phi_appr = np.matmul(self.psi,self.a)
                dphi_appr = np.matmul(self.dpsi,self.a)
                ddphi_appr = np.matmul(self.ddpsi,self.a)
                cos_phi_appr = np.cos(phi_appr)
                sin_phi_appr = np.sin(phi_appr)
                x = -self.step+np.cumsum(cos_phi_appr)*self.step
                y = np.cumsum(sin_phi_appr)*self.step
                end_time = time.time_ns()-start_time-1*1e3
                if end_time==0:
                    print("evaluation time is less then 1 ns")
                else:
                    print("evaluation time: %s ms" % (round(end_time*1e-6,3)))
                    print("time for 1 step: %s us" % (round(1e-3*end_time/self.N,3)))
            elif der_num == 1:
                #evaluation
                start_time = time.time_ns()
                time.sleep(0.000001) # sleep 1 us
                phi_appr = np.matmul(self.psi,self.a)
                dphi_appr = np.matmul(self.dpsi,self.a)
                cos_phi_appr = np.cos(phi_appr)
                sin_phi_appr = np.sin(phi_appr)
                x = -self.step+np.cumsum(cos_phi_appr)*self.step
                y = np.cumsum(sin_phi_appr)*self.step
                end_time = time.time_ns()-start_time-1*1e3
                if end_time==0:
                    print("evaluation time is less then 1 ns")
                else:
                    print("evaluation time: %s ms" % (round(end_time*1e-6,3)))
                    print("time for 1 step: %s us" % (round(1e-3*end_time/self.N,3)))
            elif der_num == 0:
                #evaluation
                start_time = time.time_ns()
                time.sleep(0.000001) # sleep 1 us
                phi_appr = np.matmul(self.psi,self.a)
                cos_phi_appr = np.cos(phi_appr)
                sin_phi_appr = np.sin(phi_appr)
                x = -self.step+np.cumsum(cos_phi_appr)*self.step
                y = np.cumsum(sin_phi_appr)*self.step
                end_time = time.time_ns()-start_time-1*1e3
                if end_time==0:
                    print("evaluation time is less then 1 ns")
                else:
                    print("evaluation time: %s ms" % (round(end_time*1e-6,3)))
                    print("time for 1 step: %s us" % (round(1e-3*end_time/self.N,3)))

            plt.subplots(2,2,figsize = (20,8))
            plt.subplot(221)
            labels = ['$\\varphi_{true}$','$\\varphi_{approx}$']
            colours = ['b','r']
            plt.plot(self.l_all_true/self.mult,phi_appr,label=labels[1],color=colours[1])
            if not flag_a_approx_is:
                plt.plot(self.l_all_true/self.mult,self.phi_true,"--",label=labels[0],color=colours[0])
            plt.grid(True)
            plt.xlabel("$l$ [mm]")
            plt.ylabel("$\\varphi(l,t=0)$")
            plt.legend(fontsize="15",loc='best')
            plt.title("phi approx and true")
            plt.subplot(222)
            labels = ['$(x,y)_{true}$','$(x,y)_{approx}$']
            colours = ['b','r']
            plt.plot(x/self.mult,y/self.mult,label=labels[1],color=colours[1])
            if not flag_a_approx_is:
                plt.plot(self.x_phi_true/self.mult,self.y_phi_true/self.mult,"--",label=labels[0],color=colours[0])
            plt.grid(True)
            plt.xlabel("$x$ [mm]")
            plt.ylabel("$y$ [mm]")
            plt.legend(fontsize="15",loc='best')
            plt.title("x,y approx and true")
            plt.subplot(223)
            labels = ['$\\frac{\partial \\varphi_{true} }{\partial l}$','$\\frac{\partial\\varphi_{approx}}{\partial l}$']
            colours = ['b','r']
            if der_num == 1 or der_num == 2:
                plt.plot(self.l_all_true,dphi_appr,label=labels[1],color=colours[1])
                if not flag_a_approx_is:
                    plt.plot(self.l_all_true,self.dphi_true,"--",label=labels[0],color=colours[0])
            plt.grid(True)
            plt.xlabel("$l$ [mm]")
            plt.ylabel("$\\frac{\partial\\varphi(l,t=0)}{\partial l}$")
            plt.legend(fontsize="15",loc='best')
            plt.title("dphi approx and true")
            plt.subplot(224)
            labels = ['$\\frac{\partial\\varphi^2_{true}}{\partial l^2}$','$\\frac{\partial\\varphi^2_{approx}}{\partial l^2}$']
            colours = ['b','r']
            if der_num == 2:
                plt.plot(self.l_all_true/self.mult,ddphi_appr,label=labels[1],color=colours[1])
                if not flag_a_approx_is:
                    plt.plot(self.l_all_true/self.mult,self.ddphi_true,"--",label=labels[0],color=colours[0])
            plt.grid(True)
            plt.xlabel("$l$ [mm]")
            plt.ylabel("$\\frac{\partial\\varphi^2(l,t=0)}{\partial l^2}$")
            plt.legend(fontsize="15",loc='best')
            plt.title("ddphi approx and true")
            plt.show()