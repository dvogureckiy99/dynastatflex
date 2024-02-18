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
Release 1.0.2
1. With very high Ne accurate angle. 
"""

class Flex_beam(object):
    def __init__(self,L=1,E=1,h=1,w=1,rho=1):
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
        self.rho = rho * 1e-9 / self.mult**3

    def Create_Simulation(self):
        """
        Creating Child Class Simulation instance and passing parameters to it 
        # Parameters
        ---------- 
        disp: bool, optional
            display data
        """
        self.Simulating = self.Simulating(self.L,self.E,self.h,self.w,self.rho,self.mult) # creating Child class instance and passing parameters to it

    def FEM(self,Ne=10,disp=False,polynom_deg=3):
        """
        Splitting a beam into finite elements
        # Parameters
        ---------- 
        Ne: int
            element's number,
        disp: bool, optional
            display data
        """
        self.polynom_deg = polynom_deg # polynomial degree
        self.Ne = Ne
        self.dl = self.L/Ne
        self.Ldl = np.linspace(0,self.L,self.Ne+1)
        try:
            self.Simulating.set_FEM_data(Ne,self.dl,self.Ldl,self.polynom_deg,disp) # inheritance of parameters by a child class
        except:
            raise ValueError("First call Create_Simulation method!") # from None
        if disp:
            display(Math("\\large \mathcal{N}_{0,1,\dots,Ne}=\\text{"+np.str_(self.Ldl)+" [mm]}"))   
            
    class Simulating(object):
        def __init__(self,L,E,h,w,rho,mult): # inheritance of parameters by a child class
            self.L = L
            self.E = E
            self.h = h
            self.w = w
            self.I = w**3*h/12
            self.EI = self.E*self.I
            self.rho = rho
            self.A = w*h
            self.mult = mult
            
        def set_FEM_data(self,Ne,dl,Ldl,polynom_deg,disp):
            self.Ne = Ne
            self.dl = dl
            self.Ldl = Ldl
            self.a_size = polynom_deg + 1  
            self.a_halfsize = int(self.a_size/2)
            if self.a_size==6:
                self.last_zeros = 2
            elif self.a_size==4:
                self.last_zeros = 1
            if self.a_size==6:
                C_val = np.array([[-6,15,-10,0,0,1],[-3,8,-6,0,1,0],[-0.5,1.5,-1.5,0.5,0,0],
                    [6,-15,10,0,0,0],[-3,7,-4,0,0,0],[0.5,-1,0.5,0,0,0]])
            elif self.a_size==4:
                C_val = np.array([[2,-3,0,1],[1,-2,1,0],[-2,3,0,0],[1,-1,0,0]])
            self.p = []
            self.dp = []
            self.ddp = []
            self.dddp = []
            if self.a_size==6:
                for i in range(self.a_size):
                    self.p.append(np.poly1d(C_val[i])) # for psi
                    self.dp.append(np.polyder(self.p[-1], m=1)) # for d^1psi/dlambda^1
                    self.ddp.append(np.polyder(self.p[-1], m=2)) # for d^2psi/dlambda^2
                    self.dddp.append(np.polyder(self.p[-1], m=3)) # for d^3psi/dlambda^3
            elif self.a_size==4:
                for i in range(self.a_size):
                    self.p.append(np.poly1d(C_val[i])) # for psi
                    self.dp.append(np.polyder(self.p[-1], m=1)) # for d^1psi/dlambda^1
                    self.ddp.append(np.polyder(self.p[-1], m=2)) # for d^2psi/dlambda^2
                    self.dddp.append(np.polyder(self.p[-1], m=3)) # for d^3psi/dlambda^3
            if disp:
                a = np.ones((1,self.a_size))[0]
                # a = np.array([0,1,12,1,1,13])
                x = np.arange(0,1+1e-3,1e-3)
                plt.subplots(2,3,figsize = (20,8))
                plt.subplot(231)
                for (obj,i) in zip(self.p,np.arange(self.a_size)):
                    y = np.polyval(obj,x)*a[i]
                    plt.plot(x,y)
                plt.grid(True)
                plt.title("$\psi_i$ - basis functions for 1 FE",fontsize=22)
                plt.subplot(232)
                for (obj,i) in zip(self.dp,np.arange(self.a_size)):
                    y = np.polyval(obj,x)*a[i]
                    plt.plot(x,y)
                plt.grid(True)
                plt.title("$\\frac{\partial \psi_i(\lambda)}{\partial \lambda}$ - basis functions for 1 FE",fontsize=22)
                plt.subplot(233)
                for (obj,i) in zip(self.ddp,np.arange(self.a_size)):
                    y = np.polyval(obj,x)*a[i]
                    plt.plot(x,y)
                plt.grid(True)
                plt.title("$\\frac{\partial^2 \psi_i(\lambda)}{\partial \lambda^2}$ - basis functions for 1 FE",fontsize=22)
                plt.subplot(234)
                y = np.zeros((1,len(x)))[0]
                for obj,i in zip(self.p,np.arange(self.a_size)):
                    y = y + np.polyval(obj,x)*a[i]
                plt.plot(x,y)
                plt.grid(True)
                plt.title("$\sum_i\psi_i$ - sum of basis functions for 1 FE",fontsize=22)
                plt.subplot(235)
                y = np.zeros((1,len(x)))[0]
                for obj,i in zip(self.dp,np.arange(self.a_size)):
                    y = y + np.polyval(obj,x)*a[i]
                plt.plot(x,y)
                plt.grid(True)
                plt.title("$\sum_i \\frac{\partial \psi_i(\lambda)}{\partial \lambda}$ - sum of basis functions for 1 FE",fontsize=22)
                plt.subplot(236)
                y = np.zeros((1,len(x)))[0]
                for obj,i in zip(self.ddp,np.arange(self.a_size)):
                    y = y + np.polyval(obj,x)*a[i]
                plt.plot(x,y)
                plt.grid(True)
                plt.title("$\sum_i \\frac{\partial^2 \psi_i(\lambda)}{\partial \lambda^2}$ - sum of basis functions for 1 FE",fontsize=22)
                plt.tight_layout()
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
            lines = np.array2string(a,formatter={'float':lambda x: "%.6f" % x}).replace('\n  ', ' ').replace('[', '').replace(']', '').splitlines()
            rv = [r'\begin{bmatrix}']
            rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
            rv +=  [r'\end{bmatrix}']
            return '\n'.join(rv)

        def __diag_mat(self,A,diag_num,flag_with_int=0):
            size1 = np.shape(A)[0]
            size2 = np.shape(A)[1]
            B = np.array([]).reshape((0,diag_num*size2))
            for d_e in range(diag_num):
                row = np.array([]).reshape((size1,0))
                for i in range(diag_num):
                    if i==d_e:
                        if not flag_with_int:
                            row = np.hstack((row, A))
                        else:
                            if not d_e:
                                row = np.hstack((row, A))
                                A_one = A
                                A_last = A
                            else:
                                A = A_last + A_one
                                row = np.hstack((row,A))
                                A_last = A
                    else:
                        row = np.hstack((row, np.zeros((size1,size2))))
                B = np.vstack((B,row))          
            return B
        def __diag_shift_mat(self,A,diag_num,flag_with_int=0):
            size1 = np.shape(A)[0]
            size2 = np.shape(A)[1]
            halfsize = int(size2/2)
            B = np.array([]).reshape((0,(diag_num+1)*halfsize))
            for d_e in range(diag_num):
                row = np.array([]).reshape((size1,0))
                if not flag_with_int:
                    row = np.hstack(( np.zeros((size1,d_e*halfsize)),row ))
                    row = np.hstack(( row, A ))
                    row = np.hstack((row, np.zeros((size1,halfsize*diag_num-(d_e+1)*halfsize)) ))
                else:
                    if not d_e:
                        row = np.hstack(( np.zeros((size1,d_e*halfsize)),row ))
                        row = np.hstack(( row, A ))
                        row = np.hstack((row, np.zeros((size1,halfsize*diag_num-(d_e+1)*halfsize)) ))
                        A_one = A
                        A_last = A
                    else:
                        A = A_last + A_one
                        row = np.hstack(( np.zeros((size1,d_e*halfsize)),row ))
                        row = np.hstack(( row, A ))
                        row = np.hstack((row, np.zeros((size1,halfsize*diag_num-(d_e+1)*halfsize)) ))
                        A_last = A
                B = np.vstack((B,row))          
            return B
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

        def __fun_static_optim(self,a_diff,disp):
            # preparing a vector for each FE, cause a_diff contain only unique values 
            a = np.array([0])
            a = np.append(a,a_diff)
            a = np.append(a,np.zeros((1,self.last_zeros))[0])
            
            self.iteration_num += 1   

            dphi_appr_power3 =  np.power(np.matmul(self.dpsi,a),3)  # [1,N]

            if self.flag_Fextxy:
                phi_appr = np.matmul(self.psi,a)  # [1,N]
                dphi_appr = np.matmul(self.dpsi,a)  # [1,N]
                ddphi_appr = np.matmul(self.ddpsi,a)  # [1,N]
                sinphiappr = np.sin(phi_appr)
                cosphiappr = np.cos(phi_appr)

                ddphi_appr_power2 =  np.power(ddphi_appr,2)  # [1,N]
                dphi_appr_power2 =  np.power(dphi_appr,2)  # [1,N]

                Fext_perp_int = -sinphiappr[-1]*self.Fxpsi[-1]+sinphiappr[0]*self.Fxpsi[0] +\
                        cosphiappr[-1]*self.Fypsi[-1]-cosphiappr[0]*self.Fypsi[0] +\
                        np.sum(np.multiply(sinphiappr.reshape(self.N_optim,1),self.Fxdpsi)*self.step_optim,axis=0) -\
                        np.sum(np.multiply(cosphiappr.reshape(self.N_optim,1),self.Fydpsi)*self.step_optim,axis=0)
                
                Fext_para_int = cosphiappr[-1]*self.Fxpsi[-1]-cosphiappr[0]*self.Fxpsi[0] +\
                        sinphiappr[-1]*self.Fypsi[-1]-sinphiappr[0]*self.Fypsi[0] -\
                        np.sum(np.multiply(cosphiappr.reshape(self.N_optim,1),self.Fxdpsi)*self.step_optim,axis=0) -\
                        np.sum(np.multiply(sinphiappr.reshape(self.N_optim,1),self.Fydpsi)*self.step_optim,axis=0)

                cost = np.concatenate([ Fext_para_int+self.EI*\
                    (np.sum(np.multiply(ddphi_appr_power2.reshape(self.N_optim,1),self.psi)*self.step_optim,axis=0)+\
                        2*(ddphi_appr[-1]*dphi_appr[-1]*self.psi[-1]-ddphi_appr[0]*dphi_appr[0]*self.psi[0])-\
                        dphi_appr_power2[-1]*self.dpsi[-1]+dphi_appr_power2[0]*self.dpsi[0]+\
                        np.sum(np.multiply(dphi_appr_power2.reshape(self.N_optim,1),self.ddpsi)*self.step_optim,axis=0)),\
                    Fext_perp_int-self.EI*(np.matmul(self.F,a)+\
                        (1/3)*(np.sum(np.multiply(dphi_appr_power3.reshape(self.N_optim,1),self.dpsi)*self.step_optim,axis=0)-\
                        dphi_appr_power3[int(self.N_optim-1)]*self.psi[int(self.N_optim-1)]+dphi_appr_power3[0]*self.psi[0])),\
                    [self.Fxext_Fx-self.EI*(np.sum(np.multiply(sinphiappr_ddphiappr,self.dpsi[:self.ind_N2,self.a_halfsize])*\
                                                        self.step_optim,axis=0)) ],\
                    [self.Fyext_Fy+self.EI*(np.sum(np.multiply(cosphiappr_ddphiappr,self.dpsi[:self.ind_N2,self.a_halfsize])*\
                                                    self.step_optim,axis=0))]   ])
            else:
                phi_appr = np.matmul(self.psi[:self.ind_N2,:6],a[:6])  # [1,N]
                ddphi_appr = np.matmul(self.ddpsi[:self.ind_N2,:6],a[:6])  # [1,N]
                sinphiappr = np.sin(phi_appr)
                cosphiappr = np.cos(phi_appr)
                sinphiappr_ddphiappr = np.multiply(sinphiappr,ddphi_appr)
                cosphiappr_ddphiappr = np.multiply(cosphiappr,ddphi_appr)

                Fextx = np.sum(np.multiply( sinphiappr,self.Fext)*\
                                           self.step_optim,axis=0)
                Fexty = np.sum(np.multiply( cosphiappr,self.Fext)*\
                                           self.step_optim,axis=0)

                cost = np.concatenate([ self.Fext_int-self.EI*(np.matmul(self.F,a)+\
                        (1/3)*(np.sum(np.multiply(dphi_appr_power3.reshape(self.N_optim,1),self.dpsi)*self.step_optim,axis=0)-\
                    dphi_appr_power3[int(self.N_optim-1)]*self.psi[int(self.N_optim-1)]+dphi_appr_power3[0]*self.psi[0])),\
                        [-Fextx-self.EI*(np.sum(np.multiply(sinphiappr_ddphiappr,self.dpsi[:self.ind_N2,self.a_halfsize])*\
                                                        self.step_optim,axis=0)) ],\
                        [Fexty+self.EI*(np.sum(np.multiply(cosphiappr_ddphiappr,self.dpsi[:self.ind_N2,self.a_halfsize])*\
                                                    self.step_optim,axis=0))]   ])
            # cost = np.sum(np.power(cost,2))
            self.phi_end = np.matmul(self.psi,a)[-1]
            if self.optim_alg == 'Nelder-Mead':
                cost = np.sum(np.power(cost,2))
                if disp:
                    print("iter={},cost={}".format(self.iteration_num,cost))
            elif self.optim_alg == 'least_squares':
                if disp:
                    print("iter={},cost={}".format(self.iteration_num,np.sum(np.power(cost,2))))
            # print("iter={}".format(self.iteration_num))
            return cost
            
        def __delta1(self,l):
            if l<0:
                return 0
            else:
                return 1

        def __delta_approx(self,l,dl,e):
            # https://math.stackexchange.com/questions/4280517/simplest-smooth-c-infty-approximation-to-diracs-delta-with-bounded-su
            # https://mathworld.wolfram.com/DeltaFunction.html
            # return e/((l-dl)**2+e**2)/np.pi
            return (1-np.tanh((l-dl)/e)**2)/2/e

        def static_preparing(self,disp=True,Fext_in=1,l_Fext=None,Fext_type='triangle',widthofFextindl=1):
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
            self.Fext_in = Fext_in
            self.l_Fext = l_Fext
            self.Fext_type = Fext_type

            self.c1 = self.E*self.I/(self.rho*self.A)
            self.c3 = 1/(self.rho*self.A)
            
            if disp:
                start_time = time.time_ns()
                time.sleep(0.000001) # sleep 1 us
            # preparing for fast computation next
            
            # self.F = self.__diag_shift_mat(self.F,self.Ne,1)
            # print(np.shape(self.F))
            self.M = np.zeros((self.a_size,self.a_size))
            # for j in range(self.a_size):
            #     for i in range(self.a_size):
            #         self.M[j][i] = sp.integrate.quad(self.__M_int,0,self.Ldl[1],args=(i,j))[0]
            # self.M = self.__diag_shift_mat(self.M,self.Ne,1)

            if disp:
                time_end = time.time_ns()-start_time-1*1e3
                print("Preparing time: %s s" % (round(time_end*1e-9,3)))
            
            self.ind_N2 = self.__search_index(self.l_all_optim,self.Ldl[2])+1

            self.psi = self.__get_psi(self.step_optim)
            self.dpsi = self.__get_dpsi(self.step_optim)
            self.ddpsi = self.__get_ddpsi(self.step_optim)
            self.dddpsi = self.__get_dddpsi(self.step_optim)
            self.psi = self.__diag_shift_mat(self.psi,self.Ne)
            self.dpsi = self.__diag_shift_mat(self.dpsi,self.Ne)
            self.ddpsi = self.__diag_shift_mat(self.ddpsi,self.Ne)
            self.dddpsi = self.__diag_shift_mat(self.dddpsi,self.Ne)
            self.index = np.array([])
            for i in range(self.Ne-1):
                self.index = np.append(self.index,(self.steps_per_fe4optim+1)+(self.steps_per_fe4optim+1)*i) 
            self.index = np.int16(self.index)
            self.psi = np.delete(self.psi, self.index,axis=0)
            self.dpsi = np.delete(self.dpsi, self.index,axis=0)
            self.ddpsi = np.delete(self.ddpsi, self.index,axis=0)
            self.dddpsi = np.delete(self.dddpsi, self.index,axis=0)
            # boundary cond
            self.ddpsi[-1] = np.zeros((1,self.a_halfsize*(self.Ne+1)))[0]

            self.F = np.zeros((self.a_halfsize*(self.Ne+1),self.a_halfsize*(self.Ne+1)))
            for j in range(self.a_halfsize*(self.Ne+1)):
                for i in range(self.a_halfsize*(self.Ne+1)):
                    self.F[j][i] = np.sum(self.ddpsi[:,i]*self.ddpsi[:,j]*self.step_optim,axis=0) +\
                        self.dddpsi[-1,i]*self.psi[-1,j]-self.dddpsi[0,i]*self.psi[0,j]-\
                        self.ddpsi[-1,i]*self.dpsi[-1,j]+self.ddpsi[0,i]*self.dpsi[0,j]

            # preparing ddFext
            if l_Fext==None:
                l_Fext = self.L/2
            else:
                l_Fext = l_Fext * 1e3 * self.mult # point of application of force
            self.flag_Fextxy = 0
            if Fext_type=='delta':
                force_appl_point = self.__search_index(self.l_all_optim,l_Fext)
                if np.shape(self.Fext_in):
                    self.flag_Fextxy = 1
                    Fx = self.Fext_in[0]
                    Fy = self.Fext_in[1]
                    Fxext = np.zeros((1,self.N_optim))[0]   
                    Fyext = np.zeros((1,self.N_optim))[0]   
                    Fxext[int(force_appl_point)]=Fx
                    Fyext[int(force_appl_point)]=Fy
                    for p in range(int(self.steps_per_fe4optim*widthofFextindl)-1):
                        Fxext[int(force_appl_point)+p+1]=Fx*(1-(p+1)/self.steps_per_fe4optim/widthofFextindl)
                        Fxext[int(force_appl_point)-p-1]=Fx*(1-(p+1)/self.steps_per_fe4optim/widthofFextindl)
                        Fyext[int(force_appl_point)+p+1]=Fy*(1-(p+1)/self.steps_per_fe4optim/widthofFextindl)
                        Fyext[int(force_appl_point)-p-1]=Fy*(1-(p+1)/self.steps_per_fe4optim/widthofFextindl)
                    self.Fxext_int = np.sum(np.cumsum( np.multiply( Fxext[:self.ind_N2],self.psi[:self.ind_N2,self.a_halfsize]) *\
                                                self.step_optim,axis=0)*self.step_optim,axis=0)
                    self.Fyext_int = np.sum(np.cumsum( np.multiply( Fyext[:self.ind_N2],self.psi[:self.ind_N2,self.a_halfsize]) *\
                                                  self.step_optim,axis=0)*self.step_optim,axis=0)
                    self.Fxext = np.sum(np.multiply( Fxext.reshape(self.N_optim,1),self.psi)*self.step_optim,axis=0)
                    self.Fyext = np.sum(np.multiply( Fyext.reshape(self.N_optim,1),self.psi)*self.step_optim,axis=0)
                else:
                    Fext = np.zeros((1,self.N_optim))[0]   
                    # dw = w/(self.step_optim*self.steps_per_fe4optim)
                    Fext[int(force_appl_point)]=self.Fext_in
                    for p in range(int(self.steps_per_fe4optim*widthofFextindl)-1):
                        Fext[int(force_appl_point)+p+1]=self.Fext_in*(1-(p+1)/self.steps_per_fe4optim/widthofFextindl)
                        Fext[int(force_appl_point)-p-1]=self.Fext_in*(1-(p+1)/self.steps_per_fe4optim/widthofFextindl)
                    
                    self.Fext = np.multiply( Fext[:self.ind_N2],self.psi[:self.ind_N2,self.a_halfsize])
                    self.Fext_int = np.sum(np.multiply( Fext.reshape(self.N_optim,1),self.dpsi)*self.step_optim,axis=0)
                
                dFext = np.zeros((1,self.N_optim))[0] 
                # dFext[int(force_appl_point)]=dw
                # for p in range(self.steps_per_fe4optim-1):
                #     dFext[int(force_appl_point)+p+1]=dw*(1-(p+1)/self.steps_per_fe4optim)
                #     dFext[int(force_appl_point)-p-1]=dw*(1-(p+1)/self.steps_per_fe4optim)
                # self.dFext = np.sum(np.multiply( dFext.reshape(self.N_optim,1),self.psi)*self.step_optim,axis=0) 

                if disp:
                    # print("distributed integral error =%e"%(np.sum(Fext*self.step_optim*self.steps_per_fe4optim)-Fext_max))
                    plt.figure(figsize = (20,4))
                    plt.subplot(1,2,1)
                    if self.flag_Fextxy:
                        plt.plot(self.l_all_optim,Fxext)
                    else:
                        plt.plot(self.l_all_optim,Fext)
                    plt.plot(self.Ldl,np.zeros((1,self.Ne+1))[0],"og")
                    plt.grid()
                    if self.flag_Fextxy:
                        plt.title("Fxext - distributed force [N/m]")
                    else:
                        plt.title("Fext - distributed force [N/m]")
                    plt.subplot(1,2,2)
                    if self.flag_Fextxy:
                        plt.plot(self.l_all_optim,Fyext)
                    else:
                        plt.plot(self.l_all_optim,dFext)
                    plt.plot(self.Ldl,np.zeros((1,self.Ne+1))[0],"og")
                    plt.grid()
                    if self.flag_Fextxy:
                        plt.title("Fyext - distributed force [N/m]")
                    else:
                        plt.title("dFext - distributed force der [N/m^2]")
                    plt.show()
                    display(Math("\\bm{F}="+self.__bmatrix(self.F[0:self.a_size,0:self.a_size])))
                    display(Math("\\bm{M}="+self.__bmatrix(self.M[0:self.a_size,0:self.a_size])))
                    # display(Math("\\bm{F}_{ext}^{'}="+self.__bmatrix(self.dFext)))
            elif Fext_type=='const':
                force_appl_point = self.__search_index(self.l_all_optim,l_Fext)
                if np.shape(self.Fext_in):
                    self.flag_Fextxy = 1
                    Fxext = np.zeros((1,self.N_optim))[0]   
                    Fyext = np.zeros((1,self.N_optim))[0]   
                    if widthofFextindl==-1:
                        Fx = self.Fext_in[0]
                        Fy = self.Fext_in[1]
                        Fxext[0:int(force_appl_point)+1]=Fx
                        Fyext[0:int(force_appl_point)+1]=Fy
                    else:
                        Fx = self.Fext_in[0]/(widthofFextindl*2*self.Ldl[1])
                        Fy = self.Fext_in[1]/(widthofFextindl*2*self.Ldl[1])
                        Fxext[int(force_appl_point)-int(self.steps_per_fe4optim*widthofFextindl):\
                            int(force_appl_point)+int(self.steps_per_fe4optim*widthofFextindl)+1]=Fx
                        Fyext[int(force_appl_point)-int(self.steps_per_fe4optim*widthofFextindl):\
                            int(force_appl_point)+int(self.steps_per_fe4optim*widthofFextindl)+1]=Fy
                    self.Fext_para = np.multiply( Fxext.reshape(self.N_optim,1),self.dpsi)
                    self.Fxpsi = np.multiply( Fxext.reshape(self.N_optim,1),self.psi)
                    self.Fypsi = np.multiply( Fyext.reshape(self.N_optim,1),self.psi)
                    self.Fxdpsi = np.multiply( Fxext.reshape(self.N_optim,1),self.dpsi)
                    self.Fydpsi = np.multiply( Fyext.reshape(self.N_optim,1),self.dpsi)
                    self.Fxext_fx = np.sum(np.multiply( Fxext[:self.ind_N2],self.psi[:self.ind_N2,self.a_halfsize])*\
                                                  self.step_optim,axis=0)
                    self.Fyext_fy = np.sum(np.multiply( Fyext[:self.ind_N2],self.psi[:self.ind_N2,self.a_halfsize])*\
                                                  self.step_optim,axis=0)
                else:
                    # dw1 = 2*w/(self.step_optim)
                    # dw2 = w/(self.step_optim)
                    w = Fext_in # force at some point
                    Fext = np.zeros((1,self.N_optim))[0] 
                    if widthofFextindl==-1:
                        Fext[0:int(force_appl_point)+1]=w
                    elif widthofFextindl==-2:
                        Fext[int(force_appl_point):]=w
                    elif widthofFextindl==-3:
                        Fext[0:int(force_appl_point)+1]=w
                        Fext[int(force_appl_point):]=-w
                    else:
                        Fext[int(force_appl_point)-int(self.steps_per_fe4optim*widthofFextindl):\
                            int(force_appl_point)+int(self.steps_per_fe4optim*widthofFextindl)+1]=w 
                    self.Fext = np.multiply( Fext[:self.ind_N2],self.psi[:self.ind_N2,self.a_halfsize])
                    self.Fext_int = -np.sum(np.multiply( Fext.reshape(self.N_optim,1),self.dpsi)*self.step_optim,axis=0)
                    dFext = np.zeros((1,self.N_optim))[0]
                    # dFext[0]=dw1 
                    # dFext[int(force_appl_point)]=-dw2
                    # self.dFext = np.sum(np.multiply( dFext.reshape(self.N_optim,1),self.psi)*self.step_optim,axis=0) 

                if disp:
                    # print("distributed integral error =%e"%(np.sum(Fext*self.step_optim*self.steps_per_fe4optim)-Fext_max))
                    plt.figure(figsize = (20,4))
                    plt.subplot(1,2,1)
                    if self.flag_Fextxy:
                        plt.plot(self.l_all_optim,Fxext)
                    else:
                        plt.plot(self.l_all_optim,Fext)
                    plt.plot(self.Ldl,np.zeros((1,self.Ne+1))[0],"og")
                    plt.grid()
                    if self.flag_Fextxy:
                        plt.title("Fxext - distributed force [N/m]")
                    else:
                        plt.title("Fext - distributed force [N/m]")
                    plt.subplot(1,2,2)
                    if self.flag_Fextxy:
                        plt.plot(self.l_all_optim,Fyext)
                    else:
                        plt.plot(self.l_all_optim,dFext)
                    plt.plot(self.Ldl,np.zeros((1,self.Ne+1))[0],"og")
                    plt.grid()
                    if self.flag_Fextxy:
                        plt.title("Fyext - distributed force [N/m]")
                    else:
                        plt.title("dFext - distributed force der [N/m^2]")
                    plt.show()
                    display(Math("\\bm{F}="+self.__bmatrix(self.F[0:self.a_size,0:self.a_size])))
                    display(Math("\\bm{M}="+self.__bmatrix(self.M[0:self.a_size,0:self.a_size])))
                    # display(Math("\\bm{F}_{ext}^{'}="+self.__bmatrix(self.dFext)))
            elif Fext_type=='triangle':
                dw = 2*self.Fext_in/self.L
                # force_appl_point = self.__search_index(self.l_all_optim,l_Fext)

                Fext = np.zeros((1,self.N_optim))[0]
                dFext = np.zeros((1,self.N_optim))[0] 
                # for (l,i) in zip(self.l_all_optim,range(self.N_optim)):
                #     Fext[i]=dw*(l-2*self.__delta1(l-l_Fext)*(l-l_Fext))
                    # dFext[i]=dw-2*self.__delta1(l-l_Fext)*dw
                for (l,i) in zip(self.l_all_optim,range(self.N_optim)):
                    Fext[i]=dw*l-self.__delta1(l-l_Fext)*(2*self.Fext_in)
                self.Fext = np.multiply( Fext[:self.ind_N2],self.psi[:self.ind_N2,self.a_halfsize])
                self.Fext_int = -np.sum(np.multiply( Fext.reshape(self.N_optim,1),self.dpsi)*self.step_optim,axis=0) 
                # self.dFext = np.sum(np.multiply( dFext.reshape(self.N_optim,1),self.psi)*self.step_optim,axis=0) 

                if disp:
                    # print("distributed integral integral error =%e"%(np.sum(Fext*self.step)-Fext_max))
                    plt.figure(figsize = (20,4))
                    plt.subplot(1,2,1)
                    plt.title("Fext - distributed force [N/m]")
                    plt.plot(self.l_all_optim,Fext)
                    plt.plot(self.Ldl,np.zeros((1,self.Ne+1))[0],"og")
                    plt.grid()
                    plt.subplot(1,2,2)
                    plt.title("dFext - distributed force derivative [N/m^2]")
                    plt.plot(self.l_all_optim,dFext)
                    plt.plot(self.Ldl,np.zeros((1,self.Ne+1))[0],"og")
                    plt.grid()
                    plt.show()
                    display(Math("\\bm{F}="+self.__bmatrix(self.F[0:self.a_size,0:self.a_size])))
                    display(Math("\\bm{M}="+self.__bmatrix(self.M[0:self.a_size,0:self.a_size])))

        def static(self,disp=True,a0=[1,2],flag_compute_a_anyway=1,optim_alg=0):
            flag_preparing_already_done = 0
            if os.path.isfile('a.npz'):
                flag_preparing_already_done = 1

            if flag_preparing_already_done:
                if disp:
                    print("Found numpy zip archive with a approx data. Checking if we can use it!")
                with np.load('a.npz') as npzfile: # for closign after using it
                    self.a_approx = npzfile['a']
                    Fext_point = npzfile['Fext_point']
                    l_Fext = npzfile['l_Fext']
                    c1 = npzfile['c1']
                    c3 = npzfile['c3']
                    EI = npzfile['EI']
                    N = npzfile['N']
                    Ne = npzfile['Ne']
                    dl = npzfile['dl']
                    step = npzfile['step']
                    Fext_type = npzfile['Fext_type']
                del npzfile

            if (not flag_preparing_already_done) or (not N==self.N) or (not Ne==self.Ne) or (not dl==self.dl)\
                  or (not step==self.step) or (not c1==self.c1) or (not c3==self.c3) or (not EI==self.EI)\
                     or (not Fext_point[1]==self.Fext_in[1]) or (not l_Fext==self.l_Fext) or\
                          (not Fext_type==self.Fext_type) or flag_compute_a_anyway:
                if flag_preparing_already_done:
                    if disp:
                        print("Checking finished. We cannot use this a approx data as some parameters mismatch. Starting optimization:")
                else:
                    if disp:
                        print("Starting optimization:")

                self.iteration_num = 0
                if np.shape(a0)[0]<3:
                    self.a_diff = np.ones((1,self.a_halfsize*(self.Ne+1)-self.last_zeros-1))[0]
                else:
                    self.a_diff = a0
                """
                bound_min = np.zeros((1,len(a0)))[0]
                bound_max = np.zeros((1,len(a0)))[0]
                bound_min[0]=-5
                bound_max[0]=5
                bound_min[1]=-10
                bound_max[1]=10
                for i in range(len(a0)-2-1):
                    if i%3 == 0:
                        bound_min[i+2]=-np.pi
                        bound_max[i+2]=np.pi
                    if i%3 == 1:
                        bound_min[i+2]=-5
                        bound_max[i+2]=5
                    if i%3 == 2:
                        bound_min[i+2]=-10
                        bound_max[i+2]=10
                bound_min[-1]=-np.pi
                bound_max[-1]=np.pi
                """
                if disp:
                    start_time = time.time()
                if optim_alg:
                    self.optim_alg = 'Nelder-Mead'
                    res = sp.optimize.minimize(self.__fun_static_optim,self.a_diff,method='Nelder-Mead',args=(disp,),\
                                           options={'maxiter':int(1e6)})
                elif optim_alg==0:
                    self.optim_alg = 'least_squares'
                    tol=1e-3
                    res = sp.optimize.least_squares(self.__fun_static_optim,self.a_diff,\
                                                    ftol=tol,gtol=tol,xtol=tol,max_nfev=1e6,method='trf',\
                                                        args=(disp,))
                if disp:
                    end_time = time.time()-start_time
                    print("status: %s"%(res.message))
                    print("status: %s"%(res.status))
                    print("evaluation time:%s ms" % (round(1e3*end_time,3)))
                    print("time on 1 iter:%s ms" % (round(1e3*end_time/self.iteration_num,0)))  
                    print("iteration number:%s" % (self.iteration_num))
                    if self.optim_alg == 'least_squares':
                        print("res cost = {}".format(res.cost))  

                for i in range(self.a_halfsize*(self.Ne+1)-self.last_zeros-1):
                    self.a_diff[i] = res.x[i]

                self.a_approx = np.array([0])
                self.a_approx = np.append(self.a_approx,self.a_diff)
                self.a_approx = np.append(self.a_approx,np.zeros((1,self.last_zeros))[0])
                
                """
                'L-BFGS-B' work long
                'Nelder-Mead' work somehow
                least_squares - worked excellent
                """
                np.savez('a.npz',\
                        c1=self.c1,EI=self.EI,c3=self.c3,\
                        N=self.N,Ne=self.Ne,step=self.step,\
                        dl=self.dl,a=self.a_approx,Fext_point=self.Fext_in,\
                        l_Fext=self.l_Fext,Fext_type=self.Fext_type)
            else:
                if flag_preparing_already_done:
                    if disp:
                        print("Checking finished. Using loaded a_approx data!")

        def __search_index(self,v,x):
            return bisect.bisect(v, x) - 1 
            
        def Ldivide(self,steps_per_fe=1,steps_per_fe4optim=1,disp=False):
            """
            Discretize beam length on piecies with some step
            # Parameters
            ----------
            phi: sympy.Function
                Sympy Function Instance with the l parameter,
            disp: bool, optional
                Display data
            """
            self.steps_per_fe = steps_per_fe
            self.l_all_true = np.linspace(0,self.L,self.Ne*steps_per_fe+1)
            self.step = self.l_all_true[1] - self.l_all_true[0]
            self.N = len(self.l_all_true)

            self.steps_per_fe4optim = steps_per_fe4optim
            self.l_all_optim = np.linspace(0,self.L,self.Ne*steps_per_fe4optim+1)
            self.N_optim = len(self.l_all_optim)
            self.step_optim = self.l_all_optim[1] - self.l_all_optim[0]

            self.l_all_dl = np.linspace(0,self.L,self.Ne+1)
            self.N_dl = self.Ne+1
            self.step_dl = self.l_all_dl[1] - self.l_all_dl[0]

            if disp:
                display(Math("\\text{number of steps in optimization=}"+np.str_(self.N_optim)))  
                display(Math("\\text{number of steps in graphic constructino=}"+np.str_(self.N)))  
                display(Math("\\text{number of points in FEM=}"+np.str_(self.Ne+1)))   

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
            for m in range(self.a_halfsize*(self.Ne+1)):
                self.a = np.append(self.a,self.__get_a(m)) 
            if disp:
                display(Math("a="+self.__bmatrix(self.a))) 

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
            if self.a_size==6:
                if m%self.a_halfsize == 0:
                    return self.fun_phi(self.Ldl[(m)//self.a_halfsize])
                if m%self.a_halfsize == 1:
                    return self.fun_dphi(self.Ldl[(m)//self.a_halfsize])
                if m%self.a_halfsize == 2:
                    return self.fun_ddphi(self.Ldl[(m)//self.a_halfsize+1])
            elif self.a_size==4:
                if m%self.a_halfsize == 0:
                    return self.fun_phi(self.Ldl[(m)//self.a_halfsize])
                if m%self.a_halfsize == 1:
                    return self.fun_dphi(self.Ldl[(m)//self.a_halfsize])

        def __psi_choser(self,e,l):
            # e from 0 to Ne-1
            if l<self.Ldl[e+1] and l>=self.Ldl[e]:
                return 1
            else:
                return 0
        def __get_psi(self,step): # psi  
            ret = np.array([]).reshape((0,self.a_size))
            L = np.arange(0,self.Ldl[1]+step/2,step)
            L /= self.Ldl[1]
            for l in L:
                l_line = np.array([])
                for i in range(self.a_size):
                    l_line = np.append(l_line,np.polyval(self.p[(i)],l))
                ret = np.vstack((ret, l_line))
            return ret
        def __get_dpsi(self,step): # dpsi
            ret = np.array([]).reshape((0,self.a_size))
            L = np.arange(0,self.Ldl[1]+step/2,step)
            L /= self.Ldl[1]
            for l in L:
                l_line = np.array([])
                for i in range(self.a_size):
                    l_line = np.append(l_line,np.polyval(self.dp[(i)],l))
                ret = np.vstack((ret, l_line))
            return ret
        def __get_ddpsi(self,step): # ddpsi
            ret = np.array([]).reshape((0,self.a_size))
            L = np.arange(0,self.Ldl[1]+step/2,step)
            L /= self.Ldl[1]
            for l in L:
                l_line = np.array([])
                for i in range(self.a_size):
                    l_line = np.append(l_line,np.polyval(self.ddp[(i)],l))
                ret = np.vstack((ret, l_line))
            return ret
        def __get_dddpsi(self,step): # dddpsi
            ret = np.array([]).reshape((0,self.a_size))
            L = np.arange(0,self.Ldl[1]+step/2,step)
            L /= self.Ldl[1]
            for l in L:
                l_line = np.array([])
                for i in range(self.a_size):
                    l_line = np.append(l_line,np.polyval(self.dddp[(i)],l))
                ret = np.vstack((ret, l_line))
            return ret

        def __M_int(self,l,i,j):
            return np.polyval(self.p[(i)],l/self.Ldl[1])*np.polyval(self.p[(j)],l/self.Ldl[1])

        def __F_int(self,l,i,j):
            return np.polyval(self.ddp[(i)],l/self.Ldl[1])*np.polyval(self.ddp[(j)],l/self.Ldl[1])

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
            dddphi = sm.diff(ddphi,l)
            self.fun_dddphi = sm.lambdify(l, dddphi, modules='numpy')
            self.dddphi_true = self.fun_dddphi(self.l_all_true)

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

                plt.subplots(1,4,figsize = (20,6))
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
                plt.ylabel("$\\varphi^{'}(l,t=0)$ [deg/mm]",fontsize=20)
                plt.grid(True)

                plt.subplot(144)
                plt.plot(self.l_all_true,np.rad2deg(self.ddphi_true))
                plt.title("given beam shape (ddphi)",fontsize=15)
                plt.xlabel("$l$ [mm]",fontsize=15)
                plt.ylabel("$\\varphi^{''}(l,t=0)$ [deg/mm^2]",fontsize=20)
                plt.grid(True)
                plt.tight_layout()
                plt.show()

        def phi_approx_preparing(self,disp):
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
                if disp:
                    print("Found numpy zip archive with preparing data: psi vectors. Checking if we can use it!")
                with np.load('psi_vectors.npz') as npzfile: # for closign after using it
                    self.psi = npzfile['psi']
                    self.dpsi = npzfile['dpsi']
                    self.ddpsi = npzfile['ddpsi']
                    self.dddpsi = npzfile['dddpsi']
                    self.psi_dl = npzfile['psi_dl']
                    self.dpsi_dl = npzfile['dpsi_dl']
                    self.ddpsi_dl = npzfile['ddpsi_dl']
                    self.dddpsi_dl = npzfile['dddpsi_dl']
                    a_size = npzfile['a_size']
                    N = npzfile['N']
                    Ne = npzfile['Ne']
                    dl = npzfile['dl']
                    step = npzfile['step']
                del npzfile
            
            if (not flag_preparing_already_done) or (not N==self.N) or (not Ne==self.Ne) or\
                  (not dl==self.dl) or (not step==self.step) or (not self.a_size==a_size):
                if flag_preparing_already_done:
                    if disp:
                        print("Checking finished. We cannot use this data as FEM or/and Ldivide parameters mismatch. Creating new one:")
                self.c1 = self.E*self.I/(self.rho*self.A)
                self.EI = self.E*self.I
                if disp:
                    start_time = time.time_ns()
                    time.sleep(0.000001) # sleep 1 us

                self.psi = self.__get_psi(self.step)
                self.dpsi = self.__get_dpsi(self.step)
                self.ddpsi = self.__get_ddpsi(self.step)
                self.dddpsi = self.__get_dddpsi(self.step)
                self.psi = self.__diag_shift_mat(self.psi,self.Ne)
                self.dpsi = self.__diag_shift_mat(self.dpsi,self.Ne)
                self.ddpsi = self.__diag_shift_mat(self.ddpsi,self.Ne)
                self.dddpsi = self.__diag_shift_mat(self.dddpsi,self.Ne)
                self.index = np.array([])
                for i in range(self.Ne-1):
                    self.index = np.append(self.index,self.steps_per_fe+1+(self.steps_per_fe+1)*i) 
                self.index = np.int16(self.index)
                self.psi = np.delete(self.psi, self.index,axis=0)
                self.dpsi = np.delete(self.dpsi, self.index,axis=0)
                self.ddpsi = np.delete(self.ddpsi, self.index,axis=0)
                self.dddpsi = np.delete(self.dddpsi, self.index,axis=0)
                self.ddpsi[-1] = np.zeros((1,self.a_halfsize*(self.Ne+1)))[0]

                self.index_dl = np.int16(np.array([]))
                for l in self.l_all_dl:
                    self.index_dl = np.append(self.index_dl,self.__search_index(self.l_all_true,l))
                self.psi_dl = self.psi[self.index_dl]
                self.dpsi_dl = self.dpsi[self.index_dl]
                self.ddpsi_dl = self.ddpsi[self.index_dl]
                self.dddpsi_dl = self.dddpsi[self.index_dl]
                self.ddpsi_dl[-1] = np.zeros((1,self.a_halfsize*(self.Ne+1)))[0]

                if disp:
                    time_end = time.time_ns()-start_time-1*1e3
                    print("Preparing time: %s s" % (round(time_end*1e-9,3)))

                np.savez('psi_vectors.npz',psi=self.psi,dpsi=self.dpsi,\
                    ddpsi=self.ddpsi,dddpsi=self.dddpsi,N=self.N,Ne=self.Ne,step=self.step,dl=self.dl,\
                    psi_dl=self.psi_dl,dpsi_dl=self.dpsi_dl,\
                    ddpsi_dl=self.ddpsi_dl,dddpsi_dl=self.dddpsi_dl,\
                    a_size=self.a_size)
            else:
                if flag_preparing_already_done:
                    if disp:
                        print("Checking finished. Using loaded data")
                
        def phi_approx(self,disp=True,der_num=0,SPACAR=False):
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
                if disp:
                    print("Optimization wasn't. Don't have an approximation!We will use a created by create_a fun!")
                flag_a_approx_is = 0
                try:
                    self.a
                except:
                    raise ValueError("Call create_a first!") from None
                self.phi_approx_preparing(disp)
            else:
                if disp:
                    print("Found an approximation. Will use it!")
                self.phi_approx_preparing(disp)
                self.a = self.a_approx
                flag_a_approx_is = 1
                
            if SPACAR:
                flag_SPACAR_done = 0
                if os.path.isfile('SPACAR.npz'):
                    flag_SPACAR_done = 1

                if flag_SPACAR_done:
                    if disp:
                        print("Found SPACAR zip archive with saved data. We'll use it!")
                    with np.load('SPACAR.npz') as npzfile: # for closign after using it
                        x_SPACAR = npzfile['x']
                        y_SPACAR = npzfile['y']
                        self.phi_SPACAR_end = npzfile['phi_end']
                    del npzfile
                else:
                    if disp:
                        print("SPACAR data didn't find.")

            if der_num == 2:
                #evaluation
                if disp:
                    start_time = time.time_ns()
                    time.sleep(0.000001) # sleep 1 us
                phi_appr = np.matmul(self.psi,self.a)
                dphi_appr = np.matmul(self.dpsi,self.a)
                ddphi_appr = np.matmul(self.ddpsi,self.a)
                dddphi_appr = np.matmul(self.dddpsi,self.a)
                cos_phi_appr = np.cos(phi_appr)
                sin_phi_appr = np.sin(phi_appr)
                x = -self.step+np.cumsum(cos_phi_appr)*self.step
                y = np.cumsum(sin_phi_appr)*self.step

                phi_appr_dl = np.matmul(self.psi_dl,self.a)
                self.phi_end = phi_appr_dl[-1]  
                dphi_appr_dl = np.matmul(self.dpsi_dl,self.a)
                ddphi_appr_dl = np.matmul(self.ddpsi_dl,self.a)
                dddphi_appr_dl = np.matmul(self.dddpsi_dl,self.a)
                index = np.int16(np.array([]))
                for l in self.Ldl:
                    if l==0:
                        index = np.append(index,0)
                    elif l==self.L:
                        index = np.append(index,self.N-1)
                    else:
                        index = np.append(index,self.__search_index(self.l_all_true,l))
                x_dl = x[index]
                y_dl = y[index]
                if disp:
                    end_time = time.time_ns()-start_time-1*1e3
                    if end_time==0:
                        print("evaluation time is less then 1 ns")
                    else:
                        print("evaluation time: %s ms" % (round(end_time*1e-6,3)))
                        print("time for 1 step: %s us" % (round(1e-3*end_time/self.N,3)))
                    if SPACAR:
                        print("phi end, fleGODynamics:{} [deg]; SPACAR:{} [deg]".format(round(np.rad2deg(phi_appr[-1]),2),\
                                                                                np.round(self.phi_SPACAR_end,2)))
                    else:
                        print("phi end, fleGODynamics:{} [deg]".format(round( np.rad2deg(phi_appr[-1]),2) ))
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

            if not flag_a_approx_is:
                error = np.sum(np.power(phi_appr-self.phi_true,2))/self.N
                derror = np.sum(np.power(dphi_appr-self.dphi_true,2))/self.N
                dderror = np.sum(np.power(ddphi_appr-self.ddphi_true,2))/self.N
                ddderror = np.sum(np.power(dddphi_appr-self.dddphi_true,2))/self.N
                print("error of phi approx={}".format( error ))
                print("error of dphi approx={}".format( derror ))
                print("error of ddphi approx={}".format( dderror ))
                print("error of dddphi approx={}".format( ddderror ))

            plt.subplots(3,2,figsize = (20,12))
            plt.subplot(3,2,1)
            labels = ['$\\varphi^{true}$','$\\varphi^{approx}$']
            colours = ['b','r']
            plt.plot(self.l_all_true/self.mult,np.rad2deg(phi_appr),label=labels[1],color=colours[1])
            plt.plot(self.Ldl/self.mult,np.rad2deg(phi_appr_dl),"og")
            if not flag_a_approx_is:
                plt.plot(self.l_all_true/self.mult,np.rad2deg(self.phi_true),"--",label=labels[0],color=colours[0])
            plt.grid(True)
            plt.xlabel("$l$ [mm]")
            plt.ylabel("$\\varphi(l,t=0)$")
            plt.legend(fontsize="15",loc='best')
            if flag_a_approx_is:
                plt.title("phi approx")
            else:
                plt.title("phi approx and true")

            plt.subplot(322)
            labels = ['$\\varphi_{l}^{true}$','$\\varphi_{l}^{approx}$']
            colours = ['b','r']
            if der_num == 1 or der_num == 2:
                plt.plot(self.l_all_true,dphi_appr,label=labels[1],color=colours[1])
                plt.plot(self.Ldl,dphi_appr_dl,"og")
                if not flag_a_approx_is:
                    plt.plot(self.l_all_true,self.dphi_true,"--",label=labels[0],color=colours[0])    
            plt.grid(True)
            plt.xlabel("$l$ [mm]")
            plt.ylabel("$\\varphi_{l}(l,t=0)$",fontsize=15)
            plt.legend(fontsize="15",loc='best')
            if flag_a_approx_is:
                plt.title("dphi approx")
            else:
                plt.title("dphi approx and true")

            plt.subplot(323)
            labels = ['$\\varphi_{ll}^{true}$','$\\varphi_{ll}^{approx}$']
            colours = ['b','r']
            if der_num == 2:
                plt.plot(self.l_all_true/self.mult,ddphi_appr,label=labels[1],color=colours[1])
                plt.plot(self.Ldl/self.mult,ddphi_appr_dl,"og")
                if not flag_a_approx_is:
                    plt.plot(self.l_all_true/self.mult,self.ddphi_true,"--",label=labels[0],color=colours[0])
            plt.grid(True)
            plt.xlabel("$l$ [mm]")
            plt.ylabel("$\\varphi_{ll}(l,t=0)$",fontsize=15)
            plt.legend(fontsize="15",loc='best')
            if flag_a_approx_is:
                plt.title("ddphi approx")
            else:
                plt.title("ddphi approx and true")

            plt.subplot(325)
            labels = ['$\\varphi_{lll}^{true}$','$\\varphi_{lll}^{approx}$']
            colours = ['b','r']
            if der_num == 2:
                plt.plot(self.l_all_true/self.mult,dddphi_appr,label=labels[1],color=colours[1])
                plt.plot(self.Ldl/self.mult,dddphi_appr_dl,"og")
                if not flag_a_approx_is:
                    plt.plot(self.l_all_true/self.mult,self.dddphi_true,"--",label=labels[0],color=colours[0])
            plt.grid(True)
            plt.xlabel("$l$ [mm]")
            plt.ylabel("$\\varphi_{lll}(l,t=0)$",fontsize=15)
            plt.legend(fontsize="15",loc='best')
            if flag_a_approx_is:
                plt.title("dddphi approx")
            else:
                plt.title("dddphi approx and true")

            plt.subplot(324)
            plt.axis("off")
            plt.subplot(326)
            plt.axis("off")

            plt.subplot(3,2,(4,6))
            labels = ['$(x,y)^{true}$','$(x,y)^{approx}$','$(x,y)^{SPACAR}$']
            colours = ['b','r','b']
            plt.plot(x/self.mult,y/self.mult,label=labels[1],color=colours[1])
            plt.plot(x_dl/self.mult,y_dl/self.mult,"og")
            if not flag_a_approx_is:
                plt.plot(self.x_phi_true/self.mult,self.y_phi_true/self.mult,"--",label=labels[0],color=colours[0])
            if SPACAR:
                plt.plot(x_SPACAR*1e3,y_SPACAR*1e3,"--",label=labels[2],color=colours[2])
            plt.grid(True)
            plt.xlabel("$x$ [mm]")
            plt.ylabel("$y$ [mm]")
            plt.legend(fontsize="15",loc='best')
            plt.axis('equal') 
            if flag_a_approx_is:
                plt.title("x,y approx")
            else:
                plt.title("x,y approx and true")

            plt.tight_layout()
            plt.show()