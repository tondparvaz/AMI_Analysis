# __________________________________________________________________________________________________ 
# Angular Momentum Integral Analysis for TBL with pressure gradient
# Armin Kianfar - Fall 2023
# University of California
# Numerical dataset: 1X>Flatplate - X=0>ZPG, X=1>B1, X=2>B2, X=3>m13, X=4>m16, X=5>m18,
#                    2>Wing, 
#                    3>Gaussian bump.
# __________________________________________________________________________________________________   
# Load Libraries
import numpy as np
import os
import matplotlib.pyplot as plt
import Derivatives as der
import Integration as integ

class AMI_Analysis:
    def __init__(self, dat_ind, xs, xe, l_ind, dlta_per, window):
        self.dat_ind = dat_ind
        self.xs = xs
        self.xe = xe
        self.l_ind = l_ind
        self.dlta_per = dlta_per
        self.window = window
    
    def initialize(self): 
        if self.dat_ind == 10:
            print('> Flat plate with zero pressure gradient.')
            self.dir = "/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/data/FP/ZPG/"
            self.odir = "ZPG/"
            self.markerInd = 30
            # Plotting parameters
            self.lin = '-'
            self.Xlbl = r'$x/L$'
            self.lgnd = 'ZPG'
            self.col = (0, 0, 0)
            self.X1 = 0.2
            self.X2 = 0.9
            self.Y1 = -7.
            self.Y2 = 5
            
            # Fluid properties
            self.nu = 0.0022222222
            self.rho = 1.0
            
            # Datapoint position
            self.xi = np.loadtxt(self.dir + 'xi.txt')
            self.nx = len(self.xi)
            # normal-tangent coordinate
            self.y = np.loadtxt(self.dir + 'y.txt')
            self.ny = len(self.y)
            # x = np.loadtxt(self.dir + 'x.txt')
            self.x = np.zeros((self.nx, 1))
            self.x[:, 0] = self.xi
            self.Lscl = np.max(self.xi)
            
        elif self.dat_ind == 11:
            print('> Flat plate B1.')
            self.dir = "/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/data/FP/B1/"
            self.odir = "B1/"
            self.markerInd = 30
            # Plotting parameters
            self.lin = '-'
            self.Xlbl = r'$x/L$'
            self.lgnd = r'FP-$\beta$1'
            self.col = (0.9290, 0.6940, 0.1250)
            self.X1 = 0.2
            self.X2 = 0.9
            self.Y1 = -7.
            self.Y2 = 5
            
            # Fluid properties
            self.nu = 0.0022222222
            self.rho = 1.0
            
            # Datapoint position
            self.xi = np.loadtxt(self.dir + 'xi.txt')
            self.nx = len(self.xi)
            # normal-tangent coordinate
            self.y = np.loadtxt(self.dir + 'y.txt')
            self.ny = len(self.y)
            # x = np.loadtxt(self.dir + 'x.txt')
            self.x = np.zeros((self.nx, 1))
            self.x[:, 0] = self.xi   
            self.Lscl = np.max(self.xi)
            
        elif self.dat_ind == 12:
            print('> Flat plate m13.')
            self.dir = "/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/data/FP/m13/"
            self.odir = "m13/"
            self.markerInd = 30
            # Plotting parameters
            self.lin = '-'
            self.Xlbl = r'$x/L$'
            self.lgnd = 'FP-m13'
            self.col = (0.4940, 0.1840, 0.5560)
            self.X1 = 0.2
            self.X2 = 0.9
            self.Y1 = -7.
            self.Y2 = 5
            
            # Fluid properties
            self.nu = 0.0022222222
            self.rho = 1.0
            
            # Datapoint position
            self.xi = np.loadtxt(self.dir + 'xi.txt')
            self.nx = len(self.xi)
            # normal-tangent coordinate
            self.y = np.loadtxt(self.dir + 'y.txt')
            self.ny = len(self.y)
            # x = np.loadtxt(self.dir + 'x.txt')
            self.x = np.zeros((self.nx, 1))
            self.x[:, 0] = self.xi  
            self.Lscl = np.max(self.xi)
            
        elif self.dat_ind == 13:
            print('> Flat plate B2.')
            self.dir = "/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/data/FP/B2/"
            self.odir = "B2/"
            self.markerInd = 30
            # Plotting parameters
            self.lin = '-'
            self.Xlbl = r'$x/L$'
            self.lgnd = r'FP-$\beta$2'
            self.col = (0.4660, 0.6740, 0.1880)
            self.X1 = 0.2
            self.X2 = 0.9
            self.Y1 = -7.
            self.Y2 = 5
            
            # Fluid properties
            self.nu = 0.0022222222
            self.rho = 1.0
            
            # Datapoint position
            self.xi = np.loadtxt(self.dir + 'xi.txt')
            self.nx = len(self.xi)
            # normal-tangent coordinate
            self.y = np.loadtxt(self.dir + 'y.txt')
            self.ny = len(self.y)
            # x = np.loadtxt(self.dir + 'x.txt')
            self.x = np.zeros((self.nx, 1))
            self.x[:, 0] = self.xi    
            self.Lscl = np.max(self.xi)
            
        elif self.dat_ind == 14:
            print('> Flat plate m16.')
            self.dir = "/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/data/FP/m16/"
            self.odir = "m16/"
            self.markerInd = 30
            # Plotting parameters
            self.lin = '-'
            self.Xlbl = r'$x/L$'
            self.lgnd = 'FP-m16'
            self.col = (0.3010, 0.7450, 0.9330)
            self.X1 = 0.2
            self.X2 = 0.9
            self.Y1 = -7.
            self.Y2 = 5
            
            # Fluid properties
            self.nu = 0.0022222222
            self.rho = 1.0
            
            # Datapoint position
            self.xi = np.loadtxt(self.dir + 'xi.txt')
            self.nx = len(self.xi)
            # normal-tangent coordinate
            self.y = np.loadtxt(self.dir + 'y.txt')
            self.ny = len(self.y)
            # x = np.loadtxt(self.dir + 'x.txt')
            self.x = np.zeros((self.nx, 1))
            self.x[:, 0] = self.xi        
            self.Lscl = np.max(self.xi)
        
        elif self.dat_ind == 15:
            print('> Flat plate m18.')
            self.dir = "/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/data/FP/m18/"
            self.odir = "m18/"
            self.markerInd = 30
            # Plotting parameters
            self.lin = '-'
            self.Xlbl = r'$x/L$'
            self.lgnd = 'FP-m18'
            self.col = (0.6350, 0.0780, 0.1840)
            self.X1 = 0.2
            self.X2 = 0.9
            self.Y1 = -7.
            self.Y2 = 5
            
            # Fluid properties
            self.nu = 0.0022222222
            self.rho = 1.0
            
            # Datapoint position
            self.xi = np.loadtxt(self.dir + 'xi.txt')
            self.nx = len(self.xi)
            # normal-tangent coordinate
            self.y = np.loadtxt(self.dir + 'y.txt')
            self.ny = len(self.y)
            # x = np.loadtxt(self.dir + 'x.txt')
            self.x = np.zeros((self.nx, 1))
            self.x[:, 0] = self.xi     
            self.Lscl = np.max(self.xi) 
            
        elif self.dat_ind == 20:
            print('> No-control wing dataset, Re_c = 400,000, is called.')
            self.dir = "/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/data/Wing0/"
            self.odir = "Wing0/"
            self.markerInd = 50
            # Plotting parameters
            self.lin = '-'
            self.Xlbl = r'$x/c$'
            self.lgnd = 'Wing'
            self.col = (0, 0.4470, 0.7410)
            self.X1 = 0.2
            self.X2 = 0.9
            self.Y1 = -7.
            self.Y2 = 5
            
            # Fluid properties
            self.ReL = 400000
            self.nu = 2.5000e-06
            self.rho = 1.0
            self.Lscl = 1.0
            self.Uscl = self.ReL*self.nu/self.Lscl
            
            # Datapoint position
            self.xi = np.loadtxt(self.dir + 'xi.txt')
            self.nx = len(self.xi)
            # normal-tangent coordinate
            self.y = np.loadtxt(self.dir + 'y.txt')
            self.ny = len(self.y)
            x = np.loadtxt(self.dir + 'x.txt')
            self.x = x.reshape((self.nx, self.ny)) 
            
        elif self.dat_ind == 21:
            print('> Suction wing dataset, Re_c = 400,000, is called.')
            self.dir = "/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/data/WingS/"
            self.odir = "WingS/"
            self.markerInd = 50
            # Plotting parameters
            self.lin = '-'
            self.Xlbl = r'$x/c$'
            self.lgnd = 'W-S'
            self.col = (0, 0.4470, 0.7410)
            self.X1 = 0.2
            self.X2 = 0.9
            self.Y1 = -7.
            self.Y2 = 5
            
            # Fluid properties
            self.ReL = 400000
            self.nu = 2.5000e-06
            self.rho = 1.0
            self.Lscl = 1.0
            self.Uscl = self.ReL*self.nu/self.Lscl
            
            # Datapoint position
            self.xi = np.loadtxt(self.dir + 'xi.txt')
            self.nx = len(self.xi)
            # normal-tangent coordinate
            self.y = np.loadtxt(self.dir + 'y.txt')
            self.ny = len(self.y)
            x = np.loadtxt(self.dir + 'x.txt')
            self.x = x.reshape((self.nx, self.ny)) 
            
        elif self.dat_ind == 22:
            print('> Blowing wing dataset, Re_c = 400,000, is called.')
            self.dir = "/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/data/WingB/"
            self.odir = "WingB/"
            self.markerInd = 50
            # Plotting parameters
            self.lin = '-'
            self.Xlbl = r'$x/c$'
            self.lgnd = 'W-B'
            self.col = (0, 0.4470, 0.7410)
            self.X1 = 0.2
            self.X2 = 0.9
            self.Y1 = -6.5
            self.Y2 = 5
            
            # Fluid properties
            self.ReL = 400000
            self.nu = 2.5000e-06
            self.rho = 1.0
            self.Lscl = 1.0
            self.Uscl = self.ReL*self.nu/self.Lscl
            
            # Datapoint position
            self.xi = np.loadtxt(self.dir + 'xi.txt')
            self.nx = len(self.xi)
            # normal-tangent coordinate
            self.y = np.loadtxt(self.dir + 'y.txt')
            self.ny = len(self.y)
            x = np.loadtxt(self.dir + 'x.txt')
            self.x = x.reshape((self.nx, self.ny)) 
                      
        elif self.dat_ind == 3:
            print('> Bump dataset, Re_L = 1,000,000, is called.')
            self.dir = "/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/data/Bump/"
            self.odir = "Bump/"
            self.markerInd = 200
            # Plotting parameters
            self.lin = '-'
            self.Xlbl = r'$x/L$'
            self.lgnd = 'Bump'
            self.col = (0.8500, 0.3250, 0.0980)
            self.X1 = -0.4
            self.X2 = 0.7
            self.Y1 = -7.5
            self.Y2 = 6
            
            # Fluid properties
            ReL = 1e6
            self.Uscl = 16.4
            self.Lscl = 0.9144
            self.nu = self.Uscl*self.Lscl/ReL
            self.rho = 1.0
            
            # Datapoint position
            self.xi = np.loadtxt(self.dir + 'xi.txt')/self.Lscl
            self.nx = len(self.xi)
            # normal-tangent coordinate
            self.y = np.loadtxt(self.dir + 'y.txt')
            self.ny = len(self.y)
            x = np.loadtxt(self.dir + 'x.txt')
            self.x = x.reshape((self.nx, self.ny)) 
               
    # Load the data
    def ReadData(self):
        U = np.loadtxt(self.dir + 'U.txt')
        self.U = U.reshape((self.nx, self.ny))
        
        V = np.loadtxt(self.dir + 'V.txt')
        self.V = V.reshape((self.nx, self.ny))
        
        P = np.loadtxt(self.dir + 'P.txt')
        self.P = P.reshape((self.nx, self.ny))
        
        uu = np.loadtxt(self.dir + 'uu.txt')
        self.uu = uu.reshape((self.nx, self.ny))
        
        vv = np.loadtxt(self.dir + 'vv.txt')
        self.vv = vv.reshape((self.nx, self.ny))
        
        ww = np.loadtxt(self.dir + 'ww.txt')
        self.ww = ww.reshape((self.nx, self.ny))
        
        uv = np.loadtxt(self.dir + 'uv.txt')
        self.uv = uv.reshape((self.nx, self.ny))

    # Detect the BL edge and edge quantities - Local reconstruction method           
    def BL_edge(self):
        Po = self.P + 0.5*(self.U**2 + self.V**2) # Total pressure
        self.Ui = np.sqrt(2*(np.amax(Po, axis = 1)[:, np.newaxis] - self.P) - self.V**2) # Inviscid velocity
        self.Uiw = self.Ui[:, 0] # Inviscid velocity at the wall
        
        # Detect delta_99 BL thickness
        self.dlta = np.empty(self.nx)
        self.Ue = np.empty(self.nx)
        self.Pe = np.empty(self.nx)
        for i in range(0, self.nx):
           ind2 = np.argmax(self.U[i, :]/self.Ui[i, :] > self.dlta_per)
           ind1 = ind2 - 1
           shift = self.y[ind2] - self.y[ind1]
           f1 = self.U[i, ind1]/self.Ui[i, ind1]
           f2 = self.U[i, ind2]/self.Ui[i, ind2]
           # delta_99 BL thickness
           self.dlta[i] =  self.y[ind1] + shift*((self.dlta_per-f1)/(f2-f1))
           # Edge quantities at delta_99
           self.Ue[i] = np.interp(self.dlta[i], self.y, self.U[i, :])
           self.Pe[i] = np.interp(self.dlta[i], self.y, self.P[i, :])   
    
    def BL_Thickness(self):
        Dd1 = (self.Ui - self.U)/self.Uiw[:, np.newaxis]
        self.d1 = integ.MP_Integral(self.y, self.dlta, Dd1, self.window)
        
        Dd2 = self.U*(self.Ui - self.U)/self.Uiw[:, np.newaxis]**2
        self.d2 = integ.MP_Integral(self.y, self.dlta, Dd2, self.window)
        
        Dd2y = self.V*(self.Ui - self.U)/self.Uiw[:, np.newaxis]**2
        self.d2y = integ.MP_Integral(self.y, self.dlta, Dd2y, self.window)

    # Applying AMI analysis and computing terms in the AMI equation
    def AMI_Budget(self):
        self.BL_edge()
        self.BL_Thickness()
        self.SkinFriction()
        
        if self.l_ind == 1:
            print("AMI length scale, ell, goes with delta_1")
            self.l = 1.75*self.d1
            
        elif self.l_ind == 2:
            print("AMI length scale, ell, goes with delta_2")
            self.l = 4.54*self.d2
        ltest = np.zeros([self.nx, 1])   
        ltest[:, 0] = self.l 
        ytest = np.zeros([1, self.ny])   
        ytest[0, :] = self.y 
        # mom_1 = 1 - np.outer(1/self.l, self.y)   
        mom_1 = np.zeros([self.nx, self.ny])
        mom_1 = 1-ytest/ltest
        Dd1l = mom_1*(self.Ui - self.U)/self.Uiw[:, np.newaxis]
        self.d1l = integ.MP_Integral(self.y, self.dlta, Dd1l, self.window)
        
        Dd2l = mom_1*self.U*(self.Ui - self.U)/self.Uiw[:, np.newaxis]**2
        self.d2l = integ.MP_Integral(self.y, self.dlta, Dd2l, self.window)
        
        # Laminar friction
        self.M1 = self.nu/(self.l*self.Uiw)
        self.M1N = self.M1*2/self.Cf
        
        # Turbulent torque by Reynolds shear stress
        Iuv = integ.MP_Integral(self.y, self.dlta, -self.uv, self.window)
        self.M2 = Iuv/((self.l*self.Uiw**2))
        
        # self.M22 = np.zeros(self.nx) # overlap region
        # for i in range(0, self.nx):
        #     # in_ind = np.argmax(self.y/self.dlta[i] > 0.1)
        #     # in_y = np.append(self.y[0:in_ind], 0.1*self.dlta[i])
        #     Iuv2 = integ.TZ_Integral(self.y, self.y[0] ,self.dlta[i], -self.uv[i, :], self.window)
        #     self.M22[i] = Iuv2/((self.l[i]*self.Uiw[i]**2))
        # Iuv2 = np.zeros(self.nx)
        # for i in range(0, self.nx):
        #     Iuv2[i] = integ.MP_Integral2(self.y, self.y[0], self.dlta[i], -self.uv[i, :], self.window)
        # self.M22 = Iuv2/((self.l*self.Uiw**2))
        
        self.M2N = self.M2*2/self.Cf
        
        # Pressure gradient
        self.DUiwDx = der.FirstDerivative(self.x[:, 0], self.Uiw)
        self.M3 = self.d1l*self.DUiwDx/self.Uiw
        self.M3N = self.M3*2/self.Cf
        
        # Total mean flux
        Dd2lDx = der.FirstDerivative(self.x, self.d2l)
        DlDx = der.FirstDerivative(self.x, self.l)
            # streamwise growth
        self.M4_1 = Dd2lDx + (self.d2l-self.d2)/self.l*DlDx + 2*self.d2l*self.DUiwDx/self.Uiw
            # mean wall-normal
        self.M4_2 = self.d2y/self.l
        self.M4 = self.M4_1 + self.M4_2
        self.M4N = self.M4*2/self.Cf
        
        # Negligible terms
            # streamwise viscous diffusion
        D2UDx2 = der.SecondDerivative(self.x, self.U)
        Ivis_x = integ.MP_Integral(self.y, self.dlta, (mom_1*self.nu*D2UDx2), self.window)
        self.M5_1 = Ivis_x/self.Uiw**2
        
            # streamwise turbuelnt torque
        DuuDx = der.FirstDerivative(self.x, self.uu)
        Iturb_x = integ.MP_Integral(self.y, self.dlta, (-mom_1*DuuDx), self.window)
        self.M5_2 = Iturb_x/self.Uiw**2
        self.M5 = self.M5_1 + self.M5_2
        self.M5N = self.M5*2/self.Cf
        
        #     # Pressure difference
        # DPdiffDx = der.FirstDerivative(self.x, self.Pe[:, np.newaxis]-self.P)
        # IPdiff = integ.MP_Integral(self.y, self.dlta, (mom_1*DPdiffDx), self.window)
        # self.M5_3 = IPdiff/self.Uiw**2
        # self.M5X = self.M5_1 + self.M5_2 + self.M5_3
        # self.M5XN = self.M5X*2/self.Cf
        
            # Non-zero wall-normal contribution
        self.M6 = -self.V[:, 0]/self.Uiw
        self.M6N = self.M6*2/self.Cf
        
        # Total AMI budget
        self.AMI = self.M1 + self.M2 + self.M3 + self.M4 + self.M5
        self.res = 0.5*self.Cf - self.AMI
        self.resN = self.res*2/self.Cf
        
    # Turb term AMI analysis
    def AMI_Turb(self):
        self.BL_edge()
        self.BL_Thickness()
        self.SkinFriction()
        self.Reynolds()
        self.in_M2 = np.zeros(self.nx) # inner layer
        self.out_M2 = np.zeros(self.nx) # outer layer
        self.vw_M2 = np.zeros(self.nx) # viscous wall region
        self.ol_M2 = np.zeros(self.nx) # overlap region
        self.ll_M2 = np.zeros(self.nx) # log-law region
        for i in range(0, self.nx):
            # in_ind = np.argmax(self.y/self.dlta[i] > 0.1)
            # in_y = np.append(self.y[0:in_ind], 0.1*self.dlta[i])
            in_Iuv = integ.TZ_Integral(self.y, self.y[0] ,0.1*self.dlta[i], -self.uv[i, :], 1)
            self.in_M2[i] = in_Iuv/((self.l[i]*self.Uiw[i]**2))
            
            vw_Iuv = integ.TZ_Integral(self.y, self.y[0] ,50*self.dnu[i], -self.uv[i, :], 1)
            self.vw_M2[i] = vw_Iuv/((self.l[i]*self.Uiw[i]**2))
        
            out_Iuv = integ.TZ_Integral(self.y,  50*self.dnu[i], self.dlta[i], -self.uv[i, :], 1)
            self.out_M2[i] = out_Iuv/((self.l[i]*self.Uiw[i]**2))
            
            ol_Iuv = integ.TZ_Integral(self.y,  50*self.dnu[i], 0.1*self.dlta[i], -self.uv[i, :], 1)
            self.ol_M2[i] = ol_Iuv/((self.l[i]*self.Uiw[i]**2))
            
            ll_Iuv = integ.TZ_Integral(self.y,  30*self.dnu[i], 0.3*self.dlta[i], -self.uv[i, :], 1)
            self.ll_M2[i] = ll_Iuv/((self.l[i]*self.Uiw[i]**2))
            
            def_dir = '/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/PythonScript/Results/'
        if self.l_ind == 1:
            dir2 = 'delta1/'
            Re = self.Re1
            Xlbl = r"$Re_{\delta_1}$"
            X1 = 1400
            X2 = 6500
            
        elif self.l_ind == 2:
            dir2 = 'delta2/'
            Re = self.Re2
            Xlbl = r"$Re_{\delta_2}$"
            X1 = 1500
            X2 = 3500
        FileName = def_dir + dir2   
        os.makedirs(FileName, exist_ok=True) 
        ylim_1 = 0
        ylim_2 = 5
        fig_lineW = 0.75
        #___________________________________________________________________________________________________________________
        plt.figure(40, figsize=(3.2,3.1))
        plt.plot(Re[self.xs:self.xe], self.in_M2[self.xs:self.xe]*1e3, self.lin, linewidth=fig_lineW, color=self.col, label=self.lgnd)
        plt.plot(Re[self.xs:self.xe], self.M2[self.xs:self.xe]*1e3, '--', linewidth=0.5, color=self.col, label=self.lgnd)
        plt.xlim(X1, X2)
        plt.ylim(ylim_1, ylim_2)
        plt.xlabel(Xlbl, fontsize=10)
        plt.ylabel(r"$\frac{1}{\ell}\int_{y^*=0}^{y^*=0.1} \frac{-\overline{u^\prime v^\prime}}{ U_{io}^2} \mathrm{d}y$ $\left(\times 10^3 \right)$", fontsize=10)
        # plt.ylabel(r"$\frac{1}{\ell}\int_{y_1}^{y_2} \frac{-\overline{u^\prime v^\prime}}{ U_{io}^2} \mathrm{d}y$ $\left(\times 10^3 \right)$", fontsize=10)        
        # legend = plt.legend(fontsize=8)
        # legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        # plt.legend(loc='upper right')  
        F2 = FileName + 'in_M2.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        # plt.close(1)  
        
        plt.figure(41, figsize=(3.2,3.1))
        plt.plot(Re[self.xs:self.xe], self.out_M2[self.xs:self.xe]*1e3, self.lin, linewidth=fig_lineW, color=self.col, label=self.lgnd)
        plt.plot(Re[self.xs:self.xe], self.M2[self.xs:self.xe]*1e3, '--', linewidth=0.5, color=self.col, label=r"log-law region")
        plt.xlim(X1, X2)
        plt.ylim(ylim_1, ylim_2)
        plt.xlabel(Xlbl, fontsize=10)
        # plt.ylabel(r"$\frac{1}{\ell}\int_{y_1}^{y_2} \frac{-\overline{u^\prime v^\prime}}{ U_{io}^2} \mathrm{d}y$ $\left(\times 10^3 \right)$", fontsize=10) 
        plt.ylabel(r"$\frac{1}{\ell}\int_{y^+=50}^{y^*=1} \frac{-\overline{u^\prime v^\prime}}{ U_{io}^2} \mathrm{d}y$ $\left(\times 10^3 \right)$", fontsize=10)
        # legend = plt.legend(fontsize=8)
        # legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque)     
        plt.legend(loc='lower left', fontsize=8)      
        F2 = FileName + 'out_M2.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        # plt.close(2)  
        
        plt.figure(42, figsize=(3.2,3.1))
        plt.plot(Re[self.xs:self.xe], self.vw_M2[self.xs:self.xe]*1e3, self.lin, linewidth=fig_lineW, color=self.col)
        plt.plot(Re[self.xs:self.xe], self.M2[self.xs:self.xe]*1e3, '--', linewidth=0.5, color=self.col, label=self.lgnd)
        plt.xlim(X1, X2)
        plt.ylim(ylim_1, ylim_2)
        plt.xlabel(Xlbl, fontsize=10)
        plt.ylabel(r"$\frac{1}{\ell}\int_{y^*=0}^{y^+=50} \frac{-\overline{u^\prime v^\prime}}{ U_{io}^2} \mathrm{d}y$ $\left(\times 10^3 \right)$", fontsize=10)
        # plt.ylabel(r"$\frac{1}{\ell}\int_{y_1}^{y_2} \frac{-\overline{u^\prime v^\prime}}{ U_{io}^2} \mathrm{d}y$ $\left(\times 10^3 \right)$", fontsize=10) 
        F2 = FileName + 'vw_M2.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        # plt.close(3)  
        
        plt.figure(43, figsize=(3.2,3.1))
        plt.plot(Re[self.xs:self.xe], self.ol_M2[self.xs:self.xe]*1e3, self.lin, linewidth=fig_lineW, color=self.col)
        plt.plot(Re[self.xs:self.xe], self.M2[self.xs:self.xe]*1e3, '--', linewidth=0.5, color=self.col, label=self.lgnd)
        plt.xlim(X1, X2)
        plt.ylim(ylim_1, ylim_2)
        plt.xlabel(Xlbl, fontsize=10)
        # plt.ylabel(r"$\frac{1}{\ell}\int_{y_1}^{y_2} \frac{-\overline{u^\prime v^\prime}}{ U_{io}^2} \mathrm{d}y$ $\left(\times 10^3 \right)$", fontsize=10) 
        plt.ylabel(r"$\frac{1}{\ell}\int_{y^+=50}^{y^*=0.1} \frac{-\overline{u^\prime v^\prime}}{ U_{io}^2} \mathrm{d}y$ $\left(\times 10^3 \right)$", fontsize=10)
        F2 = FileName + 'ol_M2.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        
        plt.figure(44, figsize=(3.2,3.1))
        plt.plot(Re[self.xs:self.xe], self.ll_M2[self.xs:self.xe]*1e3, self.lin, linewidth=fig_lineW, color=self.col)
        plt.plot(Re[self.xs:self.xe], self.M2[self.xs:self.xe]*1e3, '--', linewidth=0.5, color=self.col, label=self.lgnd)
        plt.xlim(X1, X2)
        plt.ylim(ylim_1, ylim_2)
        plt.xlabel(Xlbl, fontsize=10)
        # plt.ylabel(r"$\frac{1}{\ell}\int_{y_1}^{y_2} \frac{-\overline{u^\prime v^\prime}}{ U_{io}^2} \mathrm{d}y$ $\left(\times 10^3 \right)$", fontsize=10) 
        plt.ylabel(r"$\frac{1}{\ell}\int_{y^+=30}^{y^*=0.3} \frac{-\overline{u^\prime v^\prime}}{ U_{io}^2} \mathrm{d}y$ $\left(\times 10^3 \right)$", fontsize=10)
        F2 = FileName + 'll_M2.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        # plt.close(4)   

    # Skin friction coefficient and wall units    
    def SkinFriction(self):
        self.tauw = self.rho*self.nu*((self.U[:, 1]-self.U[:, 0])/(self.y[1]-self.y[0]))
        # self.tauw = self.rho*self.nu*((-3*self.U[:, 0]+4*self.U[:, 1]-self.U[:, 2])/(2*(self.y[1]-self.y[0])))
        self.utau = np.sqrt(self.tauw/self.rho)
        self.dnu = self.nu/self.utau
        
        self.Up = self.U/self.utau[:, np.newaxis]
        self.uup = self.uu/self.utau[:, np.newaxis]**2
        self.vvp = self.vv/self.utau[:, np.newaxis]**2
        self.wwp = self.ww/self.utau[:, np.newaxis]**2
        self.uvp = self.uv/self.utau[:, np.newaxis]**2
        
        self.Cf = 2*self.tauw/(self.rho*self.Uiw**2)
        
    def error(self):
        self.eps = (1 - self.AMI/(self.Cf/2))*100
        print("> Average error for the AMI equation:", np.mean(self.eps), "%")
    
    # Reynolds numbers  
    def Reynolds(self):
        self.BL_edge()
        self.BL_Thickness()
        self.SkinFriction()
        self.Re1 = self.d1*self.Uiw/self.nu
        self.Re2 = self.d2*self.Uiw/self.nu
        self.Retau = self.d1*self.utau/self.nu
        self.Rex = self.xi*self.Uiw/self.nu
        
    # Compute Clauser parameter    
    def Beta(self):
        self.B = -2.0*self.d1/self.Cf/self.Uiw*self.DUiwDx
        self.Bl = -2.0*self.d1l/self.Cf/self.Uiw*self.DUiwDx
        
        self.H12 = self.d1/self.d2
        self.G12 = (self.H12-1.0)/(self.H12*np.sqrt(self.Cf/2))
        
        self.H12l = self.d1l/self.d2l
        self.G12l = (self.H12l-1.0)/(self.H12l*np.sqrt(self.Cf/2))
        
        self.Reynolds()
        def_dir = '/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/PythonScript/Results/'
        if self.l_ind == 1:
            dir2 = 'delta1/'
            Re = self.Re1
            Xlbl = r"$Re_{\delta_1}$"
            X1 = 1400
            X2 = 6500
            
        elif self.l_ind == 2:
            dir2 = 'delta2/'
            Re = self.Re2
            Xlbl = r"$Re_{\delta_2}$"
            X1 = 1500
            X2 = 3500
        FileName = def_dir + dir2   
        os.makedirs(FileName, exist_ok=True) 
        
        fig_lineW = 0.75
        plt.figure(10, figsize=(3.2,3.1))
        plt.semilogy(self.Re1[self.xs:self.xe], self.B[self.xs:self.xe], self.lin, linewidth=fig_lineW, color=self.col, label=self.lgnd)
        plt.xlim(X1, X2)
        plt.ylim(1e-1, 50)
        plt.xlabel(Xlbl, fontsize=10)
        plt.ylabel(r"$\beta$", fontsize=10)
        # legend = plt.legend(fontsize=8)
        # legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        F2 = FileName + 'B.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        plt.legend(loc='lower left')  
        F2 = FileName + 'B.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        
        plt.figure(11, figsize=(3.2,3.1))
        plt.semilogy(self.Re1[self.xs:self.xe], self.Bl[self.xs:self.xe], self.lin, linewidth=fig_lineW, color=self.col, label=self.lgnd)
        plt.xlim(X1, X2)
        plt.ylim(1e-1, 50)
        plt.xlabel(Xlbl, fontsize=10)
        plt.ylabel(r"$\beta_\ell$", fontsize=10)
        # legend = plt.legend(fontsize=8)
        # legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        F2 = FileName + 'Bl.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        plt.legend(loc='lower left')  
        
        plt.figure(12, figsize=(3.2,3.1))
        plt.semilogy(self.Re1[self.xs:self.xe], self.G12[self.xs:self.xe], self.lin, linewidth=fig_lineW, color=self.col, label=self.lgnd)
        plt.xlim(X1, X2)
        plt.ylim(1e-1, 50)
        plt.xlabel(Xlbl, fontsize=10)
        plt.ylabel(r"$G$", fontsize=10)
        # legend = plt.legend(fontsize=8)
        # legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        F2 = FileName + 'G.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        plt.legend(loc='lower left') 
        
        plt.figure(13, figsize=(3.2,3.1))
        plt.semilogy(self.Re1[self.xs:self.xe], self.G12l[self.xs:self.xe], self.lin, linewidth=fig_lineW, color=self.col, label=self.lgnd)
        plt.xlim(X1, X2)
        plt.ylim(1e-1, 50)
        plt.xlabel(Xlbl, fontsize=10)
        plt.ylabel(r"$G_\ell$", fontsize=10)
        # legend = plt.legend(fontsize=8)
        # legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        F2 = FileName + 'Gl.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        plt.legend(loc='lower left') 
        # plt.close(1)  
        
        
    # Stats - History effects
    
    # Plotting - AMI budget with x
    def Plotting1(self):
        if self.l_ind == 1:
            dir2 = 'delta1/'
            
        elif self.l_ind == 2:
            dir2 = 'delta2/'
            
        def_dir = '/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/PythonScript/Results/'
        FileName = def_dir + self.odir + dir2
        os.makedirs(FileName, exist_ok=True)
        fig_lineW = 0.75
        # figsizeW = 3.2
        # figsizeH = 3.1
        # if self.dat_ind == 3:
        figsizeW = 5.8
        figsizeH = 3.1
        #___________________________________________________________________________________________________________________
        plt.figure(1, figsize=(figsizeW,figsizeH))
        # plt.plot(self.xi, self.Cf/2*1e3, 'k', linewidth=fig_lineW, label=r'$C_f/2$')
        # plt.plot(self.xi, self.M1*1e3, self.lin, linewidth=fig_lineW, color=self.col, marker='o', markersize=4, markevery=self.markerInd, label='laminar friction')
        # plt.plot(self.xi, self.M2*1e3, self.lin, linewidth=fig_lineW, color=self.col, marker='s', markersize=4,markevery=self.markerInd, label='turbulent torque')
        # plt.plot(self.xi, self.M3*1e3, self.lin, linewidth=fig_lineW, color=self.col, marker='^', markersize=4,markevery=self.markerInd, label='pressure Grad')
        # plt.plot(self.xi, self.M4*1e3, self.lin, linewidth=fig_lineW, color=self.col, marker='D', markersize=4,markevery=self.markerInd, label='mean flux')
        plt.plot(self.xi, self.res*1e3, ':', linewidth=fig_lineW, color=self.col, marker='v', markersize=4,markevery=self.markerInd, label='residual')
        # plt.plot(self.xi, self.M5*1e3, self.lin, linewidth=fig_lineW, color=self.col, marker='x',markevery=self.markerInd, label='negligible')
        if self.dat_ind == 21:
            plt.plot(self.xi, self.M6*1e3, self.lin, linewidth=fig_lineW, color='r', label='wall BC - suction')
        elif self.dat_ind == 22:
            plt.plot(self.xi, self.M6*1e3, self.lin, linewidth=fig_lineW, color='r', label='wall BC - blowing')
        # plt.plot(self.xi, self.AMI*1e3, '--k', linewidth=fig_lineW)
        # plt.plot([self.xi[0], self.xi[-1]], [0, 0], 'k', linewidth=0.5)
        plt.xlim(self.X1, self.X2)
        plt.ylim(self.Y1, self.Y2)
        plt.xlabel(self.Xlbl, fontsize=10)
        plt.ylabel(r'AMI budget $\left(\times 10^3 \right)$', fontsize=10)
        legend = plt.legend(fontsize=8)
        legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque)
        FileName = def_dir + self.odir + dir2 + "AMI_x.eps"
        plt.savefig(FileName, format='eps', bbox_inches='tight') 
        # plt.close(1)
        
        plt.figure(2, figsize=(figsizeW,figsizeH))
        # plt.plot(self.xi, self.M1N, self.lin, linewidth=fig_lineW, color=self.col, marker='o', markersize=4,markevery=self.markerInd, label='Laminar friction')
        # plt.plot(self.xi, self.M2N, self.lin, linewidth=fig_lineW, color=self.col, marker='s', markersize=4,markevery=self.markerInd, label='Turbulent torque')
        # plt.plot(self.xi, self.M3N, self.lin, linewidth=fig_lineW, color=self.col, marker='^', markersize=4,markevery=self.markerInd, label='Pressure Grad')
        # plt.plot(self.xi, self.M4N, self.lin, linewidth=fig_lineW, color=self.col, marker='D', markersize=4,markevery=self.markerInd, label='Mean flux')
        plt.plot(self.xi, self.resN, ':', linewidth=fig_lineW,color=self.col, marker='v', markersize=4,markevery=self.markerInd, label='residual')
        # plt.plot(self.xi, self.M5N, self.lin, linewidth=fig_lineW, color=self.col, marker='x',markevery=self.markerInd, label='Negligible terms')
        if self.dat_ind == 21:
            plt.plot(self.xi, self.M6N, self.lin, linewidth=fig_lineW, color='r', label='wall BC - suction')
        elif self.dat_ind == 22:
            plt.plot(self.xi, self.M6N, self.lin, linewidth=fig_lineW, color='r', label='wall BC - blowing')
        # plt.plot([self.xi[0], self.xi[-1]], [0, 0], 'k', linewidth=0.5)
        plt.xlim(self.X1, self.X2)
        plt.ylim(self.Y1, self.Y2)
        plt.xlabel(self.Xlbl, fontsize=10)
        plt.ylabel(r'Normalized AMI budget', fontsize=10)
        # plt.legend(fontsize=8)
        FileName = def_dir + self.odir + dir2 + "AMI_x_N.eps"
        plt.savefig(FileName, format='eps', bbox_inches='tight')
        # plt.close(2)
        
    def Plotting2(self):
        self.AMI_Budget()
        self.SkinFriction()
        self.Reynolds()
        def_dir = '/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/PythonScript/Results/'
        if self.l_ind == 1:
            dir2 = 'delta1/'
            Re = self.Re1
            Xlbl = r"$Re_{\delta_1}$"
            X1 = 1400
            X2 = 6500
            
        elif self.l_ind == 2:
            dir2 = 'delta2/'
            Re = self.Re2
            Xlbl = r"$Re_{\delta_2}$"
            X1 = 1500
            X2 = 3500
        FileName = def_dir + dir2   
        os.makedirs(FileName, exist_ok=True) 
        fig_lineW = 0.75
        #___________________________________________________________________________________________________________________
        plt.figure(1, figsize=(3.2,3.1))
        plt.plot(Re[self.xs:self.xe], self.M1[self.xs:self.xe]*1e3, self.lin, linewidth=fig_lineW, color=self.col, label=self.lgnd)
        plt.plot([Re[0], Re[-1]], [0, 0], 'k', linewidth=0.5)
        plt.xlim(X1, X2)
        plt.ylim(self.Y1, self.Y2)
        plt.xlabel(Xlbl, fontsize=10)
        plt.ylabel(r"$\frac{1}{Re_{\ell}}$ $\left(\times 10^3\right)$", fontsize=10)
        legend = plt.legend(fontsize=8)
        legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        plt.legend(loc='lower right')  
        F2 = FileName + 'M1.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        # plt.close(1)  
        
        plt.figure(2, figsize=(3.2,3.1))
        plt.plot(Re[self.xs:self.xe], self.M2[self.xs:self.xe]*1e3, self.lin, linewidth=fig_lineW, color=self.col, label=self.lgnd)
        plt.plot([Re[0], Re[-1]], [0, 0], 'k', linewidth=0.5)
        plt.xlim(X1, X2)
        plt.ylim(self.Y1, self.Y2)
        plt.xlabel(Xlbl, fontsize=10)
        plt.ylabel(r"$\frac{1}{\ell}\int_0^\infty \frac{-\overline{u^\prime v^\prime}}{ U_{io}^2} \mathrm{d}y$ $\left(\times 10^3 \right)$", fontsize=10)
        # legend = plt.legend(fontsize=8)
        # legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque)     
        plt.legend(loc='lower left')      
        F2 = FileName + 'M2.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        # plt.close(2)  
        
        plt.figure(3, figsize=(3.2,3.1))
        plt.plot(Re[self.xs:self.xe], self.M3[self.xs:self.xe]*1e3, self.lin, linewidth=fig_lineW, color=self.col)
        plt.plot([Re[0], Re[-1]], [0, 0], 'k', linewidth=0.5)
        plt.xlim(X1, X2)
        plt.ylim(self.Y1, self.Y2)
        plt.xlabel(Xlbl, fontsize=10)
        plt.ylabel(r'$\frac{\delta_1^\ell}{U_{io}}\frac{\mathrm{d}U_{io}}{\mathrm{d}x}$ $\left(\times 10^3 \right)$', fontsize=10)
        F2 = FileName + 'M3.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        # plt.close(3)  
        
        plt.figure(4, figsize=(3.2,3.1))
        plt.plot(Re[self.xs:self.xe], self.M4[self.xs:self.xe]*1e3, self.lin, linewidth=fig_lineW, color=self.col)
        plt.plot([Re[0], Re[-1]], [0, 0], 'k', linewidth=0.5)
        plt.xlim(X1, X2)
        plt.ylim(self.Y1, self.Y2)
        plt.xlabel(Xlbl, fontsize=10)
        plt.ylabel(r'$\frac{\mathrm{d}\delta_2^\ell}{\mathrm{d}x}+\frac{\delta_2^\ell - \delta_2}{\ell} \frac{\mathrm{d}\ell}{\mathrm{d}x}+\frac{2\delta_2^\ell}{U_{io}}\frac{\mathrm{d}U_{io}}{\mathrm{dx}}+ \frac{\delta_{2,v}}{\ell}$ $\left(\times 10^3 \right)$', fontsize=10)
        F2 = FileName + 'M4.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        # plt.close(4)  
        
        plt.figure(5, figsize=(3.2,3.1))
        plt.plot(Re[self.xs:self.xe], (self.M4[self.xs:self.xe]+self.M3[self.xs:self.xe])*1e3, self.lin, linewidth=fig_lineW, color=self.col)
        plt.plot([Re[0], Re[-1]], [0, 0], 'k', linewidth=0.5)
        plt.xlim(X1, X2)
        plt.ylim(self.Y1, self.Y2)
        plt.xlabel(Xlbl, fontsize=10)
        plt.ylabel(r'$\frac{\delta_1^\ell}{U_{io}}\frac{\mathrm{d}U_{io}}{\mathrm{d}x} + \frac{\mathrm{d}\delta_2^\ell}{\mathrm{d}x}+\frac{\delta_2^\ell - \delta_2}{\ell} \frac{\mathrm{d}\ell}{\mathrm{d}x}+\frac{2\delta_2^\ell}{U_{io}}\frac{\mathrm{d}U_{io}}{\mathrm{dx}}+ \frac{\delta_{2,v}}{\ell}$ $\left(\times 10^3 \right)$', fontsize=10)
        F2 = FileName + 'M5.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        

    # Skin friction coefficient and wall units    
    def ProfPlot(self, direc, val, Ulim, uulim1, uulim2, marker):
        self.Reynolds()
        if self.l_ind == 1:
            ReChoice = self.Re1[self.xs:self.xe]
        elif self.l_ind == 2:
            ReChoice = self.Re2[self.xs:self.xe]
            
        UiwX = np.interp(val, ReChoice, self.Uiw[self.xs:self.xe])
        lX = np.interp(val, ReChoice, self.l[self.xs:self.xe])
        utauX = np.interp(val, ReChoice, self.utau[self.xs:self.xe])
        dnuX = np.interp(val, ReChoice, self.dnu[self.xs:self.xe])
        d99X = np.interp(val, ReChoice, self.dlta[self.xs:self.xe])
        CfX = np.interp(val, ReChoice, self.Cf[self.xs:self.xe])
        
        uuX = np.zeros(self.ny)
        vvX = np.zeros(self.ny)
        wwX = np.zeros(self.ny)
        uvX = np.zeros(self.ny)
        UX = np.zeros(self.ny)
        for i in range(0, self.ny): 
            uuX[i] = np.interp(val, ReChoice, self.uu[self.xs:self.xe, i])
            vvX[i] = np.interp(val, ReChoice, self.vv[self.xs:self.xe, i])
            wwX[i] = np.interp(val, ReChoice, self.ww[self.xs:self.xe, i])
            uvX[i] = np.interp(val, ReChoice, self.uv[self.xs:self.xe, i])
            UX[i] = np.interp(val, ReChoice, self.U[self.xs:self.xe, i])
        
        ## Plotting
        def_dir = '/Users/armink/Desktop/Research_UC/AMI_MTI/APG-data/AMI_Analysis/PythonScript/Results/Stats/dlta1/'
        FileName = def_dir + direc
        os.makedirs(FileName, exist_ok=True) 
        fig_lineW = 0.75
        ylim_1 = 0
        ylim_2 = 15
        
        plt.figure(51, figsize=(3.2,3.1))
        plt.plot(self.y/lX, UX/UiwX, self.lin, linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=self.lgnd)
        # plt.plot([Re[0], Re[-1]], [0, 0], 'k', linewidth=0.5)
        plt.xlim(0, 5)
        plt.ylim(0, 1.5)
        plt.xlabel(r"$y/\ell$", fontsize=10)
        plt.ylabel(r"$\overline{u}/U_{io}$", fontsize=10)
        legend = plt.legend(fontsize=8)
        legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        # plt.legend(loc='lower right')  
        F2 = FileName + 'U_1.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight')
         
        
        plt.figure(52, figsize=(3.2,3.1))
        plt.plot(self.y/d99X, UX/UiwX, self.lin, linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=self.lgnd)
        # plt.plot([Re[0], Re[-1]], [0, 0], 'k', linewidth=0.5)
        plt.xlim(0, 1.5)
        plt.ylim(0, 1.5)
        plt.xlabel(r"$y/\delta_{99}$", fontsize=10)
        plt.ylabel(r"$\overline{u}/U_{io}$", fontsize=10)
        legend = plt.legend(fontsize=8)
        legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        # plt.legend(loc='lower right')  
        F2 = FileName + 'U_2.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        
        plt.figure(53, figsize=(3.2,3.1))
        plt.plot(self.y/lX, UX/utauX, self.lin, linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=self.lgnd)
        # plt.plot([Re[0], Re[-1]], [0, 0], 'k', linewidth=0.5)
        plt.xlim(0, 5)
        plt.ylim(0, Ulim)
        plt.xlabel(r"$y/\ell$", fontsize=10)
        plt.ylabel(r"$\overline{u}^+$", fontsize=10)
        legend = plt.legend(fontsize=8)
        legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        # plt.legend(loc='lower right')  
        F2 = FileName + 'U_3.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight')
        
        plt.figure(54, figsize=(3.2,3.1))
        plt.plot(self.y/d99X, UX/utauX, self.lin, linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=self.lgnd)
        # plt.plot([Re[0], Re[-1]], [0, 0], 'k', linewidth=0.5)
        plt.xlim(0, 1.5)
        plt.ylim(0, Ulim)
        plt.xlabel(r"$y/\delta_{99}$", fontsize=10)
        plt.ylabel(r"$\overline{u}^+$", fontsize=10)
        legend = plt.legend(fontsize=8)
        legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        # plt.legend(loc='lower right')  
        F2 = FileName + 'U_4.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight')
        
        plt.figure(55, figsize=(3.2,3.1))
        plt.semilogx(self.y/dnuX, UX/utauX, self.lin, linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=self.lgnd)
        np.savetxt('yasd.txt', self.y/dnuX, fmt='%.6f', delimiter='\t')
        np.savetxt('uasd.txt', UX/utauX, fmt='%.6f', delimiter='\t')
        # plt.plot([Re[0], Re[-1]], [0, 0], 'k', linewidth=0.5)
        plt.xlim(1e-1, 3000)
        plt.ylim(0, Ulim)
        plt.xlabel(r"$y^+$", fontsize=10)
        plt.ylabel(r"$\overline{u}^+$", fontsize=10)
        legend = plt.legend(fontsize=8)
        legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        # plt.legend(loc='lower right')  
        F2 = FileName + 'U_5.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        
        
        plt.figure(61, figsize=(3.2,3.1))
        plt.plot(self.y/lX, 1e3*uuX/UiwX**2, '-', linewidth=fig_lineW, color=self.col, label=r"$\overline{u^\prime u^\prime}$")
        plt.plot(self.y/lX, 1e3*vvX/UiwX**2, '--', linewidth=fig_lineW, color=self.col, label=r"$\overline{v^\prime v^\prime}$")
        plt.plot(self.y/lX, 1e3*wwX/UiwX**2, '-.', linewidth=fig_lineW, color=self.col, label=r"$\overline{w^\prime w^\prime}$")
        plt.plot(self.y/lX, 1e3*uvX/UiwX**2, ':', linewidth=fig_lineW, color=self.col, label=r"$\overline{u^\prime v^\prime}$")
        plt.plot([0, 5], [0, 0], 'k', linewidth=0.25)
        plt.xlim(0, 5)
        plt.ylim(uulim1, uulim2)
        plt.xlabel(r"$y/\ell$", fontsize=10)
        plt.ylabel(r"$\overline{u^\prime_i u^\prime_j}/U_{io}^2 \: (\times 10^3)$", fontsize=10)
        legend = plt.legend(fontsize=8)
        legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        plt.legend(loc='lower right')  
        F2 = FileName + 'uu_1.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        
        plt.figure(62, figsize=(3.2,3.1))
        plt.plot(self.y/d99X, 1e3*uuX/UiwX**2, '-', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=self.lgnd)
        plt.plot(self.y/d99X, 1e3*vvX/UiwX**2, '--', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=self.lgnd)
        plt.plot(self.y/d99X, 1e3*wwX/UiwX**2, '-.', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=self.lgnd)
        plt.plot(self.y/d99X, 1e3*uvX/UiwX**2, ':', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=self.lgnd)
        plt.plot([0, 5], [0, 0], 'k', linewidth=0.25)
        plt.xlim(0, 1.5)
        plt.ylim(uulim1, uulim2)
        plt.xlabel(r"$y/\delta_{99}$", fontsize=10)
        plt.ylabel(r"$\overline{u^\prime_i u^\prime_j}/U_{io}^2 \: (\times 10^3)$", fontsize=10)
        # legend = plt.legend(fontsize=8)
        # legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        # plt.legend(loc='lower right')  
        F2 = FileName + 'uu_2.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight') 
        
        plt.figure(63, figsize=(3.2,3.1))
        plt.plot(self.y/lX, uuX/utauX**2, '-', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=self.lgnd)
        plt.plot(self.y/lX, vvX/utauX**2, '--', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=self.lgnd)
        plt.plot(self.y/lX, wwX/utauX**2, '-.', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=self.lgnd)
        plt.plot(self.y/lX, uvX/utauX**2, ':', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=self.lgnd)
        plt.plot([0, 5], [0, 0], 'k', linewidth=0.25)       # plt.plot([Re[0], Re[-1]], [0, 0], 'k', linewidth=0.5)
        plt.xlim(0, 5)
        plt.ylim(uulim1, uulim2)
        plt.xlabel(r"$y/\ell$", fontsize=10)
        plt.ylabel(r"$\overline{u^\prime_i u^\prime_j}^+$", fontsize=10)
        # legend = plt.legend(fontsize=8)
        # legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        # plt.legend(loc='lower right')  
        F2 = FileName + 'uu_3.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight')
        
        plt.figure(64, figsize=(3.2,3.1))
        plt.plot(self.y/d99X, uuX/utauX**2, '-', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=r"$\overline{u^\prime u^\prime}$")
        plt.plot(self.y/d99X, vvX/utauX**2, '--', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=r"$\overline{v^\prime v^\prime}$")
        plt.plot(self.y/d99X, wwX/utauX**2, '-.', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=r"$\overline{w^\prime w^\prime}$")
        plt.plot(self.y/d99X, uvX/utauX**2, ':', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=r"$\overline{u^\prime v^\prime}$")
        plt.plot([0, 5], [0, 0], 'k', linewidth=0.25)       # plt.plot([Re[0], Re[-1]], [0, 0], 'k', linewidth=0.5)
        plt.xlim(0, 1.5)
        plt.ylim(uulim1, uulim2)
        plt.xlabel(r"$y/\delta_{99}$", fontsize=10)
        plt.ylabel(r"$\overline{u^\prime_i u^\prime_j}^+$", fontsize=10)
        # legend = plt.legend(fontsize=8)
        # legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque) 
        # plt.legend(loc='lower right')  
        F2 = FileName + 'uu_4.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight')
        
        plt.figure(65, figsize=(3.2,3.1))
        plt.semilogx(self.y/dnuX, uuX/utauX**2, '-', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=r"$\overline{u^\prime u^\prime}$")
        plt.semilogx(self.y/dnuX, vvX/utauX**2, '--', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=r"$\overline{v^\prime v^\prime}$")
        plt.semilogx(self.y/dnuX, wwX/utauX**2, '-.', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=r"$\overline{w^\prime w^\prime}$")
        plt.semilogx(self.y/dnuX, uvX/utauX**2, ':', linewidth=fig_lineW, color=self.col, marker=marker, markersize=4, markevery=self.markerInd, label=r"$\overline{u^\prime v^\prime}$")
        plt.semilogx([0, 1000], [0, 0], 'k', linewidth=0.25) 
        plt.xlim(1e-1, 3000)
        plt.ylim(uulim1, uulim2)
        plt.xlabel(r"$y^+$", fontsize=10)
        plt.ylabel(r"$\overline{u^\prime_i u^\prime_j}^+$", fontsize=10)
        # plt.legend(loc='lower right')  
        F2 = FileName + 'uu_5.eps'
        plt.savefig(F2, format='eps', bbox_inches='tight')
        

        
    # def Plotting3(self, profID):
    #     if self.dat_ind == 3:
    #         plt.xlabel(r"$y/L$", fontsize=10)
    #         lbl1 = r"$y/L=$"
    #     else:
    #         plt.xlabel(r"$y/c$", fontsize=10)
    #         lbl1 = r"$y/c=$"
    #     fig_lineW = 0.75
    #     markers = ['o', 's', 'D']
    #     n = len(profID)
    #     for i in range(0, n):
    #         ind = profID[i]
    #         markID = markers[i]
    #         lbl = lbl1 + "{:.2f}".format(self.xi[ind])
    #         plt.figure(31, figsize=(3.2,3.1))
    #         plt.semilogx(self.y/self.Lscl, self.U[ind, :]/self.Uscl, self.y/self.Lscl, self.Ui[ind, :]/self.Uscl, '--',linewidth=fig_lineW, color=self.col
    #                     , marker=markID,markevery=self.markerInd, label=lbl)
    #         plt.semilogx([self.dlta[ind]/self.Lscl, self.dlta[ind]/self.Lscl], [0.2, 1], 'r', linewidth=0.5)
    #         plt.ylabel(r"$\overline{U}/U_\infty$", fontsize=10)
        # legend = plt.legend(fontsize=8)
        # legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque)
        # plt.xlim(1e-2, 1)
        # plt.ylim(0, 1)