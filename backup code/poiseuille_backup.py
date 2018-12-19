import numpy as np;
import matplotlib.pyplot as plt;
import lbm_v3 as LBM;
import time as time;
import matplotlib.animation as animation

##pressure driven flow for D2Q9 lattice
class poiseuille(LBM.D2Q9): 
    s = 0.0;    ##initialize timestep
    h = 0.0;    ##initialize dx
       
    VX = [];
    VY = [];    
    def __init__(self, N, M, P1, P2, L, eta):
        super(poiseuille, self).__init__(N, M);
        self.P1 = P1;                   ##Pressure at inlet (1) and outlet (2)
        self.P2 = P2;
        self.L = float(L);              ##Physical length
        self.X = np.zeros(N+2);
        self.Y = np.zeros(M+2);
        self.set_FULL(False); ##Toggle full bounce back conditions. Also initialize h=dx, and initialize X and Y..
        self.set_viscosity(eta);


        
    def set_FULL(self, bool):
        N, M, L = self.N, self.M, self.L;
        if(bool):
            print "Using Full Bounce-Back";
            self.full = True;
            self.h = float(L)/(N+1);
            self.Collision = self.Collision_Full;
            self.Stream = self.Stream_Full;  
            self.X = self.h*np.arange(N+2);
            self.Y = self.h*np.arange(M+2);

        else:
            print "Using Half Bounce-Back";
            self.full = False;
            self.h = float(L)/(N); 
            self.X = self.h*(-0.5 + np.arange(N+2));
            self.Y = self.h*(-0.5 + np.arange(M+2));
         


    def set_viscosity(self, eta):       ##Update viscosity and choose timestep s accordingly
        self.eta = eta;
        self.s = (self.tau - 0.5)*(self.h**2)*(self.cs**2)/self.eta;

    def boundary(self): 
        N, M = self.N, self.M;
               ##Bounce on top and bottom walls (i.e. no-slip)
        for I in [2, 5, 6]:
            J = self.opp[I];
            for i in range(-1, N+1): 
##                print " "
##                print np.transpose(self.F[2].f);
##                print " "
##                print " "
                self.f[I][i, -1] = self.fs[J][i, -1];
                self.f[J][i, M] = self.fs[I][i, M];
##        self.bounce(self.corners);  ##Bounce corners
####        ###Zou-He no-slip walls
####        for i in range(N):
####            self.f[2][i, -1] = self.f[4][i, -1];
####            self.f[5][i, -1] = self.f[7][i, -1] - 0.5*(self.f[1][i, -1] - self.f[3][i, -1]);
####            self.f[6][i, -1] = self.f[8][i, -1] + 0.5*(self.f[1][i, -1] - self.f[3][i, -1]);
####
####            self.f[4][i, M] = self.f[2][i, M];
####            self.f[7][i, M] = self.f[5][i, M] - 0.5*(self.f[1][i, -1] - self.f[3][i, -1]);
####            self.f[8][i, M] = self.f[6][i, M] + 0.5*(self.f[1][i, -1] - self.f[3][i, -1]);
####        
####
        ######Zou-He boundary conditions to simulate pressure gradient at sides (copied from zou-he paper)
        rhoL = self.P1/self.cs**2;
        rhoR = self.P2/self.cs**2;
        vyL = 0.0;
        vyR = 0.0;
####        fL = [];
####        fR = [];
####        for I in range(9):
####            fL.append(self.F[I].fs[-1, :]);
####            fR.append(self.F[I].fs[N, :]);
        fL = np.zeros(9);
        fR = np.zeros(9);
        for j in range(-1, M+1):
            for I in range(9):
                fL[I] = self.f[I][-1, j]
                fR[I] = self.f[I][N, j]
            
            ##Left wall
            vxL = 1.0 - (1/rhoL)*(fL[0] + fL[2] + fL[4] + 2.0*(fL[3] + fL[6] + fL[7]));
            self.f[1][-1, j] = fL[3] + (2.0/3.0)*rhoL*vxL    
            self.f[5][-1, j] = fL[7] - 0.5*(fL[2] - fL[4]) + (1.0/6.0)*rhoL*vxL + 0.5*rhoL*vyL;
            self.f[8][-1, j] = fL[6] + 0.5*(fL[2] - fL[4]) + (1.0/6.0)*rhoL*vxL - 0.5*rhoL*vyL;

            ##Right wall
            vxR = -1.0 + (1/rhoR)*(fR[0] + fR[2] + fR[4] + 2.0*(fR[1] + fR[5] + fR[8])) ;
            self.f[3][N, j] = fR[1] - (2.0/3.0)*rhoR*vxR    
            self.f[7][N, j] = fR[5] + 0.5*(fR[2] - fR[4]) - (1.0/6.0)*rhoR*vxR - 0.5*rhoR*vyR;
            self.f[6][N, j] = fR[8] - 0.5*(fR[2] - fR[4]) - (1.0/6.0)*rhoR*vxR + 0.5*rhoR*vyR;

##        ##Corner nodes
##        x = 0.0;
##        self.f[1][-1, -1] = self.f[3][-1, -1];
##        self.f[5][-1, -1] = self.f[7][-1, -1];
##        self.f[2][-1, -1] = self.f[4][-1, -1];
##        for I in [0, 1, 2, 3, 4, 5, 7]:
##            x += self.f[I][-1, -1];
##        
##        self.f[6][-1, -1] = 0.5*(rhoL - x);
##        self.f[8][-1, -1] = self.f[6][-1, -1];
##
##        x = 0.0;
##        self.f[1][-1, M] = self.f[3][-1, M];
##        self.f[8][-1, M] = self.f[6][-1, M];
##        self.f[4][-1, M] = self.f[2][-1, M];
##        for I in [0, 1, 2, 3, 4, 6, 8]:
##            x += self.f[I][-1, M];
##        
##        self.f[5][-1, M] = 0.5*(rhoL - x);
##        self.f[7][-1, M] = self.f[5][-1, M];

##        x = 0.0;
##        for I in [0, 1, 2, 3, 4, 6, 8]:
##            x += self.fs[I][N, -1];
##        self.f[3][N, -1] = self.f[1][N, -1];
##        self.f[5][N, -1] = 0.5*(rhoR - x);
##        self.f[7][N, -1] = self.f[7][N, -1];   
##
##        x = 0.0;
##        for I in [0, 1, 2, 3, 4, 5, 7]:
##            x += self.fs[I][N, M];
##        self.f[3][N, M] = self.f[1][N, M];
##        self.f[6][N, M] = 0.5*(rhoR - x);
##        self.f[8][N, M] = self.f[7][N, M];
##        

    def plotFlow(self):
        h = self.h;
        X = self.X;
        Y = self.Y;
##        print Y;
        U = (self.h/self.s)*self.getImage(self.vx)
        V = (self.h/self.s)*self.getImage(self.vy)
        plt.streamplot(X, Y, U, V, density=0.6);

    def plotVorticity(self):
        h = self.h;
        X = self.X;
        Y = self.Y;
        U = (self.h/self.s)*self.vx;
        V = (self.h/self.s)*self.vy;
        omega = np.ones([self.N, self.M]);
        for i in range(-1, self.N-1):
                for j in range(-1, self.M-1):
                    omega[i+1, j+1] = (0.5/h)*( U[i+1, j] - V[i-1, j] - (V[i, j+1] - U[i, j-1]) );
        plt.imshow(np.transpose(omega), origin='lower', interpolation='nearest', extent=[0, self.L, 0, self.L*float(self.M)/self.N]);
        plt.colorbar();

    def plotProfile(self):
##        plt.scatter([10, 20, 30, 40, 50, 60, 70], [0.09, 0.15, 0.2, 0.2376, 0.276, 0.305, 0.328]);
##        plt.show();
        h = self.h;
        N, M, L = self.N, self.M, self.L;
        H = (M+1)*float(L)/(N+1);
        cp = self.h/self.s;
        Yvals = np.arange(-1, M+1);
        
        X, Y = self.X, self.Y;
##        print X;
##        print Y;
##
##        print self.L*float(M+1)/(N+1)
        U = self.vx;
        plt.plot(Y, U[N/2, :][Yvals], label="vx (y profile)");
        EXACT =  0.5*(self.cs**2)*(self.P2 - self.P1)*Y*(Y - H)/(self.eta) ;
##        EXACT = ( (self.h/self.s)*(self.P2 - self.P1)*(self.cs**2)/(2*self.eta) )*(self.L**2-4Y**2);
        plt.plot(Y, EXACT, label="exact");
        plt.title("t = {}".format(self.t));
        plt.legend();
        return np.sqrt(sum( (U[N/2, :][Yvals] - EXACT)**2 ));
        

    def heatmap(self, thing):
        h = self.h;
        X = self.X;
        Y = self.Y;
        plt.imshow(self.getImage(thing), origin='lower', interpolation='nearest', extent=[0, self.L, 0, self.L*float(self.M)/self.N]);
        plt.colorbar();

    def plotFlow2(self, i):
        h = self.h;
        X = self.X;
        Y = self.Y;
##        print Y;
        U = (self.h/self.s)*self.getImage(self.VX[i])
        V = (self.h/self.s)*self.getImage(self.VY[i])
        plt.streamplot(X, Y, U, V, density=1.4);
        plt.ylim(0, h*self.M);
        plt.xlim(0, h*self.N);
        plt.imshow(self.getImage(1-(self.fluid==0)),  extent=[0, self.L, 0, self.L*float(self.M)/self.N], alpha=0.5, cmap='gray', origin="lower", aspect='auto');
        plt.gca().set_aspect('equal');

##        plot_Circle(30, self.M/2 - 1, self.h, 15);

        
def update_anim(i, gas):
    plt.clf();
    gas.plotFlow2(i);

class poiseuille_incompressible(poiseuille):

    def updateV(self):
        self.vx = self.vx*0;
        self.vy = self.vy*0;
        for I in range(9):
            ex, ey = self.dct[I];
            self.vx += self.f[I]*ex;
            self.vy += self.f[I]*ey;

    def updateEq(self, rho, vx, vy):
        for I in range(9):
            
            ex, ey = self.dct[I];
            A = ex*vx + ey*vy;
            B = vx*vx + vy*vy;
            self.feq[I] = self.w[I]*(rho + 3.0*A + 4.5*A**2 - 1.5*B );

    def boundary(self): 
        N, M = self.N, self.M;
##        self.bounce(self.boundaries[1]);
##        self.bounce(self.boundaries[3]);
            
        ##Bounce on top and bottom walls (i.e. no-slip)
        for I in [2, 5, 6]:
            ex, ey = self.dct[I];
            J = self.opp[I];
            for i in range(-1, N+1):
                if(self.full):
                    self.f[I][i, -1] = self.fs[J][i, -1];
                    self.f[J][i, M] = self.fs[I][i, M];
                else:
                    self.f[I][i, -1] = self.fs[J][i+ex, M+ey];
                    self.f[J][i, M] = self.fs[I][i-ex, M-ey];
##        ######Zou-He No Slip Walls
##        for i in range(N):
##                self.f[2][i, -1] = self.f[4][i, -1];
##                self.f[5][i, -1] = self.f[7][i, -1] - 0.5*(self.f[1][i, -1] - self.f[3][i, -1]);
##                self.f[6][i, -1] = self.f[8][i, -1] + 0.5*(self.f[1][i, -1] - self.f[3][i, -1]);
##
##                self.f[4][i, M] = self.f[2][i, M];
##                self.f[7][i, M] = self.f[5][i, M] - 0.5*(self.f[1][i, -1] - self.f[3][i, -1]);
##                self.f[8][i, M] = self.f[6][i, M] + 0.5*(self.f[1][i, -1] - self.f[3][i, -1]);
##    
        ######Zou-He boundary conditions to simulate pressure gradient at sides (copied from zou-he paper)
        rhoL = self.P1/self.cs**2;
        rhoR = self.P2/self.cs**2;
        vyL = 0.0;
        vyR = 0.0;
####        fL = [];
####        fR = [];
####        for I in range(9):
####            fL.append(self.F[I].fs[-1, :]);
####            fR.append(self.F[I].fs[N, :]);
        fL = np.zeros(9);
        fR = np.zeros(9);
        for j in range(-1, M+1):
            for I in range(9):
                fL[I] = self.f[I][-1, j]
                fR[I] = self.f[I][N, j]
            
            ##Left wall
            vxL = rhoL - (fL[0] + fL[2] + fL[4] + 2.0*(fL[3] + fL[6] + fL[7]));
            self.f[1][-1, j] = fL[3] + (2.0/3.0)*vxL    
            self.f[5][-1, j] = fL[7] - 0.5*(fL[2] - fL[4]) + (1.0/6.0)*vxL + 0.5*vyL;
            self.f[8][-1, j] = fL[6] + 0.5*(fL[2] - fL[4]) + (1.0/6.0)*vxL - 0.5*vyL;

            ##Right wall
            vxR = -rhoR + (fR[0] + fR[2] + fR[4] + 2.0*(fR[1] + fR[5] + fR[8])) ;
            self.f[3][N, j] = fR[1] - (2.0/3.0)*vxR    
            self.f[7][N, j] = fR[5] + 0.5*(fR[2] - fR[4]) - (1.0/6.0)*vxR - 0.5*vyR;
            self.f[6][N, j] = fR[8] - 0.5*(fR[2] - fR[4]) - (1.0/6.0)*vxR + 0.5*vyR;

##        ##Corner nodes
##        x = 0.0;
##        self.f[1][-1, -1] = self.f[3][-1, -1];
##        self.f[5][-1, -1] = self.f[7][-1, -1];
##        self.f[2][-1, -1] = self.f[4][-1, -1];
##        for I in [0, 1, 2, 3, 4, 5, 7]:
##            x += self.f[I][-1, -1];
##        
##        self.f[6][-1, -1] = 0.5*(rhoL - x);
##        self.f[8][-1, -1] = self.f[6][-1, -1];
##
##        x = 0.0;
##        self.f[1][-1, M] = self.f[3][-1, M];
##        self.f[8][-1, M] = self.f[6][-1, M];
##        self.f[4][-1, M] = self.f[2][-1, M];
##        for I in [0, 1, 2, 3, 4, 6, 8]:
##            x += self.f[I][-1, M];
##        
##        self.f[5][-1, M] = 0.5*(rhoL - x);
##        self.f[7][-1, M] = self.f[5][-1, M];

##        x = 0.0;
##        for I in [0, 1, 2, 3, 4, 6, 8]:
##            x += self.fs[I][N, -1];
##        self.f[3][N, -1] = self.f[1][N, -1];
##        self.f[5][N, -1] = 0.5*(rhoR - x);
##        self.f[7][N, -1] = self.f[7][N, -1];   
##
##        x = 0.0;
##        for I in [0, 1, 2, 3, 4, 5, 7]:
##            x += self.fs[I][N, M];
##        self.f[3][N, M] = self.f[1][N, M];
##        self.f[6][N, M] = 0.5*(rhoR - x);
##        self.f[8][N, M] = self.f[7][N, M];
##
       
if __name__ == '__main__':
    number_of_frames = 500
    ###### Define parameters of channel flow
    P0 = 1.3;
    dP = 0.25;
    P1 = P0 + dP;
    P2 = P0 - dP;
    N = 120;
    M = 64;
####    N = 5;
####    M = 5;
    L = 40;
    eta = 25.0;
    
    gas = poiseuille(N, M, P1, P2, L, eta);
    h = gas.h;
    s = gas.s;
    print "s =", s

##    for SHAPE in gas.obstacles:
##        print SHAPE
##    gas.plotMask();
    rho = (1.0/gas.cs**2)*np.ones([N+2, M+2]);
    rho[-1, :] = (P1/gas.cs**2)*np.ones(M+2);
    rho[N, :] = (P2/gas.cs**2)*np.ones(M+2);
    vx = np.zeros([N+2,M+2]);
    vy = np.zeros([N+2,M+2]);

    ###### Set object mask for a half-cylinder or cylinder
    mask = np.ones([N+2, M+2]);
    r = 16;
    x0, y0 = 2*r, M/2 - 1;
    for i in range(-r, r+1):
        for j in range(-r, r+1):
##        for j in range(-r, 0):
            if( i**2 + j**2 < r**2 ):
                mask[x0+i, y0+j] = 0;
    gas.defineSolid(mask);
    gas.plotMask();
    r=r-1;
    Xcyl = h*(x0 + np.arange(-r, r+1));
    Ycyl_p = h*(y0 + np.sqrt(r**2 - np.arange(-r, r+1)**2));
##    Ycyl_p = h*y0*np.ones(np.size(Xcyl));
    Ycyl_m = h*(y0 - np.sqrt(r**2 - np.arange(-r, r+1)**2));


    ###### Do iteration and plot
    err = [];
    print gas.h;
    print gas.s;
    for I in range(9):
        gas.updateEq(rho, vx, vy);
        gas.f[I] = gas.feq[I];
    count = 0;
    start = time.time();
    while(gas.t < s*number_of_frames):
        gas.iterate(s);
        print count;
        gas.VX.append(gas.vx);
        gas.VY.append(gas.vy);

        
##        gas.plotMask();
##        plt.show();
##        gas.plotF();
##        plt.show();
        update_anim(count, gas);
        if(count%30==0):
##        if(False):
##            plt.subplot(211);
####            gas.plotVorticity();
##            gas.heatmap(gas.rho);
##            plt.plot(Xcyl, Ycyl_p, 'r');
##            plt.plot(Xcyl, Ycyl_m, 'r');
##
##            plt.subplot(212);
##            gas.plotFlow();
##            plt.plot(Xcyl, Ycyl_p, 'r');
##            plt.plot(Xcyl, Ycyl_m, 'r');
            plt.pause(0.01);
####            plt.show();
            plt.clf();
            
####        if(count%30==0):
####            gas.plotProfile();
####            plt.pause(0.01);
####            plt.clf();

##            print "rho left = ", gas.rho[-1, :];
##            print "rho right = ", gas.rho[N, :];
##            print "vx bottom = ", gas.vx[:, -1];
##            print "vx top = ", gas.vx[:, M];
##            print " " 
##            print "t = {}, error = {}, error - TOL = {}".format(gas.t, myerr, myerr - TOL);
        count+= 1;
    stop = time.time();
    print "the time taken to run was", stop-start;
##    plt.imshow(gas.getImage(gas.vx), cmap='hot', origin='lower', interpolation='nearest');
##    plt.show();
##    plt.plot(s*np.arange(np.size(err)), err);
##    plt.show();


##
##    plt.plot(gas.vx[16, :][range(-1, M-1)], label="vx (y profile)");
##    plt.plot(np.arange(32), (-0.05/(16*80))*np.arange(32)*((np.arange(32) - 32)), label="exact");
##    ##    plt.plot(gas.vx[:, 30], label="vx (x profile)");
##    plt.legend();
##    plt.show();  
    plt.show();
##    for count in range(np.size(self.obstacles)):
##        print "Obstacle {} Force Calculations".format(count+1);
##        print "Fx per unit time: ", gas.Fx[count]*s/gas.t;
##        print "Fy per unit time: ", gas.Fy[count]*s/gas.t;


    fig4=plt.figure();

    anim = animation.FuncAnimation(fig4, update_anim, number_of_frames, fargs=(gas, ));
    anim.save('cyl_Re_High.html', writer='ImageMagickWriter',fps=20)

##
        
