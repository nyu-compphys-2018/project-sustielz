import numpy as np;
import matplotlib.pyplot as plt;
import poiseuille as pois;
import time as time;


def make_poiseuille(N, M, P0, dP, L, eta, inc):
    P1 = P0 + dP;
    P2 = P0 - dP;
    if(inc):
        gas = pois.poiseuille_incompressible(N, M, P1, P2, L, eta);
    else:
        gas = pois.poiseuille(N, M, P1, P2, L, eta);
    s = gas.s;

    rho = (1.0/gas.cs**2)*np.ones([N+2, M+2]);
    rho[-1, :] = (P1/gas.cs**2)*np.ones(M+2);
    rho[N, :] = (P2/gas.cs**2)*np.ones(M+2);
    vx = np.zeros([N+2,M+2]);
    vy = np.zeros([N+2,M+2]);
    for I in range(9):
        gas.updateEq(rho, vx, vy);
        gas.f[I] = gas.feq[I];
    return gas, s;


def plotProfile(gas):
##        plt.scatter([10, 20, 30, 40, 50, 60, 70], [0.09, 0.15, 0.2, 0.2376, 0.276, 0.305, 0.328]);
##        plt.show();
        h = gas.h;
        N, M, L = gas.N, gas.M, gas.L;
        H = (M+1)*float(L)/(N+1);
        cp = h/gas.s;
        Yvals = np.arange(-1, M+1);
        
        X, Y = gas.X, gas.Y;
##        print X;
##        print Y;
##
##        print self.L*float(M+1)/(N+1)
        U = gas.vx;
        plt.plot(Y, U[N/2, :][Yvals], label="vx (y profile)");
        EXACT =  0.5*(gas.cs**2)*(gas.P2 - gas.P1)*Y*(Y - H)/(gas.eta) ;
##        EXACT = ( (self.h/self.s)*(self.P2 - self.P1)*(self.cs**2)/(2*self.eta) )*(self.L**2-4Y**2);
        plt.plot(Y, EXACT, label="exact");
        plt.title("t = {}".format(gas.t));
        plt.legend();
        return np.sqrt(sum( (U[N/2, :][Yvals] - EXACT)**2 ))/M;

def convergencePlot_incompressible(full):
    P0 = 1.0;
    dP = 0.05;
    P1 = P0 - dP;
    P2 = P0 + dP;
    N0 = 11;
    M0 = 6;
    L = 40;
    eta = 5.0;
    L2 = [];
    Nrange = np.zeros(6);
    Mrange = 1.0*Nrange;
    for i in range(6):
        
        N = N0*(i+1) - (i+1);
        M = M0*(i+1) - (i+1);
        Nrange[i] = N;
        Mrange[i] = M;
        
        gas, s = make_poiseuille(N, M, P0, dP, L, eta, False)
        gas.set_FULL(full);
        
        N = gas.N;
        M = gas.M;
    
        err = [];
        count = 0;
        TOL = 1e-5;
        myerr = 10.0;

        start = time.time();
##        while(gas.t < 1600*s):
        while(myerr - TOL> 1e-18):
            vxtemp = gas.vx;
            vytemp = gas.vy;
            gas.iterate(s);
########            if(count%30==0):
########                plt.subplot(211);
########                gas.plotProfile();
########                plt.subplot(212);
########                gas.heatmap(gas.vx);
########                
########                plt.pause(0.01);
########                plt.clf();

########                plt.subplot(211);
########                plt.imshow(gas.getImage(omega), cmap='hot', origin='lower', interpolation='nearest');
########                plt.colorbar();
########                plt.subplot(212);
########                plt.streamplot(np.arange(N+2), np.arange(M+2), gas.getImage(gas.vx), gas.getImage(gas.vy), density=1);
########                plt.pause(0.01);
########                plt.clf();


            
##                print " "
            if(count%1000==0):
                print "t = {}, error = {}, error - TOL = {}".format(gas.t, myerr, myerr - TOL);
            count+= 1;
            if((sum(sum(abs(gas.vx) + abs(gas.vy))))>1e-14):
                myerr = (sum(sum(abs(gas.vx - vxtemp) + abs(gas.vy - vytemp))))/(sum(sum(abs(gas.vx) + abs(gas.vy)))) ;
                err.append(myerr);
        print "final error = ", myerr;
        stop = time.time();
        print "the time taken to run was", stop-start;
        myL2 = plotProfile(gas);
        plt.show();
        L2.append(myL2);
##        plt.show();
##        plt.plot(Yvals, gas.vx[N/2, :][Yvals], label="vx (y profile)");
##        plt.plot(Yvals, EXACT, label="exact");
##        
##        plt.legend();
##        plt.show();  
        print "N = ", N, ", M = ", M; 
        print "L2 is ", myL2;
        print " ";
    plt.clf();
    plt.plot(np.log(Nrange*Mrange), np.log(L2));
    plt.show();

def Trapezoid(a, b, Y, N): 
    h = (b - a)/N;
    Ny = size(Y);  #Note that in general, when we call Trapezoid from Romberg, Y has more elements than N since we don't use every element at first 
    I = h/2*(Y[0] + Y[Ny-1]);
    for i in range((Ny-1)/N, Ny-1, (Ny-1)/N): #Note that step size depends on N
        I += h*Y[i];
return I;
        

def convergencePlot(full):
    ###### Parameters of the Flow ####
    P0 = 1.0;
    dP = 0.05;
    P1 = P0 - dP;
    P2 = P0 + dP;
    N0 = 4;
    M0 = 3;
    L = 40;
    eta = 5.0;
    L2 = [];

    ###### Find the "exact" solution by running code at very high resolution

    I = 9;
    ex, s = make_poiseuille(N0**I - 2, M0**I - 2, P0, dP, L, eta, False)
    ex.set_FULL(full);
    TOL = 1e-5;
    count = 0;
    myerr = 10;
    
    start = time.time();
    while(myerr - TOL> 1e-18):
        vxtemp = ex.vx;
        vytemp = ex.vy;
        ex.iterate(s);           
        if(count%100==0):
            print "t = {}, error = {}, error - TOL = {}".format(ex.t, myerr, myerr - TOL);
        count+= 1;
        if((sum(sum(abs(ex.vx) + abs(ex.vy))))>1e-14):
            myerr = (sum(sum(abs(ex.vx - vxtemp) + abs(ex.vy - vytemp))))/(sum(sum(abs(ex.vx) + abs(ex.vy)))) ;

    stop = time.time();
    print "high-resolution solution calculated."
    print "converged to steady flow within error = ", myerr;
    print "the time taken to run was", stop-start;
    print " "
    Yvals = np.arange(-1, M0**I - 2 + 1);
    EXACT = ex.vx[(N0**I)/2, :][Yvals];
    print ex.vx[(N0**I)/2, :]
    print EXACT
    myL2 = plotProfile(ex);
    plt.show();    


    ###### Find solutions and compare to high res
    Nrange = np.zeros(I-2);
    Mrange = 1.0*Nrange;
    for i in range(I-2):
        
        n = N0**(i+1) - 2;
        m = M0**(i+1) - 2;
        Nrange[i] = n;
        Mrange[i] = m;
        
        gas, s = make_poiseuille(n, m, P0, dP, L, eta, False)
        gas.set_FULL(full);
        
        N = gas.N;
        M = gas.M;
    
        err = [];
        count = 0;
        TOL = 1e-5;
        myerr = 10.0;
        
        start = time.time();
##        while(gas.t < 1600*s):
        while(myerr - TOL> 1e-18):
            vxtemp = gas.vx;
            vytemp = gas.vy;
            gas.iterate(s);
########            if(count%30==0):
########                plt.subplot(211);
########                gas.plotProfile();
########                plt.subplot(212);
########                gas.heatmap(gas.vx);
########                
########                plt.pause(0.01);
########                plt.clf();

########                plt.subplot(211);
########                plt.imshow(gas.getImage(omega), cmap='hot', origin='lower', interpolation='nearest');
########                plt.colorbar();
########                plt.subplot(212);
########                plt.streamplot(np.arange(N+2), np.arange(M+2), gas.getImage(gas.vx), gas.getImage(gas.vy), density=1);
########                plt.pause(0.01);
########                plt.clf();


            
##                print " "
            if(count%1000==0):
                print "t = {}, error = {}, error - TOL = {}".format(gas.t, myerr, myerr - TOL);
            count+= 1;
            if((sum(sum(abs(gas.vx) + abs(gas.vy))))>1e-14):
                myerr = (sum(sum(abs(gas.vx - vxtemp) + abs(gas.vy - vytemp))))/(sum(sum(abs(gas.vx) + abs(gas.vy)))) ;
                err.append(myerr);
        print "converged to steady flow within error = ", myerr;
        stop = time.time();
        print "the time taken to run was", stop-start;
        
        Yvals = np.arange(-1, M+1);
        myprofile = gas.vx[N/2, :][Yvals];
        myEXACT = EXACT[np.arange(0, M0**I - 2 + 1, M0**(I - i-1))];
        print np.arange(0, M0**I - 2 + 1, M0**(I - i-1))
        print np.size(EXACT), np.size(myEXACT);
        myL2 = np.sqrt(sum( (myprofile - myEXACT)**2 ))/M;
                                  
##        plt.show();
##        plt.plot(Yvals, gas.vx[N/2, :][Yvals], label="vx (y profile)");
##        plt.plot(Yvals, EXACT, label="exact");
##        
##        plt.legend();
##        plt.show();  
        print "N = ", N, ", M = ", M; 
        print "L2 is ", myL2;
        print " ";
    plt.clf();
    print Nrange
    print L2
    plt.plot(np.log(Nrange), np.log(L2));
    plt.show();
if __name__ == '__main__':
##    A = 10*np.arange(4**3 + 1);
##    print A
##    B = np.arange(5);
##    print B;
##    print B*4**2
##    print A[np.arange(0, (4**3) + 1, 4**1)]
##    convergencePlot_incompressible(True);
    convergencePlot(True);
        
