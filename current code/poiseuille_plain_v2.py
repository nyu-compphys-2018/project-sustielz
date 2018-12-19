import numpy as np;
import matplotlib.pyplot as plt;
import poiseuille as pois;
import time as time;
import linfit as linfit;
from scipy.integrate import simps;


def make_poiseuille(N, M, P0, dP, L, eta, inc):
    P1 = P0 + dP;
    P2 = P0 - dP;
    if(inc):
        gas = pois.poiseuille_incompressible(N, M, P1, P2, L, eta);
    else:
        gas = pois.poiseuille(N, M, P1, P2, L, eta);


    rho = (1.0/gas.cs**2)*np.ones([N+2, M+2]);
    rho[-1, :] = (P1/gas.cs**2)*np.ones(M+2);
    rho[N, :] = (P2/gas.cs**2)*np.ones(M+2);
    vx = np.zeros([N+2,M+2]);
    vy = np.zeros([N+2,M+2]);
    for I in range(9):
        gas.updateEq(rho, vx, vy);
        gas.f[I] = gas.feq[I];
    return gas, gas.s;


def plotProfile(gas):
##        plt.scatter([10, 20, 30, 40, 50, 60, 70], [0.09, 0.15, 0.2, 0.2376, 0.276, 0.305, 0.328]);
##        plt.show();
        h = gas.h;
        N, M, L = gas.N, gas.M, gas.L;
        if(gas.full):
            H = (M+2)*float(L)/(N+2);
        else:
            H = (M)*float(L)/(N);
        cp = h/gas.s;
        Yvals = np.arange(-1, M+1);
     
        X, Y = gas.X, gas.Y;
##        print X;
##        print Y;
##        print self.L*float(M+1)/(N+1)
        U = gas.vx;
        plt.plot(Y, U[N/2, :][Yvals], label="vx (y profile)");
        EXACT =  0.25*(gas.cs**2)*(gas.P2 - gas.P1)*Y*(Y - H)/(gas.eta) ;
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
    N0 = 10;
    M0 = 8;
    L = 40;
    eta = 5.0;
    L2 = [];
    I = 7;
    Nrange = np.zeros(I);
    Mrange = 1.0*Nrange;
    for i in [4]:
        
        N = N0*(i+1);
        M = M0*(i+1);
        Nrange[i] = N;
        Mrange[i] = M;
        
        gas, s = make_poiseuille(N, M, P0, dP, L, eta, True)
        gas.set_FULL(False);
        
        
        
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
            if(count%30==0):
                plt.subplot(211);
                plotProfile(gas);
                plt.subplot(212);
                gas.heatmap(gas.rho);
                
                plt.pause(0.01);
                plt.clf();


            
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



def convergencePlot(full):
    ###### Parameters of the Flow ####
    P0 = 1.0;
    dP = 0.05;
    P1 = P0 - dP;
    P2 = P0 + dP;
    N0 = 10;
    M0 = 8;
    L = 40;
    H = L*float(N0+1)/(M0+1)
    eta = 5.0;
    L2 = [];
    profiles = [];

    ###### Find the "exact" solution by running code at very high resolution

    I = 7;
    ex, s = make_poiseuille(N0*(I+2), M0*(I+2), P0, dP, L, eta, False)
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
    Yvals = np.arange(-1, M0*(I+2)+1);
##    print np.size(ex.vx[N0*I/2, :][Yvals])
    EXACT = simps(ex.vx[N0*(I+2)/2, :][Yvals], ex.Y);
    print "Area under curve exact curve is ", EXACT
    print " "
    np.save("exact_half.npy", EXACT);


    ###### Find solutions and compare to high res
    Nrange = np.zeros(I-2);
    Mrange = 1.0*Nrange;
    for i in range(I-2):
        
        n = N0*(i+1);
        m = M0*(i+1);
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
            if(count%30==0):
                plt.subplot(211);
                plotProfile(gas);
                plt.subplot(212);
                gas.heatmap(gas.rho);
                
                plt.pause(0.01);
                plt.clf();

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
        AREA = simps(gas.vx[N/2, :][Yvals], gas.Y);
        myprofile = gas.vx[N/2, :][Yvals];
        L2.append(EXACT - AREA);
        print "N = ", N, ", M = ", M;
        print "Area under curve is ", AREA
        print "Area - Exact Area is ", EXACT - AREA;
        profiles.append(myprofile);
        print " ";
    plt.show();
    plt.clf();
    np.save('L2_test.npy', L2);
    np.save('nrange_test.npy', Nrange);
    plt.scatter(np.log(Nrange), np.log(L2));
    plt.show();
    np.save('profiles_test.npy', profiles);

    for i in range(np.size(profiles)):
        plt.plot(profiles[i]);
    plt.show();
    
    
if __name__ == '__main__':
##    A = 10*np.arange(4**3 + 1);
##    print A
##    B = np.arange(5);
##    print B;
##    print B*4**2
##    print A[np.arange(0, (4**3) + 1, 4**1)]
##    convergencePlot_incompressible(True);

######
######
    convergencePlot(False);
    profiles = np.load('profiles.npy');
    L2 = np.load('L2.npy');
    N = np.load('nrange.npy');
    profilesH = np.load('profiles_test.npy');
    L2H = np.load('L2_test.npy');
    NH = np.load('nrange_test.npy');
    N0 = 10;
    M0 = 8;
    L = 40;
    H = L*float(N0+1)/(M0+1)
    M = 0.8*N;
    MH = 0.8*NH;
    print N;
    print NH;
    vmax = [];
    count = 0;
    
    for Vx in profiles:
        vmax.append(max(Vx));
        Y = np.linspace(0, H, np.size(Vx));
        plt.plot(Y, Vx, label="N = {}".format(N[count]));
        count+=1;
    plt.legend();
    plt.show();

    count = 0;
    for Vx in profilesH:
        vmax.append(max(Vx));
        Y = np.linspace(0, H, np.size(Vx));
        plt.plot(Y, Vx, label="N = {}".format(N[count]));
        count+=1;
    plt.legend();
    plt.show();
##    print vmax
##    print N[np.arange(0, count)]
##    N2 = N[np.arange(0, count)]
##    plt.scatter(np.log(M[np.arange(0, count)]), np.log(max(vmax) - vmax));
##    plt.show();
        


    
    
##    plt.scatter(M, L2);
##    plt.scatter(M, L2/M);
##    plt.scatter(M, L2/M**2);
##    plt.show();
    plt.scatter(np.log(M), np.log(L2)/M);
    fit, m, c = linfit.LinFit(np.log(M), np.log(L2));
    plt.plot(np.log(M), fit);
    print m;
    plt.scatter(np.log(MH), np.log(L2H)/MH);
    fit, m, c = linfit.LinFit(np.log(MH), np.log(L2H));
    plt.plot(np.log(MH), fit);
    print m;
    plt.show();
