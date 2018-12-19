import numpy as np;
import matplotlib.pyplot as plt;
import poiseuille_hybrid as pois;
import time as time;
import linfit as linfit;

def makeCylinder(N, M, P1, P2, L, eta, x0, y0, R):
    gas = pois.poiseuille(N, M, P1, P2, L, eta);
    mask = [];
    if(x0 - R <= 0 or x0 + R >= N or y0 - R <= 0 or y0 + R >= M):
        print "Error: Cylinder touches walls"
    else:            
        for i in range(-R, R+1):
            for j in range(-R, R+1):
                if( i**2 + j**2 < R**2 ):
                     mask.append([x0+i, y0+j]);
##    gas.addObstacle(mask);
##    plt.imshow(gas.getImage(gas.fluid));
##    plt.colorbar();
##    plot_Circle(x0, y0, 1.0, R);
##    plt.show();

    rho = (1.0/gas.cs**2)*np.ones([N+2, M+2]);
    rho[-1, :] = (P1/gas.cs**2)*np.ones(M+2);
    rho[N, :] = (P2/gas.cs**2)*np.ones(M+2);
    vx = np.zeros([N+2,M+2]);
    vy = np.zeros([N+2,M+2]);

    for I in range(9):
        gas.updateEq(rho, vx, vy);
        gas.f[I] = gas.feq[I];

    return gas;




    
if __name__ == '__main__':
    Fxt, Fyt, Fxb, Fyb = np.load('WallForces_Pvar_full.npy');
    print Fxt
    print Fyt
    print Fxb
    print Fyb
    vmax = np.load('vmax_P_full.npy')
##    print Fxt
    Pvals = np.load('Pvals_full.npy'); 
    plt.title("Drag vs Pressure");
    plt.scatter(Pvals, Fxt, label="Force calculated with MEM");
    plt.plot(Pvals, Pvals*1.5, label="Analytical Solution");
##    plt.plot(Pvals, Fxt, label="Drag (top wall)");
##    plt.plot(Pvals, Fyt, label="Lift (top wall)");
##    plt.plot(Pvals, Fxb, label="Drag (bottom wall)");
##    plt.plot(Pvals, Fyb, label="Lift (bottom wall)");
 
    plt.legend();
    plt.xlabel("P2 - P1");
    plt.ylabel("Drag/N");
    plt.show();

    plt.title("Lift vs Pressure");
    plt.scatter(Pvals, Fyt, label="Momentum Exchange Method");
    Y, m, c = linfit.LinFit(Pvals[range(10)], Fyt[range(10)]);
    plt.plot(Pvals, m*Pvals, label="Linear Fit of Incompressible Limit");
    plt.legend();
    plt.xlabel("P2 - P1");
    plt.ylabel("Lift/N");
    plt.show();

    plt.title("Force vs Max Velocity");
    plt.plot(vmax, Fxt, label="Drag (top wall)");
    plt.plot(vmax, Fyt, label="Lift (top wall)");
    plt.plot(vmax, Fxb, label="Drag (bottom wall)");
    plt.plot(vmax, Fyb, label="Lift (bottom wall)");
    plt.legend();
    plt.show();

    plt.title("Maximum Velocity versus Pressure")
    plt.plot(Pvals, vmax);
    plt.show();
    
##    FxA = np.load('Fx_Pvar.npy');
##    FyA = np.load('Fy_Pvar.npy');
##    X_A = np.load('Pvals.npy');
##    plt.subplot(121);
##    plt.plot(X_A, FxA);
##    plt.title("Drag");
##    plt.xlabel("Pressure Gradient dP");
##    plt.ylabel("Force per Unit Time");
##    plt.subplot(122);
##    plt.plot(X_A, FyA);
##    plt.title("Lift");
##    plt.xlabel("Pressure Gradient dP");
##    plt.ylabel("Force per Unit Time");
##    plt.show();
    
    ###### Define parameters of channel flow
    P0 = 1.0;
    dP = 0.05;
    P1 = P0 + dP;
    P2 = P0 - dP;
    N = 60;
    M = 40;
    L = 40;
    eta = 25.0;
    vmax = [];

    FxT = [];
    FxB = [];
    FyT = [];
    FyB = [];

    
    Fx = [];
    Fy = [];
    
    R = 10;

    Pvals = [0.005, 0.01, 0.015,0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]; ##Note: if dP>1, there will be negative density. However, not sure if this matters, since rho_0 shouldnt matter; only the difference in rho should. 
    for dP in Pvals:
        P1 = P0 + dP;
        P2 = P0 - dP;
        gas = makeCylinder(N, M, P1, P2, L, eta, N/4, M/2, R);
        h = gas.h;
        s = gas.s;
        print "h/s =", h/s

        print "h = ", h
        print "s = ", s
        ###### Do iteration and plot
        count = 0;
        start = time.time();
        while(gas.t < 1000*s):
            gas.iterate(s);
            
##            if(count%200==0):
##                plt.subplot(211);
##                gas.plotVorticity();
##                gas.heatmap(gas.vx);
##                plot_Circle(N/4, M/2, h, R);
##
##                plt.subplot(212);
##                gas.plotFlow();
##                plot_Circle(N/4, M/2, h, R);
##                plt.pause(0.01);
####            plt.show();
##            plt.clf();
            
##            if(count%30==0):
##                gas.plotProfile();
##                plt.pause(0.01);
##                plt.clf();

##            print "rho left = ", gas.rho[-1, :];
##            print "rho right = ", gas.rho[N, :];
##            print "vx bottom = ", gas.vx[:, -1];
##            print "vx top = ", gas.vx[:, M];
##            print " " 
##            print "t = {}, error = {}, error - TOL = {}".format(gas.t, myerr, myerr - TOL);
            count+= 1;
        stop = time.time();
    
        print "the time taken to run was", stop-start;
##        plt.imshow(gas.getImage(gas.vx), cmap='hot', origin='lower', interpolation='nearest');
##        plt.show();
##        plt.plot(s*np.arange(np.size(err)), err);
##        plt.show();
##
##
##
##        plt.plot(gas.vx[16, :][range(-1, M-1)], label="vx (y profile)");
##        plt.plot(np.arange(32), (-0.05/(16*80))*np.arange(32)*((np.arange(32) - 32)), label="exact");
##        ##    plt.plot(gas.vx[:, 30], label="vx (x profile)");
##        plt.legend();
##        plt.show();  

        print "Fx (per node per time) on Top Wall: ", gas.F_top[0]*s/(gas.N*gas.t);
        print "Fy (per node per time) on Top Wall: ", gas.F_top[1]*s/(gas.N*gas.t);
##        print "Lift/Drag on Top Wall: ", gas.F_top[1]/gas.F_top[0];
        print "Fx (per node per time) on Bottom Wall: ", gas.F_bot[0]*s/(gas.N*gas.t);
        print "Fy (per node per time) on Bottom Wall: ", gas.F_bot[1]*s/(gas.N*gas.t)
##        print "Lift/Drag on Bottom Wall: ", gas.F_bot[1]/gas.F_bot[0]; 
        print "0.5*dP = ", 0.5*dP#*(h**2)*(L**2)*float(M)/N;
##        print "theor / calculated = ", (0.5*dP*(h**2)*(L**2)*float(M)/N)/(gas.F_bot[0]*s/(gas.N*gas.t));

        print " "
        FxT.append(gas.F_top[0]*s/(gas.N*gas.t));
        FyT.append(gas.F_top[1]*s/(gas.N*gas.t));
        FxB.append(gas.F_bot[0]*s/(gas.N*gas.t));
        FyB.append(gas.F_bot[1]*s/(gas.N*gas.t));
##        for count in range(np.size(gas.obstacles)):
##            print "Obstacle {} Force Calculations".format(count+1);
##            print "Fx per node per unit time: ", gas.Fx[count]*s/(np.size(Fx[count])*gas.t);
##            print "Fy per node per unit time: ", gas.Fy[count]*s/(np.size(Fy[count])*gas.t);

        
        print "v_max = ", max(gas.vx[N/2, :]);
        
        print "6 pi r v_max nu =", 6*np.pi*gas.eta*gas.cs**2*R;
        vmax.append(max(gas.vx[N/4, :]));
    WallForces = [FxT, FyT, FxB, FyB];
    np.save('WallForces_Pvar_full.npy', WallForces);
##    np.save('Fx_Pvar_topcyl.npy', Fx);
##    np.save('Fy_Pvar._topcylnpy', Fy);
    np.save('Pvals_full.npy', Pvals);
    np.save('vmax_P_full.npy', vmax);

    plt.subplot(121);
    plt.title("Force vs Pressure");
    plt.plot(Pvals, FxT, label="Drag (top wall)");
##    plt.plot(Pvals, FyT, label="Lift (top wall)");
    plt.plot(Pvals, FxB, label="Drag (bottom wall)");
##    plt.plot(Pvals, FyB, label="Lift (bottom wall)");
    plt.legend();

    plt.subplot(122);
    plt.title("Max Velocity vs Pressure");
    plt.plot(Pvals, vmax);
    plt.xlabel("pressure");
    plt.ylabel("max velocity");


    
    plt.show();

