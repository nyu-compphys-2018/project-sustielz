import numpy as np;
import matplotlib.pyplot as plt;

###NOTE: I did a time test for a couple of implementations. I.e. how much longer
###       does it take to check whether a node and its neighbors are fluid or solid
#
###Calculated time taken for 40x32 poisselle channel to iterate for a set period of time

###time taken to stream + check node type + check neighbors node type: 27.8, 22.5
###Time taken to stream and check node type: 22.5
###Time taken to stream without checking node type: 19.5
###Time taken if you store a list of which nodes are fluids: 21.9
###I believe it is best to simply use the array. It is important to have functionality for embedded objects.

class D2Q9(object): ##LBM for a D2Q9 lattice. Boundaries at lattice edges must be specified separately. 
    ##Specify parameters specific to a d2q9 lattice
    dct = [(0,0), (1,0), (0,1), (-1,0), (0,-1), (1,1), (-1,1), (-1,-1), (1,-1)]; #Returns the two values corresponding to the direction
    opp = [0, 3, 4, 1, 2, 7, 8, 5, 6]; #returns the opposite direction
    w = [4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0]; #weighting function
    cs = 1.0/np.sqrt(3); ##Speed of sound in lattice units

    ##Default system parameters
    L = 1.0;
    tau = 1.0;
    
    t = 0.0;
    
    f = [];
    feq = [];
    fs = [];
    
    Fx = [];     #Force on immersed boundaries
    Fy = [];
    
    nodes = []; ##Stores a list of all points which are fluids. Excludes walls.
    solidlist = []; ##Stores a list of all points which are solids. Excludes walls. 

    


######## Initialization ########
    def __init__(self, N, M):              ##Initialize an NxM lattice, consisting of fluid with walls at edges.
        self.N = N;
        self.M = M;
        self.rho = np.zeros([N+2, M+2]);   ##walls at -1, N, -1, M
        self.vx = np.zeros([N+2, M+2]);
        self.vy = np.zeros([N+2, M+2]);
        for I in range(9):
            self.f.append(np.zeros([N+2, M+2]));
            self.feq.append(np.zeros([N+2, M+2]));
            self.fs.append(np.zeros([N+2, M+2]));

        self.fluid = np.ones([N+2, M+2]); ##Initialize the fluid mask ==1 if fluid, ==0 if solid, ==2 if wall
        self.defineSolid(self.fluid);     ##Properly initialize fluid mask, and generate lists of solid nodes and fluid nodes

        self.RIGHT = [];                 ##Initialize four walls. (Note: walls exclude corners)
        for j in range(M):
            self.RIGHT.append([N, j]);
        self.UP = [];
        for i in range(N):
            self.UP.append([i, M]);
        self.LEFT = [];                 ##Initialize four walls. (Note: walls exclude corners)
        for j in range(M):
            self.LEFT.append([-1, j]);
        self.DOWN = [];  
        for i in range(N):
            self.DOWN.append([i, -1]);
            
        self.C5, self.C6, self.C7, self.C8 = [N, M], [-1, M], [-1, -1], [N, -1] ##Init

        
######## Methods for defining where immersed solids are ########
    def defineSolid(self, MAP):       ##Given a mask, update self.fluid, and update list of solid nodes and list of fluid nodes
        N, M = self.N, self.M;
        MAP[-1, :] = 2*np.ones(M+2);      ##Make sure the walls are always walls
        MAP[N, :] = 2*np.ones(M+2);
        MAP[:, -1] = 2*np.ones(N+2);
        MAP[:, M] = 2*np.ones(N+2);
        
        self.fluid = MAP;
        
        self.nodes = [];
        self.solidlist = [];
        for i in range(self.N):
            for j in range(self.M):
                if(MAP[i, j] == 1):
                    self.nodes.append([i, j]);
                else:
                    self.solidlist.append([i, j]);

    def plotMask(self): #Plot the mask
        plt.imshow(self.getImage(self.fluid), origin='lower');
        plt.colorbar();
        plt.show();

######## Methods to update macroscopic quantities and f_eq ########
    def updateRho(self):
        self.rho = self.rho*0.0;
        for I in range(9):
            self.rho += self.f[I];
##        print self.rho

    def updateV(self):
        self.vx = self.vx*0;
        self.vy = self.vy*0;
        for I in range(9):
            ex, ey = self.dct[I];
            self.vx += (1.0/self.rho)*self.f[I]*ex;
            self.vy += (1.0/self.rho)*self.f[I]*ey;

    def updateEq(self, rho, vx, vy):
        for I in range(9):
            ex, ey = self.dct[I];
            A = ex*vx + ey*vy;
            B = vx*vx + vy*vy;
            self.feq[I] = self.w[I]*rho*( 1 + 3.0*A + 4.5*A**2 - 1.5*B );

######## Collision and Streaming ########

    def Solid_Bounce(self):
        for P in self.solidlist: 
            i, j = P;
            for I in [1, 2, 5, 6]:
                J = self.opp[I];
                ex, ey = self.dct[I];
                front = self.fluid[i+ex, j+ey];  ##We check the type of node in front and behind to determine what type of boundary 
                back = self.fluid[i-ex, j-ey];
                if(back == 2 and front == 2):
                    pass;
                elif(front == 1 and back == 0):
                    self.f[I][i+ex, j+ey] = self.fs[I][i, j];
                    self.f[I][i, j] = self.fs[J][i, j];
                    self.Fx -= self.fs[J][i,j]*ex;
                    self.Fy -= self.fs[J][i,j]*ey;
##                    print "Fx += ", self.f[J]*ex
                    
                elif(front == 0 and back == 1):
                    self.f[J][i-ex, j-ey] = self.fs[J][i, j];
                    self.f[J][i, j] = self.fs[I][i, j];

                    self.Fx += self.fs[I][i,j]*ex;
                    self.Fy += self.fs[I][i,j]*ey;
##                    print "Fx += ", self.f[I]*ex
                elif(front == 1 and back == 1):
                    self.f[I][i+ex, j+ey] = self.fs[I][i, j];
                    self.f[J][i-ex, j-ey] = self.fs[J][i, j];
                    
                    self.Fx += (self.fs[I][i,j]-self.fs[J][i,j])*ex;
                    self.Fy += (self.fs[I][i,j]-self.fs[J][i,j])*ey;
    

                    temp = self.fs[J][i, j];
                    self.f[J][i, j] = self.fs[I][i, j];
                    self.f[I][i, j] = temp;
##                    print "Fx += ", self.f[I]*ex
##                    print "i, j =", i, j;

##    def bounce(self):#, SHAPE): ##Half Bounceback. Take place after collision, before general streaming.
##        SHAPE = self.solidlist;                 ##Updates f and fs, and streams boundary nodes
##        fx = 0.0;
##        fy = 0.0;
##        for P in SHAPE: 
##            i, j = P;
##            for I in range(9):
##                J = self.opp[I];
##                ex, ey = self.dct[I];################Refine this part
####                if(i+ex>=self.N+1 or i+ex<=-2 or j+ey>=self.M+1 or j+ey<=-2):
####                    pass;
##                if(self.fluid[i+ex, j+ey] == 1):    ##Bounce only into fluid nodes
####                    if(self.fluid[i-ex, j-ey] ==1): ##Special case: both directions bounce back into fluid
####                        if(I == 1 or I == 2 or I == 5 or I == 6): ##Note: if I = 3, 4, 7 or 8, then we already performed this bounce during the opposite direction
####                            
######                    else:   ##
##                    J = self.opp[I];
##                    self.f[I][i+ex, j+ey] = self.fs[I][i, j];
##                    self.f[I][i, j] = self.fs[J][i, j];
##
##                    
####                    self.f[J][i, j] = self.fs[I][i, j];
####                    self.f[J][i+ex, j+ey] = self.fs[I][i, j]; ##Stream post-collision boundary node back into fluid 
##
##                    fx += ex*(self.f[J][i+ex, j+ey] - self.fs[I][i, j]);
##                    fy += ey*(self.f[J][i+ex, j+ey] - self.fs[I][i, j]);
####                    print "fs[{}][{} {}] - > f[{}][{} {}]".format(I, i, j, J, i, j); 
##        return fx, fy;
    


    def stream(self, I, DOMAIN): ##Stream all points (i, j) in (xvals, yvals) in the direction I
        ex, ey = self.dct[I];
        for P in DOMAIN:
            [i, j] = P;
            self.f[I][i+ex, j+ey] = self.fs[I][i, j];
           
##    def Stream(self):
##        N, M = self.N, self.M;
##        for I in range(9):      ##Stream fluid nodes in all directions
##            self.stream(I, self.nodes);
##        ##Stream all boundaries (including bounceback walls)
##        
##       
##        for I in [3, 6, 7, 2, 4]:  ##Stream right wall
##            self.stream(I, self.RIGHT);
##        for I in [4, 7, 8, 1, 3]:  ##Stream top wall
##            self.stream(I, self.UP);
##        for I in [1, 5, 8, 2, 4]:  ##Stream left wall
##            self.stream(I, self.LEFT);
##        for I in [2, 5, 6, 1, 3]:  ##Stream bottom wall
##            self.stream(I, self.DOWN);
##
##  
##
##        for I in [7, 4, 3]:
##            self.stream(I, [self.C5]) ###Stream top right corner
##        for I in [8, 4, 1]:
##            self.stream(I, [self.C6]) ###Stream top right corner
##        for I in [5, 2, 1]:
##            self.stream(I, [self.C7]) ###Stream top right corner
##        for I in [6, 2, 3]:
##            self.stream(I, [self.C8]) ###Stream top right corner
##
##        self.Solid_Bounce();

    def stream2(self, I, i, j): ##Stream all points (i, j) in (xvals, yvals) in the direction I
        ex, ey = self.dct[I];
        self.f[I][i+ex, j+ey] = self.fs[I][i, j];

    def stream_Walls(self):
        N, M = self.N, self.M;
        for I in [2, 5, 6, 1, 3]:  ##Stream bottom wall
            self.stream(I, self.DOWN);

        for I in [4, 7, 8, 1, 3]:  ##Stream top wall
            self.stream(I, self.UP);

        for I in [1, 5, 8, 2, 4]:  ##Stream left wall
            self.stream(I, self.LEFT);
            
        for I in [3, 6, 7, 2, 4]:  ##Stream right wall
            self.stream(I, self.RIGHT);



##        self.stream3(2, range(N), [-1]);   #Stream top and bottom walls
##        self.stream3(5, range(N), [-1]);
##        self.stream3(6, range(N), [-1]);
##        self.stream3(1, range(N), [-1]);
##        self.stream3(3, range(N), [-1]);
            
##        self.stream3(4, range(N), [M]);
##        self.stream3(7, range(N), [M]);
##        self.stream3(8, range(N), [M]);
##        self.stream3(1, range(N), [M]);
##        self.stream3(3, range(N), [M]);

        
##        self.stream3(1, [-1], range(M));    #Stream side walls
##        self.stream3(5, [-1], range(M));
##        self.stream3(8, [-1], range(M));
##        self.stream3(2, [-1], range(M));
##        self.stream3(4, [-1], range(M));


##        self.stream3(3, [N], range(M));
##        self.stream3(6, [N], range(M));
##        self.stream3(7, [N], range(M));
##        self.stream3(2, [N], range(M));
##        self.stream3(4, [N], range(M));

      
        ##Finally, stream corners

##        self.f[5][0, 0] = self.fs[5][-1, -1];
##        self.f[2][-1, 0] = self.fs[2][-1, -1];
##        self.f[1][0, -1] = self.fs[1][-1, -1];
        self.stream2(5, -1, -1);
        self.stream2(2, -1, -1);
        self.stream2(1, -1, -1);

##        self.f[8][0, M-1] = self.f[8][-1, M];
##        self.f[4][-1, M-1] = self.f[4][-1, M];
##        self.f[1][0, M] = self.f[1][-1, M];
        self.stream2(8, -1, M);
        self.stream2(4, -1, M);
        self.stream2(1, -1, M);

        self.stream2(6, N, -1);
        self.stream2(2, N, -1);
        self.stream2(3, N, -1);

        self.stream2(7, N, M);
        self.stream2(4, N, M);
        self.stream2(3, N, M);

    def Stream(self): ##Stream all fluid and wall nodes
        N, M = self.N, self.M;
        self.f[0] = 1.0*self.fs[0]; #f0 always streams to itself
        for P in self.nodes: ###### Stream all over fluid nodes
            i, j = P;
            for I in range(9):
                ex, ey = self.dct[I];
                self.f[I][i+ex, j+ey] = self.fs[I][i, j];   ## Stream fluid nodes
        self.stream_Walls();    ## Next, stream the walls.
        self.Solid_Bounce();    ##### Stream, bounce back + force eval on immersed boundaries




    def Collision(self, s):
        N, M = self.N, self.M;
        for I in range(9):  ##Perform Collision on Everything  
            self.fs[I] = self.f[I] + (1/self.tau)*(self.feq[I] - self.f[I]);

    def boundary(self): ##In general, the boundary will depend on the system
        pass;
    
    def iterate(self, s):
        self.updateRho();
        self.updateV();
        self.updateEq(self.rho, self.vx, self.vy);
        self.Collision(s);
        self.Stream();
        self.boundary();
        
        self.t += s;

    def getImage(self, prop):
        N, M = self.N, self.M;
        im = np.zeros([N+2, M+2]);
        for i in range(0, N+2):
            for j in range(0, M+2):
                im[i, j] = prop[i-1,j-1];

        return np.transpose(im);                
                       
    def plotF(self):
        pos = [5, 6, 2, 4, 8, 3, 1, 7, 9]; 
        for k in range(3):
            for l in range(3):
                n = 3*k + l;
                
                plt.subplot(3, 3, pos[n]);
                plt.title("f_{}".format(n));
                plt.imshow(self.getImage(self.f[n]), cmap='hot', origin='lower', interpolation='nearest');
        plt.show();

    def plotFs(self):
        pos = [5, 6, 2, 4, 8, 3, 1, 7, 9];
        for k in range(3):
            for l in range(3):
                n = 3*k + l;
                plt.subplot(3, 3, pos[n]);
                plt.title("f_{}".format(n));
                plt.imshow(self.getImage(self.fs[n]), cmap='hot',  origin='lower', interpolation='nearest');
        plt.show();    
    
if __name__ == '__main__':
    ##    test = LBM(5, 5);
##    for I in range(9):
##        for i in range(7):
##            for j in range(7):
##                test.f[I][i-1, j-1] = 7*i + j;
##                if(i<5 and j<5):
##                    test.fs[I][i, j] = 50*(j+1) + 5*i ;
##                else:
##                    test.fs[I][i, j] = 8
##                test.f[I] = 1.0*test.fs[I];
##        
##    s = 1*test.h;
##    test.plotF();
##    while(test.t<3*s):
##        test.iterate(s);
##        test.plotFs();
##        test.plotF();
####    N = 7;
####    M = 7;
####    test = np.ones([N, M]);
####    solid = np.ones([N, M]);
####    for i in range(2):
####        for j in range(2):
####            solid[i+3,j+3] = 0;
####    test = test*solid;
####    plt.imshow(test);
####    plt.show();
    print "You ran the wrong file, but you 're still doing a great job bud"


