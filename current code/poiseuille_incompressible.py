##############class poiseuille_incompressible(poiseuille):
##############
##############    def updateV(self):
##############        self.vx = self.vx*0;
##############        self.vy = self.vy*0;
##############        for I in range(9):
##############            ex, ey = self.dct[I];
##############            self.vx += self.f[I]*ex;
##############            self.vy += self.f[I]*ey;
##############
##############    def updateEq(self, rho, vx, vy):
##############        for I in range(9):
##############            
##############            ex, ey = self.dct[I];
##############            A = ex*vx + ey*vy;
##############            B = vx*vx + vy*vy;
##############            self.feq[I] = self.w[I]*(rho + 3.0*A + 4.5*A**2 - 1.5*B );
##############
##############    def boundary(self): 
##############        N, M = self.N, self.M;
################        self.bounce(self.boundaries[1]);
################        self.bounce(self.boundaries[3]);
##############            
##############        ##Bounce on top and bottom walls (i.e. no-slip)
##############        for I in [2, 5, 6]:
##############            ex, ey = self.dct[I];
##############            J = self.opp[I];
##############            for i in range(-1, N+1):
##############                if(self.full):
##############                    self.f[I][i, -1] = self.fs[J][i, -1];
##############                    self.f[J][i, M] = self.fs[I][i, M];
##############                else:
##############                    print "hey dog"
##############                    self.f[I][i, -1] = self.fs[J][i+ex, M+ey];
##############                    self.f[J][i, M] = self.fs[I][i-ex, M-ey];
##############        ######Zou-He boundary conditions to simulate pressure gradient at sides (copied from zou-he paper)
##############        rhoL = self.P1/self.cs**2;
##############        rhoR = self.P2/self.cs**2;
##############        vyL = 0.0;
##############        vyR = 0.0;
##################        fL = [];
##################        fR = [];
##################        for I in range(9):
##################            fL.append(self.F[I].fs[-1, :]);
##################            fR.append(self.F[I].fs[N, :]);
##############        fL = np.zeros(9);
##############        fR = np.zeros(9);
##############        for j in range(-1, M+1):
##############            for I in range(9):
##############                fL[I] = self.f[I][-1, j]
##############                fR[I] = self.f[I][N, j]
##############            
##############            ##Left wall
##############            vxL = rhoL - (fL[0] + fL[2] + fL[4] + 2.0*(fL[3] + fL[6] + fL[7]));
##############            self.f[1][-1, j] = fL[3] + (2.0/3.0)*vxL    
##############            self.f[5][-1, j] = fL[7] - 0.5*(fL[2] - fL[4]) + (1.0/6.0)*vxL + 0.5*vyL;
##############            self.f[8][-1, j] = fL[6] + 0.5*(fL[2] - fL[4]) + (1.0/6.0)*vxL - 0.5*vyL;
##############
##############            ##Right wall
##############            vxR = -rhoR + (fR[0] + fR[2] + fR[4] + 2.0*(fR[1] + fR[5] + fR[8])) ;
##############            self.f[3][N, j] = fR[1] - (2.0/3.0)*vxR    
##############            self.f[7][N, j] = fR[5] + 0.5*(fR[2] - fR[4]) - (1.0/6.0)*vxR - 0.5*vyR;
##############            self.f[6][N, j] = fR[8] - 0.5*(fR[2] - fR[4]) - (1.0/6.0)*vxR + 0.5*vyR;
##############
##############
