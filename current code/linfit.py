from numpy import loadtxt, linspace
from matplotlib.pyplot import plot, xlabel, ylabel, title, show

def LinFit(X, Y):
	#Calculate relevant quantities
	N = len(X)
	Ex = sum(X)/N
	Ey = sum(Y)/N
	Exx = sum(X*X)/N
	Exy = sum(X*Y)/N

	#Calculate best fit line
	m = (Exy -Ex*Ey)/(Exx -Ex**2)
	c = (Exx*Ey -Ex*Exy)/(Exx - Ex**2)

	Y_fit = m*X + c
	return Y_fit, m, c;

