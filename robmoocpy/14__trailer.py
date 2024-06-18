from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py
from sympy import *
from sympy.diffgeom import *


def L(F,g,i=1):
    if size(g)==2: return Matrix([[L(F,g[0],i)],[L(F,g[1],i)]])
    if i==1: return LieDerivative(F,g)
    return L(F,L(F,g,i-1))
    
def ψ(x1,x2): return x2,-(1*(x1**2)-1)*x2-x1    

x1,x2=symbols("x1 x2")
C=CoordSystem('C',Patch('P',Manifold('M',2)),[x1,x2])
x1,x2 = C.coord_functions()
E = C.base_vectors() 
Fx = x2*E[0]-((1*(x1**2)-1)*x2+x1)*E[1]
Hx = x1-cos(x2)
print(Hx)
print(L(Fx,Hx))
hx = lambdify((x1,x2), Hx)
print(hx(1,2))


x = np.array([[0],[0],[0],[1],[1]])
sc=3
ax=init_figure(-sc,sc,-sc,sc)
x1,x2,x3,x4,x5 = x[0:5,0]    
draw_field(ax,ψ,-sc,sc,-sc,sc,0.3)
draw_tank_trailer(x1,x2,x3,x4,x5); 
pause(1)


