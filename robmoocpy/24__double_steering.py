from roblib import *
import sympy as sp

def draw_robot(x):
   p1,p2,ψ,s1,s2,s3=list(x[0:6,0])
   M0=array([[ -1  ,1],[0,0]])
   M0=add1(M0)
   W=array([[-0.5,0.5],[0,0],[1,1]])
   R1=tran2H(p1,p2)@rot2H(ψ)
   M1=R1@M0
   δ1=arctan2(s1,s3)
   δ2=arctan2(s2,s3)
   W1=R1@tran2H(1,0)@rot2H(δ1)@W
   W2=R1@tran2H(-1,0)@rot2H(δ2)@W
   plot2D(M1,'blue',1)
   plot2D(W1,'green',1)
   plot2D(W2,'black',1)

def Aψ(ψ):
    return array([[-0.5*sin(ψ),-0.5*sin(ψ), cos(ψ)], [0.5*cos(ψ),0.5*cos(ψ),sin(ψ)],[1,-1,0]])

def f(x,u):
    u1,u2,u3,u4=list(u[0:4,0])
    ψ,s1,s2,s3=list(x[2:7,0])
    v1=sqrt(s1**2+s3**2)
    v2=sqrt(s2**2+s3**2)
    ds = array([[(s1/v1)*u1+s3*u2], [(s2/v2)*u3+s3*u4],[(s3/v1)*u1-s1*u2]])
    s=array([[s1],[s2],[s3]])
    return(vstack((Aψ(ψ)@s,ds)))


a=array([4,1])
t=sp.symbols('t')
x,y=3*sp.cos(t),3*sp.sin(2*t)
dxd = sp.lambdify(t,sp.diff(x,t))
dyd = sp.lambdify(t,sp.diff(y,t))


ax = init_figure(-5,5,-5,5)
x=array([[3],[-4],[0],[1],[1],[1]]) #x,y,ψ,s1,s2,s3
dt=0.01
for t in arange(0,1,dt):
    clear(ax)
    draw_robot(x)
    draw_disk(ax,a,0.1,"blue")
    T = arange(0,2*pi,0.01)
    plot(3*cos(T), 3*sin(2*T),color='magenta')
    pause(0.001)
    u=array([[0],[0],[0],[0]])
    x = x + dt*f(x,u)
pause(1)

