from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py


def f(x,u):  # x,y,z,ψ,vx,vy,vz,w
    x1,x2,x3,ψ,vx,vy,vz,ω=x[0:8,0]
    u1,u2,u3=u[0:3,0]
    return np.array([[vx],[vy],[vz],[ω],[u1*np.cos(ψ)],[u1*np.sin(ψ)],[u3],[u2]])

def draw(x,w1,w2,w3):
    x1,x2,x3,ψ,vx,vy,vz,ω=x[0:8,0]
    u1,u2,u3=u[0:3,0]
    clean3D(ax,-30,30,-30,30,0,60)
    draw_axis3D(ax,0,0,0,np.eye(3,3),10)
    R=eulermat(0,0,ψ)
    M=tran3H(x1,x2,x3)@eulerH(0,0,ψ) @ auv3H()
    draw3H(ax,M,'blue',True)
    U1=5*R@np.array([[u1],[0],[0]])
    draw_arrow3D(ax,x1,x2,x3,*U1[0:3,0],"red")
    U2=5*R@np.array([[0],[0],[u2]])
    draw_arrow3D(ax,x1,x2,x3,*U2[0:3,0],"green")
    U3=5*R@np.array([[0],[0],[u3]])
    draw_arrow3D(ax,x1,x2,x3,*U3[0:3,0],"red")
    ax.scatter(w1,w2,w3,color='magenta')  #target
    pause(0.001)

def w(t):
    return 30*sin(0.1*t),10*cos(0.1*t),10*(1+cos(0.3*t))

def build_ρ():
    x1,x2,x3,ψ,vx,vy,vz,ω,c1,b1=symbols("x1 x2 x3 ψ vx vy vz ω c1 b1")
    C=CoordSystem('C',Patch('P',Manifold('M',10)),[x1,x2,x3,ψ,vx,vy,vz,ω,c1,b1])


ax = figure3D()
x = np.array([[0],[0],[0],[0],[1],[0],[0],[0]])
dt = 0.02
c1,b1=0.1,0.1    
for t in arange(0,3,dt):
    x1,x2,x3,ψ,vx,vy,vz,ω=x[0:8,0]
    u=np.array([[sin(t)],[cos(t)],[sin(t)]])    
    x = x + dt * f(x,u)
    draw(x,1,2,3)

