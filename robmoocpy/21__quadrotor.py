from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py
ax=figure3D()

m,g,b,d,l=10,9.81,2,1,1
I=array([[10,0,0],[0,10,0],[0,0,20]])
dt = 0.01  
B=array([[b,b,b,b],[-b*l,0,b*l,0],[0,-b*l,0,b*l],[-d,d,-d,d]])


def clock_quadri(p,R,vr,wr,w):
    w2=w*abs(w)
    τ=B@w2.flatten()
    p=p+dt*R@vr
    vr=vr+dt*(-adjoint(wr)@vr+inv(R)@array([[0],[0],[g]])+array([[0],[0],[-τ[0]/m]]))
    R=R@expw(dt*wr)
    wr=wr+dt*(inv(I)@(-adjoint(wr)@I@wr+τ[1:4].reshape(3,1)))
    return p,R,vr,wr

    
def control(p,R,vr,wr):
    return array([[6],[5],[5],[5]]) 


p = array([[0], [0], [-5]])  #x,y,z (front,right,down)
R = eye(3)
vr = array([[10], [0], [0]])
wr = array([[0], [0], [0]])
α=array([[0,0,0,0]]).T #angles for the blades

for t in arange(0,1,dt):
    w=control(p, R, vr, wr)
    p, R, vr, wr = clock_quadri(p, R, vr, wr, w)
    clean3D(ax, -25, 25, -25, 25, 0, 25)
    draw_quadrotor3D(ax, p, R, α, 5 * l)
    α = α + dt * 30 * w
    pause(0.001)
pause(1)


