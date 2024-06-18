from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

def draw_rob(x,col):
    x = x.flatten()
    lx,ly,ψ = x[0],x[1],x[2]
    R = Rlatlong(lx,ly) @ eulermat(0,0,ψ)
    draw_robot3D(ax,latlong2cart(ρ,lx,ly),R,col,1) 
    
def f(x,u):
    x = x.flatten()
    lx,ly,ψ = x[0],x[1],x[2]
    return array([[cos(ψ)/(ρ*cos(ly))], [sin(ψ)/ρ], [u]])


ρ = 30 
ax=figure3D()
x   = array([[-2],[0],[0.3]])
dt = 0.1

for t in arange(0,0.1*10,dt):
    clean3D(ax,-ρ,ρ,-ρ,ρ,-ρ,ρ)
    draw_earth3D(ax,ρ,eye(3),'gray')    
    u = 0.1 * randn(1) [0]
    x = x + dt*f(x,u)    
    draw_rob(x,"blue")
    #draw_earth3D(ax,ρ,eye(3),'gray')
    pause(0.001)

pause(1)
