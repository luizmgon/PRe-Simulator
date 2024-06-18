from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

def f(x,u):
    xr,yr,θr,vr=x.flatten()
    u1,u2=u.flatten()
    return (array([[vr*cos(θr)],[vr*sin(θr)],[u1],[u2]]))



ax=init_figure(-30,30,-30,30)

dt = 0.1
x = array([[0],[1],[pi/3],[1]])
u = array([[1],[1]])


for t in arange(0,10,dt) :
    clear(ax)
    draw_tank(x)  	
    x = x+dt*f(x,u)
pause(1)