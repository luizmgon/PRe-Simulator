from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

def f(x,u):
    x    = x.flatten()
    return array([[5*cos(x[2])],[5*sin(x[2])],[u]])


x    = array([[15],[20],[1]])
dt   = 0.1
ax=init_figure(-30,30,-30,30)
for t in arange(0,1,dt):
    clear(ax)
    draw_disk(ax,array([[0],[0]]),10,"cyan")
    u = 0.5
    draw_tank(x,'red')
    x = x+dt*f(x,u)            

