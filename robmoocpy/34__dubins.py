from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

def f(x,u):
    x = x.flatten()
    θ = x[2]
    return array([[cos(θ)],[sin(θ)],[u]])

def control(x):
    u=0
    return u
    
    
x   = array([[0],[0],[0.1]])
dt  = 0.1
ax=init_figure(-10,10,-10,10)

for t in arange(0,3,dt):
    clear(ax)
    u = control(x)
    x = x + dt*f(x,u)    
    draw_tank(x,'red',0.3) 
    
