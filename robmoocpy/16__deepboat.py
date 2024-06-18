import tensorflow as tf
from roblib import *

def draw_boat(X,t,col='darkblue'):
    if (t-t//1>0.01): return()
    draw_tank(X,col,0.25,0.3)
    plot([-30,30], [0,0], 'red')
    pause(0.001)

def control(X):
    x,y,θ,v=tolist(X)
    u=array([[tanh(1-v)],[-0.1]])
    return u


#**************************************************************************************************

def SimuDubins():
    def f_dubins(X,u):
        u1,u2=tolist(u)
        θ,v=tolist(X[2:4])
        dX=array([[v*cos(θ)],[v*sin(θ)],[u1-u2],[u1+u2]])
        return dX
    X=array([[-30],[-10],[-1],[0]])  #x,y,theta,v
    for t in arange(0,10,dt):
        draw_boat(X,t)
        u=control(X) #x,y,θ,v
        X=X+dt*f_dubins(X,u)

dt=0.1
ax=init_figure(-35,20,-20,20)

SimuDubins()

pause(10)
