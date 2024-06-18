from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py
from numpy import *


def f(x,u):
    x,u=x.flatten(),u.flatten()
    xdot = array([[x[3]*cos(x[2])],[x[3]*sin(x[2])],[u[0]],[u[1]]])
    return(xdot)

def control(x,w,dw,ddw):
    u=array([[0],[0]]) #TO DO
    return u    
    

ax=init_figure(-50,50,-50,50)
m   = 20
X   = 10*randn(4,m)
a,dt = 0.1,0.1

for t in arange(0,3,dt):
    clear(ax)
    for i in range(m):        
        w = zeros((2,1)) #TO DO
        dw = zeros((2,1))  #TO DO
        ddw = zeros((2,1))#TO DO
        x=X[:,i].reshape(4,1)
        u       = control(x,w,dw,ddw)
        x=X[:,i].reshape(4,1)
        draw_tank(x,'b')
        x=x+f(x,u)*dt        
        X[:,i]  = x.flatten()
        plot([w[0][0]],[w[1][0]],'r+')


