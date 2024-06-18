from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

def f(x,u):
    x=x.flatten()
    θ = x[2]
    return array([[cos(θ)], [sin(θ)], [u]])

ax=init_figure(-60,60,-60,60)
m   = 20
X   = 20*randn(3,m)
dt  = 0.2

for t in arange(0,10,dt):
    clear(ax)
    for i in range(m):
        xi=X[:,i].flatten()
        xi=xi.reshape(3,1)
        draw_tank(xi,'b')
        u=0
        xi=xi+f(xi,u)*dt        
        X[:,i]  = xi.flatten()        



