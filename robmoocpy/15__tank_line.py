from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

def f(X,u):
    θ=X[2,0]
    return array([[cos(θ)], [sin(θ)],[u]])

X=array([[-20],[-10],[4]])
u=1
dt= 0.1
a,b = array([[-30],[-4]]), array([[30],[6]])
ax=init_figure(-40,40,-40,40)

for t in arange(0,5,dt):
    clear(ax)
    draw_tank(X,'darkblue')
    plot2D(hstack((a,b)),'red')
    plot2D(a,'ro')
    plot2D(b,'ro')    
    X   = X+dt*f(X,u)

    
    