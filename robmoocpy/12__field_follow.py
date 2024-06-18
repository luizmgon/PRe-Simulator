from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py


    
    
def draw(x):
    draw_tank(x,'darkblue',0.3)
    a,b = array([[-30],[0]]), array([[30],[0]])
    draw_segment(a,b,'red',2)
    
def f(x,u):
    θ=x[2,0]
    return array([[cos(θ)], [sin(θ)],[u]])
           

x=array([[-2],[-2],[3]])
dt= 0.05
s=10

def f1(x1,x2):        
        return x2,-x1
 
ax=init_figure(-s,s,-s,s)
draw_field(ax,f1,-s,s,-s,s,1)


for t in arange(0,8,dt):
    draw(x)
    pause(0.1)
    u=-0.3
    x=x+dt*f(x,u)

