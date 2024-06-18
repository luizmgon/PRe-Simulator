from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

def f(x1,x2):        
    return -(x1**3+x2**2*x1-x1+x2),-(x2**3+x1**2*x2-x1-x2)


xmin,xmax,ymin,ymax=-2.5,2.5,-2.5,2.5 
ax=init_figure(xmin,xmax,ymin,ymax)
draw_field(ax,f,xmin,xmax,ymin,ymax,0.3)    
dt=0.05
x=array([[0],[1]])
for t in arange(0,10,dt):
    x1,x2=x[0,0],x[1,0]
    dx1,dx2=f(x1,x2)
    x=x+dt*array([[dx1],[dx2]])
    ax.scatter(x1,x2,1.6,color='red')

pause(1)
    
    






