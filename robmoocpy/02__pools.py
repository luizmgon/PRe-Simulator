from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

def draw_pools(x):
    x=x.flatten()
    plot([0,0],[10,1],'black',linewidth=2)    
    plot([-7,23],[0,0],'black',linewidth=5)    
    plot([16,16],[1,10],'black',linewidth=2)    
    plot([4,4,6,6],[10,1,1,10],'black',linewidth=2)    
    plot([10,10,12,12],[10,1,1,10],'black',linewidth=2)    
    P=array([[0,x[0]],[0,1],[-6,0],[22,0],[16,1],[16,x[2]],[12,x[2]],[12,1]
            ,[10,1],[10,x[1]],[6,x[1]],[6,1],[4,1],[4,x[0]]])
    draw_polygon(ax,P,'blue')
    P=array([[1,10],[1,x[0]],[1+0.1*u[0,0],x[0]],[1+0.1*u[0,0],10]])
    draw_polygon(ax,P,'blue')
    P=array([[13,10],[13,x[2]],[13+0.1*u[1,0],x[2]],[13+0.1*u[1,0],10]])
    draw_polygon(ax,P,'blue')


def f(x,u):
    return(array([[0.1],[0.3],[1.5]]))


dt = 0.05
x = array([[4],[5],[2]])
u = array([[1],[2]])
ax=init_figure(-10,25,-2,12)

for t in arange(0,1,dt) :
    clear(ax)
    draw_pools(x)
    x = x + dt*f(x,u)  
pause(3)
