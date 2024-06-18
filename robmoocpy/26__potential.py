from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

def f(x,u):
    x,u  = x.flatten(), u.flatten()
    v,θ = x[2],x[3]    
    return array([[v*cos(θ)],[v*sin(θ)],[u[0]],[u[1]]])
    

def f1(x1,x2):  
    return -x1,-x2
    
    
x    = array([[4,-3,1,2]]).T #x,y,v,θ
dt   = 0.1
s=5
ax=init_figure(-s,s,-s,s)
for t in arange(0,5,dt):
    clear(ax)
    phat = array([[1],[2]])
    qhat = array([[3],[4]])        
    draw_disk(ax,qhat,0.3,"magenta")
    draw_disk(ax,phat,0.2,"green")
    u=array([[0],[0.3]])
    x=x+dt*f(x,u)    
    draw_tank(x[[0,1,3]],'red',0.2) # x,y,θ
    draw_field(ax,f1,-s,s,-s,s,0.4)

pause(1)    


