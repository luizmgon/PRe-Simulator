
from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

def f(x,u):
    x,u=x.flatten(),u.flatten()
    xdot = array([[x[3]*cos(x[2])],[x[3]*sin(x[2])],
                  [u[0]],[u[1]],[x[3]]])
    return(xdot)

def control(x,w,dw):
    u=array([[0],[0]]) #TO DO
    return u    
    
    
ax=init_figure(-30,30,-30,30)
xa  = array([[10],[0],[1],[1],[0]])
m= 6
X=array([4*arange(0,m),zeros(m),ones(m),3*ones(m),zeros(m)])
Lx,Ly = 20,5
e   = np.linspace(0.,2*pi,30)
p   = array([[Lx*cos(e)],[Ly*sin(e)]])
S   = zeros((5,1))
dt  = 0.05
for t in arange(0,5,dt):
    clear(ax)
    wa  = array([[0],[0]]) #TODO
    dwa = array([[0],[0]]) #TODO
    ua  = control(xa,wa,dwa)    
    plot(wa[0][0],wa[1][0],'ro')
    plot(p[0][0],p[1][0])
    draw_tank(xa,'blue')
    xa  = xa + dt*f(xa,ua)
    for i in range(m):
        ui      = array([[0],[0]])  #TODO
        x=X[:,i].reshape(5,1)
        draw_tank(x,'black')
        x=x+f(x,ui)*dt        
        X[:,i]  = x.flatten()            
pause(1)


