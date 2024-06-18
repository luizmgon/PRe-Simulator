from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py


def f(x,u):
    x=x.flatten()
    u=u.flatten()
    return array([[x[3]*cos(x[2])],[x[3]*sin(x[2])],[u[0]],[u[1]]])
        
def control (x,w,dw):
    x=x.flatten()
    A = array([[-x[3]*sin(x[2]), cos(x[2])],
                  [x[3]*cos(x[2]), sin(x[2])]])
    y = array([[x[0]],[x[1]]])
    dy = array([[x[3]*cos(x[2])],[x[3]*sin(x[2])]])
    return inv(A) @ ((w - y) + 2*(dw - dy))

    
    
def setpoint(t): 
    w = array([[5+cos(t)],[5+sin(t)]])
    return w

def dsetpoint(t): 
    dw = array([[-sin(t)],[cos(t)]])
    return dw

    
ax=init_figure(-1,11,-1,11)

P = array([[1,1,1,1,2, 3,4,5,4,8,10,8],
           [1,4,7,9,10,8,6,4,0,0,0,8]])
n = len(P[0])-1
plot(P[0], P[1], 'or')
dt = 0.1
k=0


A1=array([[2,0],[4,2],[2,7]])
A2=array([[7,2],[8,3],[3,10]])
draw_polygon(ax,A1,'green')
draw_polygon(ax,A2,'green')

x = array([[0,0,0,1]]).T
for t in arange(0,10,dt):
    w = setpoint(t)
    dw = dsetpoint(t)
    u = control(x,w,dw)
    plot(w[0],w[1], 'm.')
    x = x+f(x,u)*dt
    draw_tank(x,'darkblue',0.2)
pause(10)
