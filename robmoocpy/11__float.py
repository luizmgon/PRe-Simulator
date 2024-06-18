from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

def draw_buoy(x):
    clear(ax) 
    x=x.flatten()
    plot([-10,10],[0,0],'black',linewidth=1)    
    d=x[0]
    P=array([[-ech,-1.8*ech],[ech,-1.8*ech],[ech,0],[-ech,0]])
    draw_polygon(ax,P,'blue')
    plot([   0,   L,  L,  L/2,   L/2,   L/2,  0,  0],
         [-L-d,-L-d, -d,   -d,   2-d,    -d, -d,-L-d],'black',linewidth=3)
    b=-x[2]     
    P=array([[0,-L-d+L],[L,-L-d+L],[L,-L/2-L*b/2-d],[0,-L/2-L*b/2-d]])
    draw_polygon(ax,P,'white')
    
        
ech=5
x = array([[3],[0],[0]])
L=1 #length of the cube
ax=init_figure(-ech,ech,-1.8*ech,0.2*ech)
draw_buoy(x)
pause(3)

