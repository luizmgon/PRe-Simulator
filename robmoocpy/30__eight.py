from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

def φ0(p1,p2):        
    return -(p1**3+p2**2*p1-p1+p2),-(p2**3+p1**2*p2-p1-p2)
    

dt,s= 0.1,5       
ax=init_figure(-s,s,-s,s)
draw_field(ax,φ0,-s,s,-s,s,0.51)
    

pause(1)    
    
    



