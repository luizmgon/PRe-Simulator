from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

s=2
ax=init_figure(-s,s,-s,s)
clear(ax)
draw_disk(ax,array([[0],[0]]),0.2,"blue",1,1)

x=array([[1.22],[0],[0],[1]])
draw_disk(ax,array([[x[0]],[x[1]]]),0.1,"red",1)
pause(1)    





 
