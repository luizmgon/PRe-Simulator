from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

    


r=10
a,b,ech = array([[-25,0,pi/2]]).T, array([[25,0,pi/2]]).T, 40      #simu 1
ax=init_figure(-ech,ech,-ech,ech)
clear(ax)
draw_tank(a,"black")
draw_tank(b,"blue")

draw_arc(array([[0],[5]]),array([[4],[6]]),r,'red')

pause(1)
