from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

def draw(x,col):
    θ = x[0,0]/r
    a=array([[r*cos(θ)],[r*sin(θ)],[θ+pi/2]])
    draw_tank(a,col)
        
    
L=100
r=L/(2*pi)
ax=init_figure(-20,20,-20,20)
clear(ax)
draw_disk(ax,array([[0],[0]]),r+3,'lightblue')
draw_disk(ax,array([[0],[0]]),r-3,'white')
xa=array([[25],[0]])   
draw(xa,'black')
pause(1)


