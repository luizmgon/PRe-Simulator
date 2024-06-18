from roblib import *  # available at https://www.ensta-bretagne.fr/jaulin/roblib.py

ax=init_figure(-5,15,-5,15)
m=20
p=10*rand(2,m)
plot(p[0,:],p[1,:],'ob')
pause(1)

