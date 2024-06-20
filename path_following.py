import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from scipy.integrate import odeint

from roblib import *
from utils import *
from path_generator import *


def init_path():
    # Define your waypoints (x, y)
    waypoints = 5*np.array([
        [0, 0],
        [1, 3],
        [2, 1],
        [3, 4],
        [4, 2]
    ])

    # Separate x and y coordinates
    x = waypoints[:, 0]
    y = waypoints[:, 1]

    # Create a cubic spline interpolation
    cs = CubicSpline(x, y)
    
    n_points = 150

    # Generate more points for a smoother curve
    x_smooth = np.linspace(min(x), max(x), n_points).reshape((n_points,1))
    y_smooth = cs(x_smooth).reshape((n_points,1))
    
    
    # path = np.zeros((n_points,5))
    path = np.hstack((x_smooth,y_smooth))
    path[:,0] = x_smooth[:,0]
    path[:,1] = y_smooth[:,0]
    
    
    s = np.zeros((n_points,1))
    theta_path = np.zeros((n_points,1))
    c_path = np.zeros((n_points,1))
    
    for i in range(1,len(path)):
        Xi = path[i,:2]
        Xi_1 = path[i-1,:2]
        ds = np.linalg.norm(Xi-Xi_1)
        s[i:] += ds
        
        theta_path[i,0] = np.arctan2(Xi[1]-Xi_1[1],Xi[0]-Xi_1[0])
        
        diff_ang = sawtooth(theta_path[i,0]-theta_path[i-1,0])
        c_path[i,0] = diff_ang/ds
    
    path = np.hstack((path,s,theta_path,c_path))
    

    return path



def PF_Kin_Computation(X, V_B, s1, y1, ThetaRef, c, u_des,respect_saturation = True,condi = False):
    k_delta = 1
    R = 0.1
    L = 0.25
    wmax = 10*np.pi
    umax = (R*wmax-L/4)/(1+L/4)
    u_des = 1
    rs = 2
    
    Theta = sawtooth(X[2] - ThetaRef)
    
    k1 = 1
    k2 = 3.5
    
    dot_s = V_B[0,0] * np.cos(Theta) + k1 * s1
    Delta = -(np.pi / 2) * np.tanh(k_delta * y1)
    
    ##### Disrespect saturation
    if not respect_saturation:
        # print(y1)
        # y1*=-1
        # Delta *= -1
        # c*=-1
        dot_y1 = -c*dot_s*s1 + V_B[0,0]*np.sin(Theta)
        dot_delta = -(np.pi / 2) * (1 - (np.tanh(k_delta * y1)) ** 2) * k_delta * dot_y1
        Err_angulaire = sawtooth(Theta - Delta)

        u_2 = 1*dot_delta - k2 * Err_angulaire + 1 * c * dot_s
        # u_2 = 10*sawtooth(ThetaRef-X[2])
        u_1 = u_des
    
    # if condi:
    #     # print('diff :',diff_angle)
    #     print(Delta)

    ##### Respect saturation
    else:
        Delta_p = -(np.pi/2)*(1-(np.tanh(k_delta*y1))**2)*k_delta
        # Delta_p = -np.pi/2 * k_delta * (1/np.cosh(k_delta*y1))**2
        K1 = k1/(1+(s1*c*(1-s1*Delta_p))**2)
        

        diff_angle = sawtooth(Theta-Delta)
        fr = -k2*np.tanh(diff_angle) + K1*s1*c*(1-s1*Delta_p)
        gr = c*np.cos(Theta)*(1-s1*Delta_p) + Delta_p*np.sin(Theta)
        v = (u_des/umax)*(R*wmax - L*(0.25+fr**2))/(1+L*(0.25+gr**2))
        r = fr + gr*v
        # if condi:
        #     r = 10*sawtooth(ThetaRef-X[2])
        #     print(r)
        u_1,u_2 = v,r
    
    
    return u_1, u_2, dot_s, Delta

def SplinePath(s, Xrob, Yrob, XYDpath):
    i = 0
    XYDpath = XYDpath.T
    while XYDpath[2, i] <= s:
        i += 1
        if i >= len(XYDpath.T)-1:
            # print('End of path reached')
            break
    
    Xrab = XYDpath[0, i-1] + (XYDpath[0, i] - XYDpath[0, i-1]) * (s - XYDpath[2, i-1]) / (XYDpath[2, i] - XYDpath[2, i-1])
    Yrab = XYDpath[1, i-1] + (XYDpath[1, i] - XYDpath[1, i-1]) * (s - XYDpath[2, i-1]) / (XYDpath[2, i] - XYDpath[2, i-1])
    
    # Xrab = XYDpath[0, i-1]
    # Yrab = XYDpath[1, i-1]
    
    Diff_angulaire = sawtooth(XYDpath[3, i] - XYDpath[3, i-1])

    ThetaRef = sawtooth(XYDpath[3, i-1] + (Diff_angulaire) * (s - XYDpath[2, i-1]) / (XYDpath[2, i] - XYDpath[2, i-1]))

    c = XYDpath[4, i-1] + (XYDpath[4, i] - XYDpath[4, i-1]) * (s - XYDpath[2, i-1]) / (XYDpath[2, i] - XYDpath[2, i-1])

    s1 = -np.cos(ThetaRef) * (Xrab - Xrob) - np.sin(ThetaRef) * (Yrab - Yrob)
    y1 = np.sin(ThetaRef) * (Xrab - Xrob) - np.cos(ThetaRef) * (Yrab - Yrob)
        
    return Xrab, Yrab, ThetaRef, c, s1, y1

def sig_OA(P_R,C_OA,Psi_PF,rs):
    dist2obs = np.linalg.norm(P_R-C_OA)
    # if sigmaOA:
    #     return True
    
    # if sigmaOA:
    #     delta = sawtooth(Psi_PF-cap_robot)
    #     cond_direction = np.cos(delta)>=0
    # else: cond_direction = True
    
    #Check ψPF ∈ SOA:
    psiPF_in_circle = False
    tmp = P_R + 0.01*np.array([[np.cos(Psi_PF),np.sin(Psi_PF)]]).T
    d = dist2obs
    if dist2obs<2*rs:
        if dist2obs<=rs:
            psiPF_in_circle=True
        else:
            while np.linalg.norm(tmp-C_OA)<d:
                d = np.linalg.norm(tmp-C_OA)
                if d <= rs:
                    psiPF_in_circle = True
                tmp += 0.1*np.array([[np.cos(Psi_PF),np.sin(Psi_PF)]]).T
    
    # if not (psiPF_in_circle and cond_direction):
    #     return False
    return dist2obs<2*rs and psiPF_in_circle

def S_OA(alphak,sig_OA):
    if 0<=alphak<=np.pi and sig_OA:
        return 2
    if -np.pi<=alphak<=0 and sig_OA:
        return 1
    else :
        return 0

def psi_OA(PR,COA,SOA):
    vec = COA-PR
    psi = np.arctan2(vec[1,0],vec[0,0])
    return sawtooth(psi + (-1)**(SOA-1)*np.pi/2)

def f(X,t,u):
    # print(u.shape)
    x,y,th = X
    u1,u2 = u.flatten()
    return np.array([[u1*cos(th),u1*sin(th),u2]]).flatten()

if __name__=='__main__':
    # path = init_path()
    path,path_img = process_path()

    xmin = np.min(path[:,0]) - 10
    ymin = np.min(path[:,1]) - 10
    xmax = np.max(path[:,0]) + 10
    ymax = np.max(path[:,1]) + 10

    # X = np.array([200,20,0])
    X = np.array([path[0,0],path[0,1],0])
    traj = X[:2].reshape((2,1))
    dX = np.zeros_like(X)
    V_B = np.zeros((3,1))

    k_delta = 1
    R = 0.1
    L = 0.25
    wmax = 100*np.pi
    umax = (R*wmax-L/4)/(1+L/4)
    u_des = 10
    rs = 2


    obstacles = generate_obstacles(6,path,L=20,random=True)
    COA = None
    SOA = None
    sigmaOA=False
    psiOA = None

    s = 0
    tmax = 200
    dt = 0.01

    # ax = init_figure(xmin,xmax,ymin,ymax)
    ax = init_figure(0,path_img.shape[0],0,path_img.shape[1])
    live = False

    hist_delta = []
    hist_y = []
    for t in np.arange(0,tmax,dt):
        xr,yr,thetar = X.flatten()
        ########### Node spline path ###########
        Xrab, Yrab, ThetaRef, c, s1, y1 = SplinePath(s,xr,yr,path)    
        ############################################
        
        ########### node pf ###########
        u1, u2, dot_s, Delta = PF_Kin_Computation(X,V_B,s1,y1,ThetaRef,c,u_des)
        s += dt*dot_s
        ############################################
        
        PsiPF = sawtooth(ThetaRef + Delta)
        
        ########### Node lidar ###########
        PR = X[:2].reshape((2,1))
        # dk, pt = distance_point_to_obstacle(X,obstacles)
        obs_in_range,pt,ak,dk,lidar = sim_lidar(X,obstacles,num_rays=16)
        if COA is None or not obs_in_range:
            COA = pt
        elif dk<np.linalg.norm(COA-PR):
            COA = pt
        PR_COA = COA-PR
        sigmaOA = sig_OA(PR,COA,PsiPF,rs)
        ############################################
        
        
        if sigmaOA:
            if SOA is None:
                # alpha = sawtooth(np.arctan2(PR_COA[1,0],PR_COA[0,0])-X[2])
                SOA = S_OA(ak,sigmaOA)

            psiOA = psi_OA(PR,COA,SOA)

            y1_obs = ((-1)**(SOA-1))*(np.linalg.norm(PR_COA)-rs)
            # y1_obs = -(np.linalg.norm(PR_COA)-rs)
            # print(y1_obs)
            u1, u2, _, Delta = PF_Kin_Computation(X,V_B,0,y1_obs,psiOA,0*1/rs,u_des,respect_saturation=True,condi=True)
            hist_y.append(y1_obs)

        else:
            SOA = None
            hist_y.append(y1)
        
        hist_delta.append(Delta)
        
        X = odeint(f, X, np.arange(0,3*dt/2,dt), args=(np.array([[u1,u2]]), ))[-1]
        X[2] = sawtooth(X[2])
        
        
        dX = f(X,t,np.array([[u1,u2]]))
        Rot = np.array([[cos(X[2]), -sin(X[2]), 0],
                        [sin(X[2]),  cos(X[2]), 0],
                        [0,          0,         1]])
        V_B = np.linalg.inv(Rot)@(dX.reshape((3,1)))
        
        traj = np.hstack((traj,X[:2].reshape((2,1))))
        
        if live:      
            clear(ax)
            ax.set_xlim(X[0] - 10 / 2, X[0] + 10 / 2)
            ax.set_ylim(X[1] - 10 / 2, X[1] + 10 / 2)
            
            # draw_obstacles(ax,obstacles)
            draw_lidar(ax,X,lidar,obstacles)
            ax.plot(path[:,0], path[:,1])
            ax.scatter(Xrab,Yrab, label='rabbit')
            
            ax.plot([X[0],X[0]+20*np.cos(PsiPF)],[X[1],X[1]+20*np.sin(PsiPF)],label=r'$\psi_{PF}$')
            if sigmaOA:
                ax.plot([X[0],X[0]+20*np.cos(psiOA)],[X[1],X[1]+20*np.sin(psiOA)],label=r'$\psi_{OA}$')
                
            plot_circle(ax,X[0],X[1],2*rs)
            
            plot_circle(ax,COA[0,0],COA[1,0],rs)
            ax.scatter(COA[0,0],COA[1,0],label='Point le plus proche')

            draw_tank(X,r=0.2,w=0.4)
            ax.legend()
            pause(dt)
        
    # clear(ax)
    draw_obstacles(ax,obstacles)
    ax.imshow(path_img)
    plt.gca().invert_yaxis()
    ax.plot(path[:,0], path[:,1],label='Path')
    ax.plot(traj[0,:], traj[1,:],label='Robot trajectory')
    ax.set_title("Path following with Safe Maneuvering Zone")
    ax.legend()

    # plt.figure()
    # plt.plot(traj[0,:], traj[1,:],label='Robot trajectory')
    # plt.show()
    # fig2 = figure()
    # ax_delta = fig2.add_subplot(111)
    # ax_delta.xmin=0
    # ax_delta.xmax=50/dt
    # ax_delta.ymin=-np.pi
    # ax_delta.ymax=np.pi

    # fig3 = figure()
    # ax_y = fig3.add_subplot(111)
    # ax_y.xmin=0
    # ax_y.xmax=50/dt
    # ax_y.ymin=-5
    # ax_y.ymax=5

    # ax_delta.plot(hist_delta)
    # ax_delta.set_title("Delta")
    # ax_y.plot(hist_y)
    # ax_y.set_title("y")

    plt.show()
