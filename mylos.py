def command_mylos(n_state, v_state, target):
    global previous_theta
    x,y,theta = n_state
    xd, yd = target[:2]

    theta_d = float(theta) - float(np.arctan2(yd - y, xd - x))

    if(theta_d >= np.pi): theta_d -= 2*np.pi
    elif(theta_d <= -np.pi): theta_d += 2*np.pi

    k = -2

    variacao = theta_d - previous_theta

    ac = np.array([[0],[0],[k*theta_d + k*(variacao/0.1)]])
    previous_theta = theta_d
    torque = M @ ac + C(v_state) @ v_state + D(v_state) @ v_state

    return torque

def get_target(target, n_state, traj):
    if(traj == 'l'):
        return get_chemin_target_line(target, n_state)
    else: 
        global idx
        idx, res = get_chemin_target_circle(idx, target, n_state)
        return res