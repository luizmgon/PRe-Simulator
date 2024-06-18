import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline



M = np.array([[80.2815,     0,       0],
             [       0, 157.5,      0],
             [       0,    0, 41.9327]], dtype=float)

M_model = np.array([[80.2815,     0,       0],
             [       0, 157.5,      11],
             [       0,    11, 41.9327]], dtype=float)

tracking_point_distance = 1

H = np.array([[1,0,0],
              [0,0,0],
              [0,1/tracking_point_distance,1]], dtype=float)
T = np.array([[1,0,0],
              [0,1,tracking_point_distance],
              [0,0,1]], dtype=float)

def dJ(n, v):
    phi = float(n[2])
    phi_dot = float(v[2])
  
    return phi_dot * np.array([
        [-np.sin(phi), -np.cos(phi), 0],
        [np.cos(phi), -np.sin(phi), 0],
        [0, 0, 0],
    ], dtype=float)

def J(n):
    phi = float(n[2])
    ret = np.array([[np.cos(phi), -np.sin(phi), 0],
                    [np.sin(phi),  np.cos(phi), 0],
                    [0, 0, 1]], dtype=float)
    
    return ret

def C(v_state):
    u,v,r = v_state.flatten()
    Xudot = -5.2815
    Yvdot = -82.5

    c = np.array([[0         ,         0,  Yvdot * v],
                  [0         ,         0, -Xudot * u],
                  [-Yvdot * v, Xudot * u,          0]], dtype=float) 
    return c
      
def D(v):
    Xu = -77.5544
    Yv = -157.5000
    Nr = -41.9327
    return -np.array([[Xu, 0, 0],
                      [0, Yv, 0],
                      [0, 0, Nr]], dtype=float)

def command_h(n_state, v_state, target, d_target, state_error_integral, dt, speed_error_integral, previous_speed_error):

    invT = np.linalg.inv(T)
    invJ = np.linalg.inv(J(n_state))
    
    state_error = target - n_state
    state_error_integral += state_error * dt

    if(np.all(np.abs(state_error[:2])< 0.5)): 
        # print("quietei")
        torque = np.zeros((3,1))
        return torque, state_error_integral, speed_error_integral, previous_speed_error


    Kp_state = np.array([[1], [1], [1]])
    Ki_state = np.array([[0], [0], [0]])

    n_correction = d_target + Kp_state * state_error + Ki_state * state_error_integral

    error_repere_body = invT @ invJ @ state_error

    vref = H @ invT @ invJ @ n_correction

    speed_error =  vref - v_state
    speed_error_integral += speed_error * dt

    speed_error_derivative = 0
    if(previous_speed_error is not None):
        speed_error_derivative = (speed_error - previous_speed_error) / dt

    Kp_speed = np.array([[1], [1], [1]])
    Ki_speed = np.array([[0], [0], [0]])
    Kd_speed = np.array([[0.5], [1], [0.8]])

    ac = Kp_speed * speed_error + Ki_speed * speed_error_integral + Kd_speed * speed_error_derivative

    Cv = C(v_state)
    Dv = D(v_state)

    torque = M @ ac + Cv @ v_state + Dv @ v_state
    # torque[0] = 0
    torque[1] = 0

    return torque, state_error_integral, speed_error_integral, speed_error

def command_sous_diag(n_state, v_state, k):

    x,y,phi = n_state.flatten()
    u,v,r = v_state.flatten()


    Cv = C(v_state)
    Dv = D(v_state)

    dv = (-Cv[1,2]*r - Dv[1,1]*v)/M[1,1]

    A = np.array([[np.cos(phi), -u*np.sin(phi)-v*np.cos(phi)],
                  [np.sin(phi), u*np.cos(phi)-v*np.sin(phi)]])
    
    B = [[-dv*np.sin(phi)], [dv*np.cos(phi)]]
   
    inv = np.linalg.inv(A)
    dif = (k - B)
    ac = inv @ dif

    ac[1] = -10*(r - ac[1])

    ac = np.array([[float(ac[0])],[dv],[float(ac[1])]])
    torque = M @ ac + Cv @ v_state + Dv @ v_state
    torque = M @ ac + Cv @ v_state + Dv @ v_state

    return torque

def command_sous(n_state, v_state, k):

    x,y,phi = n_state.flatten()
    u,v,r = v_state.flatten()


    Cv = C(v_state)
    Dv = D(v_state)

    alpha = -M_model[1,2] / M_model[1,1]

    # dv = (-Cv[1,2]*r - Dv[1,1]*v)/M[1,1]

    A = np.array([[np.cos(phi), -alpha * np.sin(phi)],
                  [np.sin(phi), alpha * np.cos(phi)]])
    
    B = [[-u*r*np.sin(phi) - v*r*np.cos(phi) + (np.sin(phi)/M[1,1]) * (Cv[1,2]*r + Dv[1,1] * v)], 
         [u*r*np.cos(phi) - v*r*np.sin(phi) + (np.cos(phi)/M[1,1]) * (-Cv[1,2]*r - Dv[1,1]*v)]]
   
    k[0] = k[0] / 10
    k[1] = k[1] / 10

    inv = np.linalg.inv(A)
    p = np.max(inv)
    dif = (k - B)
    ac = inv @ dif

    dr = float(ac[1])
    dv = (-M_model[1,2] * dr - Cv[1,2] * r -Dv[1,1] * v) / M_model[1,1]

    ac = np.array([[float(ac[0])],[dv],[float(ac[1])]])

    torque = M_model @ ac + Cv @ v_state + Dv @ v_state

    print(torque)



    return torque

def command_los(n_state, v_state, target, previous_theta_d):
    
    x,y,theta = n_state.flatten()
    xd, yd = target.flatten()[:2]


    desired = np.arctan2(yd - y, xd - x)
    theta_d = theta - desired

    if(theta_d >= np.pi): theta_d -= 2*np.pi
    elif(theta_d <= -np.pi): theta_d += 2*np.pi

    k = -12
    # print(theta_d)

    variacao = 0
    if(previous_theta_d != None):
        variacao = theta_d - previous_theta_d
    
    previous_theta_d = theta_d

    ac = np.array([[0],[0],[k*theta_d + 5*k*(variacao)]])

    # ac = np.array([[0],[0],[k*theta_d]])
    torque = M @ ac + C(v_state) @ v_state + D(v_state) @ v_state

    return previous_theta_d, torque

def command_full(n_state, v_state, target, d_target, dd_target):
    dn_state = J(n_state) @ v_state
    c = dd_target + 2 * (d_target - dn_state) + (target - n_state) #Fblin
    # c =  dd_target +(d_target - dn_state) +5*np.sign((d_target - dn_state) + (target - n_state)) #slid1
    # c =  10*np.sign((d_target - dn_state) + (target - n_state)) #slid2


    ac = np.linalg.inv(J(n_state)) @ (-dJ(n_state, v_state) @ v_state + c)
    torque = M @ ac + C(v_state) @ v_state + D(v_state) @ v_state
    # print(torque)

    return torque

def command_auv(v_state, error_state, u_target, path_points, s_values, mu_values, phi_values):
    # State
    u,v,r = v_state.flatten()
    s1, y1, phi, s = error_state.flatten()

    # if(s1 > 5): s1 = 5
    # if(s1 <-5): s1 = -5
    # if(y1 > 5): y1 = 5
    # if(y1 <-5): y1 = -5

    # Gains
    k1 = 1
    k2 = 1
    k3 = 1
    k4 = 1
    k5 = 1

    # Masses
    m = 75
    m_u = M_model[0,0]
    m_ur = M_model[1,2] #FIXME
    m_v = M_model[1,1]

    Yv = -157.5
    Xuu = Yv/10
    Xvv = -5.2815
    Nv = -30
    Nvv = 11
    Nr = -30
    Y_vv = -300

    d_v = -Yv * u * v - Y_vv * v * np.abs(v)

    # Speed derivatives
    du = -k4*(u - u_target)
    dd_u = -k4 * du
    dv = (-m_ur*u*r - d_v ) / m_v
    v_t = np.sqrt(u**2 + v**2)
    d_vt = (u*du + v*dv) / v_t

    # Damping terms
    d_u = -Xuu * u**2 - Xvv * v**2
    d_r = -Nv*u*v - Nvv * v * np.abs(v) - Nr * u *r
    d_dv = -Yv * (du * v + u * dv) - 2*Y_vv * v * dv


    ds = v_t * np.cos(phi) + k2 * s1

    # ds = np.clip(ds, -5, 5)

    if(s == 0 and ds < 0): ds = 0

    # xs, ys, phif, curv, dcurv_dt, g_c = path_desc(s, ds,s_values, mu_values, phi_values)
    xs, ys, _, phif, curv, g_c = path_interrogation(s, path_points)

    dcurv_dt = g_c * ds

    #Beta
    beta = np.arctan2(v, u)
    d_beta = (dv*u - du*v) / (u**2 + v**2)

    d_phi = r + d_beta - curv * ds

    # Error state derivatives
    d_s1 = -ds*(1-curv*y1) + v_t*np.cos(phi)
    d_y1 = -curv * ds * s1 + v_t * np.sin(phi)
    dd_s = d_vt * np.cos(phi) - v_t * np.sin(phi) * d_phi + k2 * d_s1
    # dd_s = np.clip(dd_s, -5, 5)
    dd_y1 = -dcurv_dt * ds * s1 - curv * dd_s * s1 - curv * ds * d_s1 + d_vt * np.sin(phi) + v_t*np.cos(phi) * d_phi

    # Delta derivatives
    k_delta = 0.5
    phi_a = np.pi/2
    delta = -phi_a * np.tanh(k_delta * y1)
    diff = phi - delta
    print(diff)
    d_delta = (-phi_a * k_delta * d_y1)/  (np.cosh(k_delta * y1))**2
    dd_delta = -phi_a * k_delta*(dd_y1/(np.cosh(k_delta * y1)**2) - 2 * k_delta * d_y1**2 * np.tanh(k_delta * y1) / (np.cosh(k_delta*y1))**2)

    rd = d_delta - d_beta - k1 * sawtooth(phi - delta) + curv * ds

    f_alpha = dd_delta - k1*(d_phi - d_delta) + curv*dd_s + g_c * ds 
    + dd_u*v/v_t**2 + 2 * d_vt * d_beta / v_t + (u/v_t**2)*((m_ur/m_v)*du * r + d_dv/m_v)
    
    dr = (f_alpha - k3*(r - rd) - k5*sawtooth(phi - delta))/(1 - (m_ur/m_v)*(np.cos(beta)**2))

    Fu = m_u * du + d_u
    N_r = m * dr + d_r

    # if(Fu >= 0):
    #     Fu = max(Fu, 1000)
    # else: Fu = min(Fu, -1000)
    # if(N_r >= 0):
    #     N_r = max(N_r, 1000)
    # else: N_r = min(N_r, -1000)

    Fu = np.clip(Fu, -1000, 1000)
    N_r = np.clip(N_r, -1000, 1000)


    torque = np.array([[Fu], [0], [N_r]])

    return torque, ds

def draw_p():
    def dmu_ds(s, mu, a2, a4, b1):
        return 1 / np.sqrt((2*a2*mu + 4*a4*mu**3)**2 + b1**2)

    def euler_method(func, s0, mu0, a2, a4, b1, h, mu_max):
        s_values = [s0]
        mu_values = [mu0]

        dxs_dmu = 2*a2*mu0 + 4*a4*mu0**3
        dys_dmu = b1

        phi0 = np.arctan2(dys_dmu, dxs_dmu)
        phi_values = [phi0]
        
        s = s0
        mu = mu0
        
        while(mu < mu_max):
            mu += func(s, mu, a2, a4, b1) * h
            s += h
            dxs_dmu = 2*a2*mu + 4*a4*mu**3
            dys_dmu = b1

            phi = np.arctan2(dys_dmu, dxs_dmu)
            s_values.append(s)
            mu_values.append(mu)
            phi_values.append(phi)
        
        return s_values, mu_values, phi_values

    def calculate_xy(mu_values, a2, a4, b1):
        x_values = [a2 * mu**2 + a4 * mu**4 for mu in mu_values]
        y_values = [b1 * mu for mu in mu_values]
        return x_values, y_values

    # Parâmetros
    a2 = -0.02
    a4 = 3 * 10**-6
    b1 = 0.2
    s0 = 0
    mu0 = -90
    h = 0.01  # Tamanho do passo
    mu_max = 90

    # Execução do método de Euler
    s_values, mu_values, phi_values = euler_method(dmu_ds, s0, mu0, a2, a4, b1, h, mu_max)
    # # Cálculo de x e y em função de mu
    x_values, y_values = calculate_xy(mu_values, a2, a4, b1)

    return x_values, y_values, s_values, mu_values, phi_values


def get_mu(s_in,s_values, mu_values):

    indice = 0
    s_atual = s_values[indice]
    while(s_atual < s_in):
        indice += 1
        s_atual = s_values[indice]

    return mu_values[indice]

def get_path_points():
    # x = np.array([0, 0,  0,  0,  0,  0,  0])
    # x = np.array([0, -10, -20, -10, 0, -10, -20])
    x = np.array([2, -3, -8, -3, 2, -3, -8])
    y = np.array([  0,   5, 10, 15, 20,  25, 30])

    cs = CubicSpline(y, x)

    y = np.linspace(0, y[-1], 10000)
    x = cs(y)

    distances = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
    s = np.zeros_like(x)
    s[1:] = np.cumsum(distances)

    # Calcular as derivadas usando diferenças finitas
    dx = np.gradient(x)
    dy = np.gradient(y)
    ddx = np.gradient(dx)
    ddy = np.gradient(dy)

    phi_f = np.arctan2(dy, dx)

    # Calcular a curvatura
    curvature = dx * ddy - dy * ddx / (dx**2 + dy**2)**1.5

    # Calcular as diferenças em s e c
    ds = np.diff(s)
    dc = np.diff(curvature)

    # Calcular o gradiente de c em relação a s
    g_c = dc / ds
    g_c = np.insert(g_c, 0, 0)

    return np.vstack((x, y, s, phi_f, curvature, g_c, dx, dy))

    if(True):
            # Plotar o caminho, curvatura, phi_f e g_curvature
        plt.figure(figsize=(12, 10))

        # Plotar o caminho
        plt.subplot(4, 1, 1)
        plt.plot(x, y, 'b-', label='Caminho')
        plt.scatter(x, y, color='red', label='Pontos de Controle')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Caminho, Curvatura, phi_f e Gradiente de Curvatura')
        plt.legend()
        plt.grid(True)

        # Plotar a curvatura
        plt.subplot(4, 1, 2)
        plt.plot(s, curvature, 'g-', label='Curvatura')
        plt.xlabel('Distância ao Longo do Caminho (s)')
        plt.ylabel('Curvatura')
        plt.legend()
        plt.grid(True)

        # Plotar phi_f
        plt.subplot(4, 1, 3)
        plt.plot(s, phi_f, 'r-', label='phi_f')
        plt.xlabel('Distância ao Longo do Caminho (s)')
        plt.ylabel('phi_f')
        plt.legend()
        plt.grid(True)


        # Plotar g_curvature
        plt.subplot(4, 1, 4)
        plt.plot(s, g_c, 'm-', label='Gradiente de Curvatura (g_curvature)')
        plt.xlabel('Distância ao Longo do Caminho (s)')
        plt.ylabel('Gradiente de Curvatura')
        plt.legend()
        plt.grid(True)

        # plt.tight_layout()
        plt.show()

def path_interrogation(s_in, path_points):

    indice = 0
    s_values = path_points[2]
    s_atual = s_values[indice]
    while(s_atual < s_in):
        indice += 1
        s_atual = s_values[indice]

    # print(indice)

    res = path_points[:,indice]

    return res

def sawtooth(phi):
    while phi >  np.pi:
        phi -= 2 * np.pi
    while phi < -np.pi:
        phi += 2 * np.pi
    return phi

def initiate_error_state(n_state, v_state, s0, path_points):
    u,v,_ = v_state.flatten()
    xs, ys, _ , phif = path_interrogation(s0, path_points)[:4]
    
    R = np.array([[np.cos(phif), np.sin(phif), 0],
                    [-np.sin(phif),  np.cos(phif), 0],
                    [0, 0, 1]], dtype=float)
    
    new_state = n_state + np.array([[-xs],[-ys],[-phif + np.arctan2(v , u)]])

    s1,y1,phi = (R @ new_state).flatten()
    phi = sawtooth(phi)
    error_state = np.array([[s1], [y1], [phi], [s0]])
    return error_state

def update_error_state(error_state,v_state, n_state, ds, dt_ctr, path_points):
    u,v,_ = v_state.flatten()

    s = error_state[3,0] + ds * dt_ctr

    if(s < 0): s = 0

    xs, ys, _ , phif = path_interrogation(s, path_points)[:4]
    
    R = np.array([[np.cos(phif), np.sin(phif), 0],
                    [-np.sin(phif),  np.cos(phif), 0],
                    [0, 0, 1]], dtype=float)
    
    new_state = n_state + np.array([[-xs],[-ys],[-phif + np.arctan2(v , u)]])
    s1,y1,phi = (R @ new_state).flatten()
    phi = sawtooth(phi)

    return np.array([[s1], [y1], [phi], [s]])

def path_desc(s, ds,s_values, mu_values, phi_values ):

    mu = get_mu(s,s_values, mu_values)

    #Path description
    a2 = -0.02
    a4 = 3 * 10**-6
    b1 = 0.2

    xs = a2 * mu**2 + a4 * mu**4
    ys = b1 * mu 

    dxs_dmu = 2*a2*mu + 4*a4*mu**3
    ddxs_dmu = 2*a2 + 12 * a4 * mu**2
    dddxs_dmu = 24 * a4 * mu

    dys_dmu = b1

    phif = np.arctan2(dys_dmu, dxs_dmu)
    dphif_dmu = -dys_dmu * ddxs_dmu/((dxs_dmu**2+(dys_dmu)**2))

    ddphif_dmu = (dys_dmu * (2 * dxs_dmu * (ddxs_dmu)**2 - dddxs_dmu * ((dxs_dmu)**2 + dys_dmu) )) / (dxs_dmu**2 + dys_dmu**2)**2

    dmuds = 1/(np.sqrt(dxs_dmu**2 + dys_dmu**2))

    curv = dphif_dmu * dmuds

    ddmuds_dmu = -0.5*(dxs_dmu**2 + dys_dmu**2)**(-3/2) * 2*dxs_dmu*ddxs_dmu

    dcurv_dmu = ddphif_dmu * dmuds + dphif_dmu * ddmuds_dmu
    g_c = dcurv_dmu * dmuds

    dcurv_dt = (ddphif_dmu * dmuds + dphif_dmu * ddmuds_dmu) * dmuds * ds

    return xs, ys, phif, curv, dcurv_dt, g_c


get_path_points()