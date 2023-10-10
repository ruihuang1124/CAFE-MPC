def approx_equal(a, b):
    if abs(a-b) < 1e-6:
        return True
    return False

def approx__leq(a, b):
    if (a < b) or approx_equal(a,b):
        return True
    return False

def approx_geq(a, b):
    if (a > b) or approx_equal(a,b):
        return True
    return False

from cgi import print_arguments
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import lcm

import sys
sys.path.append('../../')
from lcmtypes.python.wbTraj_lcmt import wbTraj_lcmt

def plot_gait_schedule(tau):
    " Visualize the gait schedule given a state trajectory as input"
    plt.rcdefaults()
    fig, ax = plt.subplots()

    time = [s.time for s in tau]
    foot_idx = [3, 2, 1, 0]
    foot_names = ["HL","HR","FL","FR"]
    colormap = ["y", "k"]
    for foot in range(4):
        ax.barh(foot_idx[foot], time[1] - time[0], left=[0] + time[0:-1], 
                        color=[colormap[s.contact[foot]] for s in tau])    
    ax.set_yticks([0,1,2,3],foot_names)
    # ax.set_xlim([time[0], time[-1]])
    ax.set_xlabel("time (s)")
    plt.show()

def plot_com_pos(time, pos_tau):
    _, axs = plt.subplots(1, 3)    
    axs[0].plot(time, [p[0] for p in pos_tau])
    axs[0].set_xlabel('time (s)')
    axs[0].set_ylabel('x (m)')

    axs[1].plot(time, [p[1] for p in pos_tau])
    axs[1].set_xlabel('time (s)')
    axs[1].set_ylabel('y (m)')

    axs[2].plot(time, [p[2] for p in pos_tau])
    axs[2].set_xlabel('time (s)')
    axs[2].set_ylabel('z (m)')
    plt.show()

def plot_eul(time, eul_tau):
    _, axs = plt.subplots(1, 3)    
    axs[0].plot(time, [eul[0] for eul in eul_tau])
    axs[0].set_xlabel('time (s)')
    axs[0].set_ylabel('yaw (rad)')

    axs[1].plot(time, [eul[1] for eul in eul_tau])
    axs[1].set_xlabel('time (s)')
    axs[1].set_ylabel('pitch (rad)')

    axs[2].plot(time, [eul[2] for eul in eul_tau])
    axs[2].set_xlabel('time (s)')
    axs[2].set_ylabel('roll (rad)')
    plt.show()

def plot_com_vel(time, vel_tau):
    _, axs = plt.subplots(1, 3)
    axs[0].plot(time, [vel[0] for vel in vel_tau])
    axs[0].set_xlabel('time (s)')
    axs[0].set_ylabel('vx (m/s)')

    axs[1].plot(time, [vel[1] for vel in vel_tau])
    axs[1].set_xlabel('time (s)')
    axs[1].set_ylabel('vy (m/s)')

    axs[2].plot(time, [vel[2] for vel in vel_tau])
    axs[2].set_xlabel('time (s)')
    axs[2].set_ylabel('vz (m/s)')
    plt.show()

def plot_jnt_position(time, jnt_pos_tau, leg):
    _, ax = plt.subplots()
    ax.plot(time, [jnt_pos[3*leg] for jnt_pos in jnt_pos_tau])
    ax.plot(time, [jnt_pos[3*leg+1] for jnt_pos in jnt_pos_tau])
    ax.plot(time, [jnt_pos[3*leg+2] for jnt_pos in jnt_pos_tau])
    ax.set_xlabel("time")
    ax.set_ylabel("Joint anlges (radian)")
    ax.legend(["abad", "hip", "knee"])
    plt.show()


LEG_INDEX = {'FR': 0, 'FL': 1, 'HR': 2, 'HL': 3}
def plot_swing_height(time, z_tau):
    """ Plot foot positions w.r.t. time
    """
    _, axs = plt.subplots(2, 2)

    for legname, legid in LEG_INDEX.items():
        axs[int(legid/2), legid % 2].plot(time, [z[legid] for z in z_tau])
        axs[int(legid/2), legid % 2].set_xlabel('time (s)')
        axs[int(legid/2), legid % 2].set_ylabel('z (m)')
        axs[int(legid/2), legid % 2].set_title(legname)
    plt.show()

COLOR = ['r', 'b', 'm', 'g']
def plot_foot_positions(time, pf_tau):
    _, axs = plt.subplots(2, 2)

    for legname, legid in LEG_INDEX.items():
        pf_x = [pf[3*legid] for pf in pf_tau]
        pf_y = [pf[3*legid+1] for pf in pf_tau]
        pf_z = [pf[3*legid+2] for pf in pf_tau]
        
        axs[int(legid/2), legid % 2].plot(pf_x, pf_z)
    plt.show()

def plot_footPosition_and_CoM(pf_tau, p_tau):
    _, axs = plt.subplots(2, 2)

    for legname, legid in LEG_INDEX.items():
        pf_x = [pf[3*legid] for pf in pf_tau]
        pf_y = [pf[3*legid+1] for pf in pf_tau]
        pf_z = [pf[3*legid+2] for pf in pf_tau]
        x_com = [p[0] for p in p_tau]
        y_com = [p[1] for p in p_tau]
        z_com = [p[2] for p in p_tau]

        axs[int(legid/2), legid % 2].plot(pf_x, pf_z)
        axs[int(legid/2), legid % 2].plot(x_com, z_com)
    plt.show()

def animate_footPositions_and_CoM(leg, pf_tau, p_tau):
    fig, ax = plt.subplots()
    legName = next((name for name, id in LEG_INDEX.items() if id == leg), None)
    
    pf_x_array = [pf[3*leg] for pf in pf_tau]        
    pf_z_array = [pf[3*leg+2] for pf in pf_tau]
    x_com_array = [p[0]+0.2 for p in p_tau]
    z_com_array = [0.01 for p in p_tau]
    
    ax.plot(pf_x_array, pf_z_array, "r-")
    ax.plot(x_com_array, z_com_array, "b-")
    point_pf = ax.plot(pf_x_array[0], pf_z_array[0], "ro")[0]
    point_com = ax.plot(x_com_array[0], z_com_array[0], "bo")[0]

    def update(frame):
        pf_x = pf_x_array[frame]
        pf_z = pf_z_array[frame]
        x_com = x_com_array[frame]
        z_com = z_com_array[frame]

        point_pf.set_xdata(pf_x)
        point_pf.set_ydata(pf_z)
        point_com.set_xdata(x_com)
        point_com.set_ydata(z_com)
        return(point_pf, point_com)
    
    ani = animation.FuncAnimation(fig=fig, func=update, frames=len(pf_x_array), interval=100)
    plt.show()

def flip_hip_knee_direction(jnt_pos):
    jnt_pos[np.array([2,5,8,11])] *= -1.0
    jnt_pos[np.array([1,4,7,10])] *= -1.0 
    return jnt_pos

def flip_left_and_right(vec12):
    vec12_flipped = vec12.copy()
    vec12_flipped[np.array([0,1,2])] = vec12[np.array([3,4,5])]
    vec12_flipped[np.array([3,4,5])] = vec12[np.array([0,1,2])]
    vec12_flipped[np.array([6,7,8])] = vec12[np.array([9,10,11])]
    vec12_flipped[np.array([9,10,11])] = vec12[np.array([6,7,8])]
    return vec12_flipped

def flip_contact_left_right(contact):
    c_flipped = contact.copy()
    c_flipped[0] = contact[1]
    c_flipped[1] = contact[0]
    c_flipped[2] = contact[3]
    c_flipped[3] = contact[2]
    return c_flipped

def rearrange_traj_for_viz(wbtraj_lcmt):
    for k in range(wbtraj_lcmt.sz):                
        qJ = flip_hip_knee_direction(np.array(wbtraj_lcmt.qJ[k]))                
        wbtraj_lcmt.qJ[k] = list(flip_left_and_right(qJ))        

def write_traj_to_file(time, pos, eul, vel, eulrate, pf, vf, jnt, contact):
    base = np.hstack((eul, pos, eulrate, vel))   
    
    # Rearrange leg order and reverse rotation for urdf file used in MPC
    for k in range(len(jnt)):
        jnt_k = jnt[k]
        jnt_k = flip_hip_knee_direction(jnt_k)
        jnt_k = flip_left_and_right(jnt_k)

        pf_k = flip_left_and_right(pf[k])
        vf_k = flip_left_and_right(vf[k])
        c_k = flip_contact_left_right(contact[k])

        jnt[k] = jnt_k
        pf[k] = pf_k
        vf[k] = vf_k
        contact[k] = c_k
    
    np.savetxt("data/time.csv", np.asarray(time), delimiter=",", fmt='%8.4f')
    np.savetxt("data/body_state.csv", base, delimiter=",", fmt='%8.4f')
    np.savetxt("data/ee_pos.csv", np.asarray(pf), delimiter=",", fmt='%8.4f')
    np.savetxt("data/jnt.csv", np.asarray(jnt), delimiter=",", fmt='%8.4f')
    np.savetxt("data/contact.csv", np.asarray(contact), delimiter=",", fmt='%u')
    np.savetxt("data/ee_vel.csv", np.asarray(vf), delimiter=",", fmt='%8.4f')

def publish_trajectory_lcm(pos_tau, eul_tau, vel_tau, eulrate_tau, jnt_tau, contact_tau):
    lcm_ = lcm.LCM()
    wbtraj_lcmt = wbTraj_lcmt()
    traj_sz = len(pos_tau)
    wbtraj_lcmt.sz = len(pos_tau)
    for k in range(traj_sz):
        wbtraj_lcmt.pos.append(list(pos_tau[k]))
        wbtraj_lcmt.eul.append(list(eul_tau[k]))
        wbtraj_lcmt.vWorld.append(list(vel_tau[k]))
        wbtraj_lcmt.eulrate.append(list(eulrate_tau[k]))
        wbtraj_lcmt.qJ.append(list(jnt_tau[k]))
        wbtraj_lcmt.qJd.append([0]*12)
        wbtraj_lcmt.torque.append([0]*12)
        wbtraj_lcmt.contact.append(list(contact_tau[k]))
    
    rearrange_traj_for_viz(wbtraj_lcmt)

    lcm_.publish("visualize_wb_traj", wbtraj_lcmt.encode())