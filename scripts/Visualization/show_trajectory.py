
import time
from turtle import pd
import numpy as np


from scripts.Visualization.mini_cheetah import MiniCheetah
from scripts.Visualization.animator import Animator
from scripts.Reference_python.utils import write_traj_to_file

import lcm
from lcmtypes.python.wbTraj_lcmt import wbTraj_lcmt
import threading

urdf_filename =  "urdf/mini_cheetah_simple_correctedInertia.urdf"
robot = MiniCheetah(urdf_filename)    
animator = Animator(robot)


def visualize_motion_lcm_handler(channel, data):
    print("received visualization lcm message")    
    msg = wbTraj_lcmt.decode(data)

    print("Size of the received trajectory: ", msg.sz)    
    
    t = np.array(msg.time) 
    animator.plot_eul(t, np.array(msg.eul))
    # animator.plot_eulrate(t, np.array(msg.eulrate))
    animator.plot_pos(t, np.array(msg.pos))
    animator.plot_vel(t, np.array(msg.vWorld))
    # animator.plot_joint_torques(t, np.array(msg.torque))
    # # animator.plot_joint_speed(t, np.array(msg.qJd))
    # animator.plot_joint_angles(t, np.array(msg.qJ))
    animator.plot_contact(t, np.array(msg.contact))    
    # animator.plot_defect(t, np.array(msg.defect))
    # animator.plot_roll_abad_rate(t, np.array(msg.qJd), np.array(msg.eulrate))
    # animator.plot_centroid_ang_momentum(t, np.array(msg.hg))
    # animator.plot_centroid_ang_momentum_derivative(t, np.array(msg.dhg))
    animator.show_plots()    
    
    EE_pos = []
    for k in range(msg.sz):
        pf = []
        for leg in range(4):
            qleg = msg.qJ[k][3*leg:3*(leg+1)]
            pf.append(robot.fk(msg.pos[k], msg.eul[k], leg, qleg))
        EE_pos.append(np.hstack(pf))    
    EE_vel = [np.zeros(12) for _ in range(msg.sz)]
    write_traj_to_file(msg.time, msg.pos, msg.eul, 
                       msg.vWorld, msg.eulrate, 
                       EE_pos, EE_vel, 
                       msg.qJ, msg.qJd,
                       msg.contact)

lc = lcm.LCM()
subscription = lc.subscribe("visualize_wb_traj", visualize_motion_lcm_handler)

try:
    while True:        
        lc.handle()
except KeyboardInterrupt:
    pass