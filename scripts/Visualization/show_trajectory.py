
import time
from turtle import pd
import numpy as np


from scripts.Visualization.mini_cheetah import MiniCheetah
from scripts.Visualization.animator import Animator

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
    # animator.plot_pos(t, np.array(msg.pos))
    animator.plot_vel(t, np.array(msg.vWorld))
    # animator.plot_joint_speed(t, np.array(msg.qJd))
    animator.plot_joint_angles(t, np.array(msg.qJ))
    animator.plot_contact(t, np.array(msg.contact))    
    animator.plot_defect(t, np.array(msg.defect))
    animator.plot_roll_abad_rate(t, np.array(msg.qJd), np.array(msg.eulrate))
    animator.plot_centroid_ang_momentum(t, np.array(msg.hg))
    animator.plot_centroid_ang_momentum_derivative(t, np.array(msg.dhg))
    animator.show_plots()    
    

lc = lcm.LCM()
subscription = lc.subscribe("visualize_wb_traj", visualize_motion_lcm_handler)

try:
    while True:        
        lc.handle()
except KeyboardInterrupt:
    pass