
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
animator.initialization()


def visualize_motion_lcm_handler(channel, data):
    print("received visualization lcm message")    
    msg = wbTraj_lcmt.decode(data)

    print("Size of the received trajectory: ", msg.sz)    
    
    # t = np.array(msg.time) 
    # animator.plot_eul(t, np.array(msg.eul))
    # animator.plot_eulrate(t, np.array(msg.eulrate))
    # animator.plot_pos(t, np.array(msg.pos))
    # animator.plot_vel(t, np.array(msg.vWorld))
    # animator.plot_joint_speed(t, np.array(msg.qJd))    
    # animator.plot_defect(t, np.array(msg.defect))
    # animator.plot_centroid_ang_momentum(t, np.array(msg.hg))
    # animator.plot_centroid_ang_momentum_derivative(t, np.array(msg.dhg))
    # animator.show_plots()
    for k in range(msg.sz):
        eul_k = np.array(msg.eul[k])
        rpy_k = eul_k[[2,1,0]]
        quat_k = np.array(robot.pb.getQuaternionFromEuler(rpy_k))
        pos_k = np.array(msg.pos[k])             
        qJ_k = np.array(msg.qJ[k]) 

        animator.set_pose(pos_k, quat_k, qJ_k)        
        if k > msg.wb_sz -1:
            animator.hide_legs()
        else:
            animator.recover_legs()
        time.sleep(0.001)
    

lc = lcm.LCM()
subscription = lc.subscribe("visualize_wb_traj", visualize_motion_lcm_handler)

render_thread = threading.Thread(target=animator.render)
render_thread.start()
try:
    while True:        
        lc.handle()
except KeyboardInterrupt:
    pass