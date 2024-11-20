from dataclasses import dataclass
from gait_schedule import Stance, Trot
from gait_schedule import Bound
from gait_schedule import FlyTrot
from gait_schedule import Pace
from gait_schedule import FlyPace
from gait_schedule import Pronk
from reference_management import ReferenceManager
import utils
import numpy as np
from mini_cheetah_pybullet import MiniCheetah, ArcDog


# Desired Trajectories
xinit, yinit, zinit = 0.0, 0.0, 0.36
vx_des, vy_des, z_des = 0.5, 0.0, 0.36
swingHeight = 0.15

periodic_plan_horizon = 30.0
transition_time = 2.5
dt = 0.01
N = round(periodic_plan_horizon/dt) + 1

# Desired Gait
periodicGait = Trot
endGait = Stance
endGait.switchingTimes = np.array([0.0, 0.15])
N = N + round(endGait.switchingTimes[-1]/dt)


# Setup the planners
reference_planner = ReferenceManager()
reference_planner.setPeriodicGait(periodicGait)
reference_planner.setEndGait(endGait)
reference_planner.setPlanningHorizon(periodic_plan_horizon)
reference_planner.setInitialCoMPosition(xinit, yinit, zinit)
reference_planner.setCoMTargetAndTransitionTime(vx_des, vy_des, z_des, transition_time)
reference_planner.setSwingHeight(swingHeight)
reference_planner.computeReferenceTrajectoryOnce()

# Create a pybullet model for ik computation
# urdf_filename =  "../../urdf/mini_cheetah_simple_correctedInertia.urdf"
urdf_filename =  "../../urdf/arcdog_simple_correctedInertia.urdf"
robot = ArcDog(urdf_file=urdf_filename)

pos_tau, vel_tau = [], []
z_tau, pf_tau, vf_tau = [], [], []
pfoot_tau = []
contact_tau = []
eul_tau, eulrate_tau = [],[]
time = []
jnt_tau = []
jntvel_tau = []
for k in range(N):
    t = k * dt
    pos = reference_planner.getCoMPositionAtTime(t)
    vel = reference_planner.getCoMVelAtTime(t)
    z = np.array([0.0, 0.0, 0.0, 0.0])
    contact = reference_planner.getContactStatusAtTime(t)
    pf = []
    vf = []
    for l in range(4):
        if contact[l] == 0:
            pf.append(reference_planner.getSwingFootPositionAtTime(l, t))
            vf.append(reference_planner.getSwingFootVelocityAtTime(l, t))
        else:
            pf.append(reference_planner.getFootholdLocationAtTime(l, t))
            vf.append(np.array([0,0,0]))
        z[l] = pf[l][2]
    eul = np.array([0,0,0])
    jnt_pos = robot.ik(pos, eul, np.hstack(pf))

    pos_tau.append(pos)
    vel_tau.append(vel)
    eul_tau.append(eul)
    eulrate_tau.append(np.array([0,0,0]))
    jnt_tau.append(jnt_pos)    
    jntvel_tau.append(np.zeros(12))
    pf_tau.append(np.hstack(pf))
    vf_tau.append(np.hstack(vf))
    contact_tau.append(contact)    
    time.append(t)

utils.write_traj_to_file(time, pos_tau, eul_tau, vel_tau, eulrate_tau, 
                         pf_tau, vf_tau, jnt_tau, jntvel_tau, contact_tau)
utils.publish_trajectory_lcm(time, pos_tau, eul_tau, vel_tau, eulrate_tau, 
                             jnt_tau, jntvel_tau, contact_tau)

# utils.plot_com_pos(time, pos_tau)
# utils.plot_com_vel(time, vel_tau)
# utils.plot_swing_height(time, z_tau)
# utils.plot_foothold_locations(time, pfoot_tau)
# utils.plot_foot_positions(time, pf_tau)
# utils.plot_footPosition_and_CoM(pf_tau, pos_tau)
# utils.animate_footPositions_and_CoM(0, pf_tau, pos_tau)
# utils.plot_jnt_position(time, jnt_tau, 0)


