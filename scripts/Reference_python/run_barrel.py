
import utils
import numpy as np
from mini_cheetah_pybullet import MiniCheetah
from barrel_roll import BarrelRoll

# Create a pybullet model for ik computation
robot = MiniCheetah()

EE_pos = robot.fk(np.zeros(3), np.zeros(3), 0, np.array([0.0, -1.2, 2.4]))
print(EE_pos)

zd_stand = .1464
zd_air = 0.35

barrel_planner = BarrelRoll()
barrel_planner.buildSchedule()
barrel_planner.setStandingHeight(zd_stand)
barrel_planner.setBarrelRollHeight(zd_air)

dt = 0.01
plan_horizon = barrel_planner.getOverallDuration()
N = round(plan_horizon/dt) + 1

pos_tau, vel_tau = [], []
z_tau, pf_tau, vf_tau = [], [], []
pfoot_tau = []
contact_tau = []
eul_tau, eulrate_tau = [],[]
time = []
jnt_tau = []
for k in range(N):
    t = k * dt
    pos = barrel_planner.getCoMPosition(t)
    vel = barrel_planner.getCoMVelocity(t)     
    eul = barrel_planner.getEulerAngle(t)
    contact = barrel_planner.getContactFlagsAtTime(t)   
    pf = []
    vf = []
    for l in range(4):        
            pf.append(barrel_planner.getFootPosition(l, t))
            vf.append(barrel_planner.getFootVelocity(l, t))
    jnt_pos = robot.ik(pos, eul, np.hstack(pf))
    pos_tau.append(pos)
    vel_tau.append(vel)
    eul_tau.append(eul)
    eulrate_tau.append(np.array([0,0,0]))
    contact_tau.append(contact)
    pf_tau.append(np.hstack(pf))
    vf_tau.append(np.hstack(vf))
    jnt_tau.append(jnt_pos)
    time.append(t)

# utils.write_traj_to_file(time, pos_tau, eul_tau, vel_tau, eulrate_tau, pf_tau, vf_tau, jnt_tau, contact_tau)
utils.publish_trajectory_lcm(pos_tau, eul_tau, vel_tau, eulrate_tau, jnt_tau)

# utils.plot_com_pos(time, pos_tau)
# utils.plot_com_vel(time, vel_tau)
# utils.plot_swing_height(time, z_tau)
# utils.plot_foothold_locations(time, pfoot_tau)
# utils.plot_foot_positions(time, pf_tau)
# utils.plot_footPosition_and_CoM(pf_tau, pos_tau)
# utils.animate_footPositions_and_CoM(0, pf_tau, pos_tau)
# utils.plot_jnt_position(time, jnt_tau, 0)
# utils.plot_eul(time, eul_tau)


