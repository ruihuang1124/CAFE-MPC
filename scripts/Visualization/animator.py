import time
import numpy as np
import pybullet as pb
from pybullet_utils.bullet_client import BulletClient
import pybullet_data as pd
import matplotlib.pyplot as plt
import threading

lock = threading.Lock()

class Animator:
    def __init__(self, robot):   
        self.robot = robot
        self.robot_urdf = robot.urdf_file       
        self.ground_urdf = "plane.urdf" # defautl ground 
        self.robotID = None
        self.ground = None
        self.additional_path = None                
        self.LINK_RGBAs = None

    def initialization(self):
        self.pb = BulletClient(pb.GUI)
        self.pb.configureDebugVisualizer(pb.COV_ENABLE_SINGLE_STEP_RENDERING,1)
        self.pb.configureDebugVisualizer(pb.COV_ENABLE_GUI,0)
        self.pb.setAdditionalSearchPath(pd.getDataPath())
        self.pb.setGravity(0,0,0)

        # self.ground = self.pb.loadURDF(self.ground_urdf)
        self.robotID = self.pb.loadURDF(self.robot_urdf)
        self.pb.loadURDF(self.ground_urdf)
        self.set_pose(self.robot.INIT_POS, self.robot.INIT_QUAT, self.robot.DEFAULT_JOINT_POSE)
        self.pb.configureDebugVisualizer(self.pb.COV_ENABLE_SINGLE_STEP_RENDERING,1)        

        # store the RGBA visual feature for all links (excluding the base)        
        self.LINK_RGBAs = []        
        visualShapes = self.pb.getVisualShapeData(self.robotID)
        for visualShape in visualShapes:
            self.LINK_RGBAs.append(visualShape[7])

    def addGround(self,ground_urdf):
        self.ground_urdf = ground_urdf
    
    def hide_legs(self):
        num_joints = self.pb.getNumJoints(self.robotID) 
        # Floating base is not accouted by joint index in PyBullet. 
        for jointIndex in range(num_joints):
            self.pb.changeVisualShape(self.robotID, jointIndex, rgbaColor=[0.0,0.0,0.0,0.0])

    def recover_legs(self):
        num_joints = self.pb.getNumJoints(self.robotID) 
        for jointIndex in range(num_joints):
            self.pb.changeVisualShape(self.robotID, jointIndex, rgbaColor=self.LINK_RGBAs[jointIndex])


    def render(self):        
        while True:
            self.pb.configureDebugVisualizer(self.pb.COV_ENABLE_RENDERING,1)            
            time.sleep(0.1)

    def update_camera(self):
        base_pos = np.array(self.pb.getBasePositionAndOrientation(self.robotID)[0])
        [yaw, pitch, dist] = self.pb.getDebugVisualizerCamera()[8:11]
        self.pb.resetDebugVisualizerCamera(3, 30.0, 52.0, base_pos)        
        
    
    def show_config(self, root_pos, root_quat, jnt_pose):
        self.set_pose(root_pos, root_quat, jnt_pose)
    
    def show_motion_sequence(self, base_pos_traj, base_quat_traj, jnt_angle_traj):
        for k in range(len(base_pos_traj)):
            self.set_pose(base_pos_traj[k], base_quat_traj[k], jnt_angle_traj[k])
            time.sleep(0.01)
        
    def set_pose(self, root_pos, root_quat, jnt_pose):
        num_joints = self.pb.getNumJoints(self.robotID) 
        self.pb.resetBasePositionAndOrientation(self.robotID, root_pos, root_quat)
        for j in range(num_joints):
            j_info = self.pb.getJointInfo(self.robotID, j)
            j_state = self.pb.getJointStateMultiDof(self.robotID, j)

            j_pose_idx = j_info[3]
            j_pose_size = len(j_state[0])
            j_vel_size = len(j_state[1])

            if (j_pose_size > 0):
                j_pose_idx = j_pose_idx - 7 # shift joint index by 7 to exclude base position and orientation
                j_pose = jnt_pose[j_pose_idx:(j_pose_idx + j_pose_size)]
                j_vel = np.zeros(j_vel_size)
                self.pb.resetJointStateMultiDof(self.robotID, j, j_pose, j_vel)

    def show_plots(self):
        plt.show()

    def plot_pos(self, time, pos):        
        fig, axs = plt.subplots(3, 1)
        axs[0].plot(time, pos[:, 0])
        # axs[0].set_xlabel('time (s)')
        axs[0].set_ylabel('x (m)')
        axs[0].grid()

        axs[1].plot(time, pos[:, 1])
        # axs[1].set_xlabel('time (s)')
        axs[1].set_ylabel('y (m)')
        axs[1].grid()

        axs[2].plot(time, pos[:, 2])
        axs[2].set_xlabel('time (s)')
        axs[2].set_ylabel('z (m)')
        axs[2].grid()
        # plt.show()
    
    def plot_vel(self, time, vel):
        fig, axs = plt.subplots(3, 1)
        axs[0].plot(time, vel[:, 0])
        # axs[0].set_xlabel('time (s)')
        axs[0].set_ylabel('$v_x$ (m/s)')
        axs[0].grid()

        axs[1].plot(time, vel[:, 1])
        # axs[1].set_xlabel('time (s)')
        axs[1].set_ylabel('$v_y$ (m/s)')
        axs[1].grid()

        axs[2].plot(time, vel[:, 2])
        axs[2].set_xlabel('time (s)')
        axs[2].set_ylabel('$v_z$ (m/s)')
        axs[2].grid()
        # plt.show()
    
    def plot_eul(self, time, eul):
        fig, axs = plt.subplots(3, 1)
        axs[0].plot(time, eul[:, 0])
        # axs[0].set_xlabel('time (s)')
        axs[0].set_ylabel('yaw (rad)')
        axs[0].grid()

        axs[1].plot(time, eul[:, 1])
        # axs[1].set_xlabel('time (s)')
        axs[1].set_ylabel('pitch (rad)')
        axs[1].grid()

        axs[2].plot(time, eul[:, 2])
        axs[2].set_xlabel('time (s)')
        axs[2].set_ylabel('roll (rad)')
        axs[2].grid()
        # plt.show()
    
    def plot_eulrate(self, time, eulrate):
        fig, axs = plt.subplots(3, 1)
        axs[0].plot(time, eulrate[:, 0])
        # axs[0].set_xlabel('time (s)')
        axs[0].set_ylabel('yaw rate (rad/s)')
        axs[0].grid()

        axs[1].plot(time, eulrate[:, 1])
        # axs[1].set_xlabel('time (s)')
        axs[1].set_ylabel('pitch rate (rad/s)')
        axs[1].grid()

        axs[2].plot(time, eulrate[:, 2])
        axs[2].set_xlabel('time (s)')
        axs[2].set_ylabel('roll rate (rad/s)')
        axs[2].grid()
        # plt.show()
    
    def plot_centroid_ang_momentum(self, time, hg):
        fig, axs = plt.subplots(3, 1)
        t_ = time[0:-1]
        hg_ = hg[0:-1, :]
        axs[0].plot(t_, hg_[:, 0])
        # axs[0].set_xlabel('time (s)')
        axs[0].set_ylabel('hg_x')
        axs[0].grid()

        axs[1].plot(t_, hg_[:, 1])
        # axs[1].set_xlabel('time (s)')
        axs[1].set_ylabel('hg_y')
        axs[1].grid()

        axs[2].plot(t_, hg_[:, 2])
        axs[2].set_xlabel('time (s)')
        axs[2].set_ylabel('hg_z ')
        axs[2].grid()
        # plt.show()
    
    def plot_centroid_ang_momentum_derivative(self, time, dhg):
        fig, axs = plt.subplots(3, 1)
        time_ = time[0:-1]
        dhg_ = dhg[0:-1, :]
        axs[0].plot(time_, dhg_[:, 0])
        # axs[0].set_xlabel('time (s)')
        axs[0].set_ylabel('dhg_x')
        axs[0].grid()

        axs[1].plot(time_, dhg_[:, 1])
        # axs[1].set_xlabel('time (s)')
        axs[1].set_ylabel('dhg_y')
        axs[1].grid()

        axs[2].plot(time_, dhg_[:, 2])
        axs[2].set_xlabel('time (s)')
        axs[2].set_ylabel('dhg_z ')
        axs[2].grid()
        # plt.show()
    
    def plot_defect(self, time, defects):
        fig, ax = plt.subplots()
        ax.plot(time, defects)
        ax.set_xlabel('time (s)')
        ax.set_ylabel('defect')
        ax.grid()
        # plt.show()

    def plot_roll_abad_rate(self, time, qd, eulrate):        
        _, ax = plt.subplots()
        ax.plot(time, eulrate[:, 2])
        ax.plot(time, qd[:, 0])
        ax.plot(time, qd[:, 3])
        ax.plot(time, qd[:, 6])
        ax.plot(time, qd[:, 9])
        ax.set_xlabel('time (s)')
        ax.set_ylabel('Rate of change (rad/s)')
        ax.grid()
        ax.legend(["Roll", "FL", "FR", "HL", "HR"])


    def plot_joint_angles(self, time, q):
        """ Plot joint angles versus time
        time: time step in Nx1 numpy array 
        q: joint angles in Nxn numpy array 
        """
        LEG_INDEX = {'FL': 0, 'FR': 1, 'HL': 2, 'HR': 3}
        fig, axs = plt.subplots(2, 2)
        for legname, legid in LEG_INDEX.items():
            axs[int(legid/2), legid % 2].plot(time, q[:, 3*legid:3*legid+3])
            axs[int(legid/2), legid % 2].set_xlabel('time (s)')
            axs[int(legid/2), legid % 2].set_ylabel('q (rad)')
            axs[int(legid/2), legid % 2].set_title(legname)
            axs[int(legid/2), legid % 2].grid()
            axs[0, 0].legend(['abad', 'hip', 'knee']) 

    def plot_joint_speed(self, time, qd):
        """ Plot joint torques versus time
        time: time step in Nx1 numpy array 
        q: joint angles in Nxn numpy array 
        """
        LEG_INDEX = {'FL': 0, 'FR': 1, 'HL': 2, 'HR': 3}
        fig, axs = plt.subplots(2, 2)
        for legname, legid in LEG_INDEX.items():
            axs[int(legid/2), legid % 2].plot(time, qd[:, 3*legid:3*legid+3])
            axs[int(legid/2), legid % 2].set_xlabel('time (s)')
            axs[int(legid/2), legid % 2].set_ylabel('qd (rad/s)')
            axs[int(legid/2), legid % 2].set_title(legname)
            axs[int(legid/2), legid % 2].grid()
        axs[0, 0].legend(['abad', 'hip', 'knee']) 

    def plot_joint_torques(self, time, torque):
        """ Plot joint torques versus time
        time: time step in Nx1 numpy array 
        torque: joint torque in Nxn numpy array 
        """
        LEG_INDEX = {'FL': 0, 'FR': 1, 'HL': 2, 'HR': 3}
        fig, axs = plt.subplots(2, 2)
        for legname, legid in LEG_INDEX.items():
            axs[int(legid/2), legid % 2].plot(time, torque[:, 3*legid:3*legid+3])            
            axs[int(legid/2), legid % 2].set_ylabel('torque (Nm)')
            axs[int(legid/2), legid % 2].set_title(legname)
            axs[int(legid/2), legid % 2].grid()
        axs[0, 0].legend(['abad', 'hip', 'knee'])
        axs[1, 0].set_xlabel('time (s)')
        axs[1, 1].set_xlabel('time (s)')
        # plt.show()

    def plot_GRF(self, time , GRF):
        LEG_INDEX = {'FL': 0, 'FR': 1, 'HL': 2, 'HR': 3}
        fig, axs = plt.subplots(2, 2)
        for legname, legid in LEG_INDEX.items():
            axs[int(legid/2), legid % 2].plot(time[0:-2], GRF[0:-2, 3*legid:3*legid+3])
            axs[int(legid/2), legid % 2].set_xlabel('time (s)')
            axs[int(legid/2), legid % 2].set_ylabel('GRF (N)')
            axs[int(legid/2), legid % 2].set_title(legname)
            axs[0, 0].legend(['Fx', 'Fy', 'Fz'])
        # plt.show()        

    def plot_contact(self, time, contact):
        LEG_INDEX = {'FL': 0, 'FR': 1, 'HL': 2, 'HR': 3}
        fig, axs = plt.subplots(2, 2)
        for legname, legid in LEG_INDEX.items():
            axs[int(legid/2), legid % 2].plot(time, contact[:, legid])
            axs[int(legid/2), legid % 2].set_xlabel('time (s)')            
            axs[int(legid/2), legid % 2].set_title(legname)
            axs[int(legid/2), legid % 2].grid()
        # plt.show()

    def plot_feasibility(self, time, feas):
        plt.plot(time, feas)
        # plt.show()