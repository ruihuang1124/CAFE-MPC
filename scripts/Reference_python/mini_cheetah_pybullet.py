
from tkinter.messagebox import NO
import numpy as np    
from pybullet_utils.bullet_client import BulletClient
import pybullet
import pybullet_data as pd

import math


HIP_OFFSETS = np.array([[0.19, -0.049, 0],
                          [0.19, 0.049, 0],
                          [-0.19, -0.049, 0],                          
                          [-0.19, 0.049, 0]])

# Rotation direction of Knee and hip follows mini_cheetah/mini_cheetah.urdf
# Warning: this convention is opposite the convention in mini_cheetah_simple_correctedInertia.urdf
DEFAULT_JOINT_POSE = [0, -0.8, 1.6, 
                        0, -0.8, 1.6, 
                        0, -0.8, 1.6, 
                        0, -0.8, 1.6]   

class MiniCheetah:     
    INIT_POS = [0,0,0.25]
    INIT_QUAT = [0,0,0,1]
    EE_ID = [3, 7, 11, 15] # FR, FL, HR, HL
                           # Note: 
                           # This convention is different in mini_cheetah_simple_correctedInertia.urdf, which is FL, FR, HL, HR
    JOINT_DAMPING = [0.1, 0.05, 0.01,
                     0.1, 0.05, 0.01,
                     0.1, 0.05, 0.01,
                     0.1, 0.05, 0.01]  
    def __init__(self, urdf_file=None):        
        # damping used by numerical IK solver
        self.pb = BulletClient(pybullet.DIRECT)
        self.pb.setAdditionalSearchPath(pd.getDataPath())
        self.urdf_file = 'mini_cheetah/mini_cheetah.urdf'   # urdf file stored in pybullet_data
        self.robot = self.pb.loadURDF(self.urdf_file)
        self.set_pose(self.INIT_POS, self.INIT_QUAT, DEFAULT_JOINT_POSE)

    def ik(self, 
           base_pos = None, 
           base_eul = None, 
           foot_pos = None):
        """
            Compute inverse kinmatics
            args: 
                base_pos: base position nadrray of dim 3
                base_quat: base orientation in quaterion
                foot_pos: foot positions in world frame ndarray of dim 12            
        """
        # reset base position and orientation to solve IK in world frame
        base_quat = self.eul2quat(base_eul)
        self.pb.resetBasePositionAndOrientation(self.robot, base_pos, base_quat)
        joint_lim_low, joint_lim_high = self.get_joint_limits()
        joint_angles = self.pb.calculateInverseKinematics2(self.robot, self.EE_ID,
                                                    foot_pos.reshape([4,3]),
                                                    jointDamping = self.JOINT_DAMPING,
                                                    lowerLimits=joint_lim_low,
                                                    upperLimits=joint_lim_high,
                                                    solver=0)
        joint_angles = np.array(joint_angles)
        return joint_angles.reshape(12)

    def fk(self,
           base_pos = None,
           base_eul = None,
           leg = None, qleg = None):
        pf_in_hip = self.leg_fk(leg, qleg)
        pf_in_body = HIP_OFFSETS[leg] + pf_in_hip
        quat = self.eul2quat(base_eul)
        rot = np.reshape(self.pb.getMatrixFromQuaternion(quat), (3,3))
        pf_in_world = np.matmul(rot, pf_in_body) + base_pos

        return pf_in_world
    
    def getSideSign(self,leg_id):
        """Get if the leg is on the right (-1) or the left (+) of the robot"""
        sideSigns = [-1,1,-1,1]
        return sideSigns[leg_id]

    def leg_fk(self, leg_id, leg_angles):
      l1 = 0.062
      l2 = 0.209
      l3 = 0.195
      l4 = 0.004
      sideSign = self.getSideSign(leg_id)

      s1 = math.sin(leg_angles[0])
      s2 = math.sin(leg_angles[1])
      s3 = math.sin(leg_angles[2])

      c1 = math.cos(leg_angles[0])
      c2 = math.cos(leg_angles[1])
      c3 = math.cos(leg_angles[2])

      c23 = c2*c3 - s2*s3
      s23 = s2*c3 + c2*s3

      p = np.zeros(3)
      p[0] = l3 * s23 + l2 * s2
      p[1] = (l1+l4) * sideSign * c1 + l3 * (s1 * c23) + l2 * c2 * s1
      p[2] = (l1+l4) * sideSign * s1 - l3 * (c1 * c23) - l2 * c1 * c2
      return p.copy()

    def eul2quat(self, eul):
        rpy = np.array([eul[2], eul[1], eul[0]])
        return np.array(self.pb.getQuaternionFromEuler(rpy))    
    
    def print_link_jnt_info(self):
        num_joints = self.pb.getNumJoints(self.robot)
        for j in range(num_joints):
            jnt_info = self.pb.getJointInfo(self.robot, j)
            jnt_name = jnt_info[1]
            print("joint {} name".format(j), jnt_name)

    def get_joint_limits(self):
        num_joints = self.pb.getNumJoints(self.robot)
        joint_limit_low = []
        joint_limit_high = []

        for i in range(num_joints):
            joint_info = self.pb.getJointInfo(self.robot, i)
            joint_type = joint_info[2]

            if (joint_type == self.pb.JOINT_PRISMATIC or joint_type == self.pb.JOINT_REVOLUTE):
                joint_limit_low.append(joint_info[8])
                joint_limit_high.append(joint_info[9])

        return joint_limit_low, joint_limit_high  
    
    def set_pose(self, root_pos, root_quat, jnt_pose):
        num_joints = self.pb.getNumJoints(self.robot) 
        self.pb.resetBasePositionAndOrientation(self.robot, root_pos, root_quat)
        for j in range(num_joints):
            j_info = self.pb.getJointInfo(self.robot, j)
            j_state = self.pb.getJointStateMultiDof(self.robot, j)

            j_pose_idx = j_info[3]
            j_pose_size = len(j_state[0])
            j_vel_size = len(j_state[1])

            if (j_pose_size > 0):
                j_pose_idx = j_pose_idx - 7 # shift joint index by 7 to exclude base position and orientation
                j_pose = jnt_pose[j_pose_idx:(j_pose_idx + j_pose_size)]
                j_vel = np.zeros(j_vel_size)
                self.pb.resetJointStateMultiDof(self.robot, j, j_pose, j_vel)