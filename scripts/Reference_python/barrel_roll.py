from dataclasses import dataclass
import numpy as np
import math

import quad_mode_definition as quad_mode
from gait_schedule import ModeSchedule
from foothold_plan import DEFAULT_FOOTHOLDS
from interpolation import lerp

# from mini_cheetah_pybullet import DEFAULT_JOINT_POSE
from mini_cheetah_pybullet import MiniCheetah

from utils import (approx__leq, approx_geq)

DEFAULT_JOINT_POSE = [np.array([0, -1.2, 2.4]), 
                      np.array([0, -1.2, 2.4]), 
                      np.array([0, -1.2, 2.4]), 
                      np.array([0, -1.2, 2.4])]  

def projectile_pos(T,t):
    h = 1./8.*9.81*T*T
    a = -4.0*h/(T**2)
    z = a*t*(t-T)
    return z
def projectile_vel(T,t):
    h = 1./8.*9.81*T*T
    a = -4.0*h/(T**2)
    zd = a*(2*t - T)
    return zd

@dataclass()
class BarrelRollGait:
    modeSeqStr = ["Stance", "FL-HL", "Fly", "Stance"] # The first stance is when the left legs start pushing off the ground
    modeSequence = quad_mode.stringSeq2modeNumSeq(modeSeqStr)
    switchingTimes = np.array([0.0, 0.1, 0.16, 0.50, 2.0])

class BarrelRoll:
    def __init__(self) -> None:
        self.modeSeq_ = BarrelRollGait.modeSequence
        self.switchingTimes_ = BarrelRollGait.switchingTimes
        self.modeSchedule_ = ModeSchedule(self.modeSeq_, self.switchingTimes_)
        
        self.t_barrel_start_ = 0.0
        self.t_barrel_end_ = BarrelRollGait.switchingTimes[3]
        self.barrel_dur_ = self.t_barrel_end_ - self.t_barrel_start_
        self.zd_stand_ = 0.0 
        self.zd_barrel_ = 0.0
        self.pCoM_landing_ = np.array([0, 0, self.zd_stand_])

        self.robot_ = MiniCheetah()

        # Schedule of each independent leg
        self.legContactStatus_ = [[],[],[],[]]  # List of 4 lists, each representing a contact status for one foot
        self.legSwitchingTimes_ = [[],[],[],[]]    # Switching times for each leg    
    
    def buildSchedule(self):
        self.buildLegContactSchedule_()

    def setStandingHeight(self, zd):
        self.zd_stand_ = zd
    
    def setBarrelRollHeight(self, zd):
        self.zd_barrel_ = zd
    
    def setCoMLandingPosition(self, pCoM_landing):
        self.pCoM_landing_ = pCoM_landing

    def getOverallDuration(self):
        return self.switchingTimes_[-1]

    def getCoMPosition(self, t):
        # update z position
        p_z = self.getZPosition_(t)

        # update x and y position
        p_xy = self.getXYPositions_(t)

        pCoM = np.hstack((p_xy, p_z))
        return pCoM

    def getCoMVelocity(self, t):
        # Assume zero reference velocity for now
        # Change if needed
        vCoM = np.array([0.0, 0.0, 0.0])
        # if t <= self.barrel_dur_:
        #     vCoM[1] = self.pCoM_landing_[1]/self.barrel_dur_
        # vCoM[2] = self.getZVelocity_(t)
        return vCoM

    def getEulerAngle(self, t):
        yaw_t = 0.0
        pitch_t = 0.0
        roll_t = self.getRollAngle_(t)        
        return np.array([yaw_t, pitch_t, roll_t])

    def getEulerRate(self, t):
        yaw_rate = 0.0
        pitch_rate = 0.0
        # roll_rate = self.getRollRate_(t)
        roll_rate = 0.0
        return np.array([yaw_rate, pitch_rate, roll_rate])
    
    def getFootPosition(self, leg, time):        
        pCoM = self.getCoMPosition(min(time, self.t_barrel_end_))
        if time <= self.t_barrel_end_:
            contactFlags = self.getContactFlagsAtTime(time)
            if contactFlags[leg]:            
                pf = pCoM + DEFAULT_FOOTHOLDS[leg]
                pf[2] = 0.0
            else:
                eul = self.getEulerAngle(time)
                qleg = DEFAULT_JOINT_POSE[leg]
                pf = self.robot_.fk(pCoM, eul, leg, qleg)
        else:
            pf = pCoM + DEFAULT_FOOTHOLDS[leg]
            pf[2] = 0.0
        return pf    
    
    def getFootVelocity(self, leg, time):
        return np.array([0,0,0])
    
    def getContactFlagsAtTime(self, time):
        modeIndex = self.modeSchedule_.getModeIndexAtTime(time)
        modeSequence = self.modeSchedule_.getModeSequence()
        mode = modeSequence[modeIndex]

        contactFlag = quad_mode.modeNumber2stanceLegs(mode)
        return contactFlag

    def getRollAngle_(self, t):
        # Linearly interpolates roll flip
        t = min(t, self.t_barrel_end_)
        roll_t = (t - self.t_barrel_start_)/(self.barrel_dur_) * 2.0 * math.pi     
        return roll_t
    def getRollRate_(self, t):
        roll_rate_t = 0.0
        if (t >=0 and t <= self.t_barrel_end_):
            roll_rate_t = 2.0 * math.pi /(self.barrel_dur_)
        return roll_rate_t
    
    def getZVelocity_(self, t):
        vzt = 0.0
        if (t >=0 and t <= self.t_barrel_end_):
            vzt = projectile_vel(self.barrel_dur_, t)        
        return vzt
    
    def getZPosition_(self, t): 
        # half_barrel_time = self.barrel_dur_/2.0
        # t_barrel_middle = self.t_barrel_start_ + half_barrel_time
        # z_diff = self.zd_barrel_ - self.zd_stand_
        # if (t >= 0.0 and t < t_barrel_middle):
        #     zt = self.zd_stand_ + (t - self.t_barrel_start_)/half_barrel_time * z_diff
        # elif (t >= t_barrel_middle and t <= self.t_barrel_end_):
        #     zt = self.zd_barrel_ - (t - t_barrel_middle)/half_barrel_time * z_diff
        # else:
        #     zt = self.zd_stand_
        if (t >=0 and t <= self.t_barrel_end_):
            zt = projectile_pos(self.barrel_dur_, t) + self.zd_stand_
        else:
            zt = self.zd_stand_
        return zt
    
    def getXYPositions_(self, t):        
        t_span = self.barrel_dur_/1.5
        t_rel = min(t, t_span) - self.t_barrel_start_
        x = lerp(0.0, self.pCoM_landing_[0], t_rel/t_span)
        y = lerp(0.0, self.pCoM_landing_[1], t_rel/t_span)
        return np.array([x,y])
    
    def buildLegContactSchedule_(self):
        # Initialize schedule of each indepedent leg
        init_mode = self.modeSeq_[0]
        init_switchingTime = self.switchingTimes_[0]
        init_contact = quad_mode.modeNumber2stanceLegs(init_mode)
        num_modes = len(self.modeSeq_)

        for l in range(4):
            self.legContactStatus_[l].append(init_contact[l])
            self.legSwitchingTimes_[l].append(init_switchingTime)
        
        for i in range(1,num_modes):
            mode = self.modeSeq_[i]
            contact = quad_mode.modeNumber2stanceLegs(mode)
            switchingTime = self.switchingTimes_[i]
            for l in range(4):
                if contact[l] != self.legContactStatus_[l][-1]:
                    self.legContactStatus_[l].append(contact[l])
                    self.legSwitchingTimes_[l].append(switchingTime)
            
        # For each leg, make sure it ends up with the same time
        finalTime = self.switchingTimes_[-1]
        for l in range(4):
            self.legSwitchingTimes_[l].append(finalTime)

