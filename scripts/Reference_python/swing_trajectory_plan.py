import quad_mode_definition as qm
import gait_schedule as gs
from utils import (approx__leq, approx_geq, approx_equal)
import interpolation as interp
import numpy as np

class FootSwingHeightTrajectory():
    def __init__(self):        
        self.h_ = 0.05
        self.zp_ = 0.0
        self.zv_ = 0.0
        self.ps_ = np.array([0.0,0.0,0.0])
        self.pf_ = np.array([0.0,0.0,0.0])
        self.p_ = np.array([0.0,0.0,0.0])
        self.v_ = np.array([0.0,0.0,0.0])
        

    def set_apex_height(self, h):
        self.h_ = h    
    
    def set_start_position(self, ps):
        self.ps_ = ps

    def set_end_position(self, pf):
        self.pf_ = pf

    def get_position(self):
        return self.p_.copy()
    
    def get_velocity(self):
        return self.v_.copy()

    def computeTrajectory(self, phase, swingTime):
        """ Compute foot swing trajectory with a Bezier curver
        Input phase: where are we in the swing [0, 1]
              swingTime: time duration of a swing        
        """
        if approx_equal(phase, 0.0):
            phase = 0.0
        if approx_equal(phase, 1.0):
            phase = 1.0
            
        self.p_ = interp.CubicBezier(self.ps_, self.pf_, phase)
        self.v_ = interp.CubicBezierFirstDerivative(self.ps_, self.pf_, phase)/swingTime

        # use two Bezier curves for the vertial motion
        if phase < 0.5:
            zp = interp.CubicBezier(self.ps_[2], self.pf_[2] + self.h_, phase*2)
            zv = interp.CubicBezierFirstDerivative(self.ps_[2], self.pf_[2]+self.h_, phase*2)/(0.5*swingTime)
        else:
            zp = interp.CubicBezier(self.pf_[2]+self.h_, self.pf_[2], phase*2 - 1)
            zv = interp.CubicBezierFirstDerivative(self.pf_[2]+self.h_, self.pf_[2], phase*2 -1)/(0.5*swingTime)
        
        self.p_[2] = zp
        self.v_[2] = zv
        

class SwingTrajectoryPlanner:
    def __init__(self, gaitSchedule, footPlanner) -> None:
        self.gaitSchedule_ = gaitSchedule
        self.footPlaner_ = footPlanner
        self.footSwingTrajectories_ = [FootSwingHeightTrajectory() for l in range(4)]        

    def setSwingHeight(self, height) -> None:
        for l in range(4):
            self.footSwingTrajectories_[l].set_apex_height(height)
    
    def computeSwingTrajectoryAtTime(self, leg, time) -> None:
        legContactStatus, legSwitchingTimes = self.gaitSchedule_.getLegContactSchedule()
        
        i = self.gaitSchedule_.getLegModeIndexAtTime(leg, time)        
        contactMode = legContactStatus[leg][i]       

         # If stance, warn the user and return 0
        if contactMode:
            return 0   

        # Get the foothold position at the end of this swing
        pf_next = self.footPlaner_.getFootholdLocation(leg, time)

        # Get the foothold position before swing
        pf_prev = self.footPlaner_.getPreviousFootholdLocation(leg, time)

        self.footSwingTrajectories_[leg].set_start_position(pf_prev)
        self.footSwingTrajectories_[leg].set_end_position(pf_next)

        # Calculate swing phase and swing period of time
        swingStartTime = legSwitchingTimes[leg][i]
        swingFinalTime = legSwitchingTimes[leg][i + 1]
        swingPeriod = swingFinalTime - swingStartTime
        if abs(swingPeriod) < 1e-6:
            raise Exception("No proper swing period was found")
        swingPhase = (time - swingStartTime)/swingPeriod

        self.footSwingTrajectories_[leg].computeTrajectory(swingPhase, swingPeriod) 


    def getSwingPosition(self, leg):
        return self.footSwingTrajectories_[leg].get_position()

    def getSwingVelocity(self, leg):
        return self.footSwingTrajectories_[leg].get_velocity()