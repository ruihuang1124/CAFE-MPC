import numpy as np
from gait_schedule import GaitSchedule
from body_trajectory_plan import CoMTrajectoryPlanner

# Default foothold locations w.r.t. the CoM
DEFAULT_FOOTHOLDS = [np.array([0.22, -0.10, 0.]), 
                     np.array([0.22, 0.10, 0.]),
                     np.array([-0.18, -0.10, 0.]), 
                     np.array([-0.18, 0.10, 0.])]
KSCALE = 1
class FootholdPlanner:
    def __init__(self, 
                 gaitSchedule : GaitSchedule,
                 coMPlanner : CoMTrajectoryPlanner):
        self.gaitSchedule_ = gaitSchedule
        self.coMPlanner_ = coMPlanner
        self.pf_ = [[],[],[],[]]

    def computeFootholdLocations(self):
        """ 
        @brief: 
               Update the foothold locations for all contact modes of each leg.
               This function only needs to be called ONCE provided that the gaitschedule and the CoMPlan is unchanged.
               Foothold location for swing leg is determined by Raibert Heuristics.
               Foothold location for stance leg remains the same as the previous swing phase
        """
        
        legContactStatus, legSwitchingTimes = self.gaitSchedule_.getLegContactSchedule()
        
        for l in range(4):
            numLegModes = len(legContactStatus[l])
            self.pf_[l] = [DEFAULT_FOOTHOLDS[l].copy() for _ in range(numLegModes)]

            # For swing leg, foothold location is determined by the end of swing
            for i in range(1, numLegModes):
                legContact = legContactStatus[l][i]
                if legContact == 0:                    
                    touchDownTime = legSwitchingTimes[l][i+1]
                    stancePeriod = 0.2 # Default stance time
                    if i < numLegModes - 2:
                        stancePeriod = legSwitchingTimes[l][i+2] - touchDownTime

                    comPos = self.coMPlanner_.getCoMPosition(touchDownTime)
                    comVel = self.coMPlanner_.getCoMVelocity(touchDownTime)
                    x = comPos[0]
                    y = comPos[1]
                    vx = comVel[0]
                    vy = comVel[1]
                    x_offset = min(vx * KSCALE* stancePeriod/2.0, 0.2) + DEFAULT_FOOTHOLDS[l][0]
                    y_offset = min(vy * KSCALE* stancePeriod/2.0, 0.2) + DEFAULT_FOOTHOLDS[l][1]
                    pf_x = x + x_offset
                    pf_y = y + y_offset

                    self.pf_[l][i] = np.array([pf_x, pf_y, 0.0])
                     

            # For stance leg, foothold location is determined by the end of the previous swing phase
            for i in range(1, numLegModes):
                legContact = legContactStatus[l][i]
                if legContact == 1:
                    self.pf_[l][i] = self.pf_[l][i-1]

    def getFootholdLocation(self, leg, t) -> np.array:
        """
        @brief:
                Get the foothold location at time t for one leg
                Remember to call computeFootholdLocations before using this function
        """
        legModeIndex = self.gaitSchedule_.getLegModeIndexAtTime(leg, t)
        return self.pf_[leg][legModeIndex].copy()

    def getPreviousFootholdLocation(self, leg, t) -> np.array:
        """
        @brief: 
                If the requested time is a swing phase, this function retuns
                the foothold location before swing starts                
        """
        legModeIndex = self.gaitSchedule_.getLegModeIndexAtTime(leg, t)
        return self.pf_[leg][legModeIndex-1].copy()
