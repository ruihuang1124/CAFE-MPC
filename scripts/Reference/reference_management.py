from gait_schedule import GaitSchedule
from body_trajectory_plan import CoMTrajectoryPlanner
from swing_trajectory_plan import SwingTrajectoryPlanner
from foothold_plan import FootholdPlanner

import numpy as np

class ReferenceManager:
    def __init__(self) -> None:
        self.gaitSchedule_ = GaitSchedule()
        self.comPlanner_ = CoMTrajectoryPlanner()        
        self.footholdPlanner_ = FootholdPlanner(self.gaitSchedule_, self.comPlanner_)
        self.swingPlanner_ = SwingTrajectoryPlanner(self.gaitSchedule_, self.footholdPlanner_)
        self.planHorizon_ = 0.5
    
    def setPeriodicGait(self, periodicGait) -> None:
        self.gaitSchedule_.setPeriodicGait(periodicGait)
    
    def setInitialCoMPosition(self, xinit, yinit, zinit) -> None:
        self.comPlanner_.setInitialCondition(xinit, yinit, zinit)

    def setCoMTargetAndTransitionTime(self, vx_des, vy_des, z_des, transition_time) -> None:
        self.comPlanner_.setTargets(vx_des, vy_des, z_des) 
        self.comPlanner_.setTransitionTime(transition_time)

    def setSwingHeight(self, swingHeight):
        self.swingPlanner_.setSwingHeight(swingHeight)
    
    def setPlanningHorizon(self, planHorizon) -> None:
        self.planHorizon_ = planHorizon

    def computeReferenceTrajectoryOnce(self) -> None:
        self.gaitSchedule_.buildSchedule(self.planHorizon_ * 2.0)
        self.footholdPlanner_.computeFootholdLocations()
    
    def getModeSchedule(self):
        return self.gaitSchedule_.getModeSchedule()

    def getContactStatusAtTime(self, t):
        return self.gaitSchedule_.getContactFlagsAtTime(t)
    
    def getCoMPositionAtTime(self, t) -> np.array:
        return self.comPlanner_.getCoMPosition(t)
    
    def getEulerAngleAtTime(self, t) -> np.array:
        return self.comPlanner_.getEuler(t)

    def getCoMVelAtTime(self, t) -> np.array:
        return self.comPlanner_.getCoMVelocity(t)
    
    def getEulerRateAtTime(self, t) -> np.array:
        return self.comPlanner_.getEulRate(t)

    def getSwingFootPositionAtTime(self, leg, t):
        self.swingPlanner_.computeSwingTrajectoryAtTime(leg, t)
        return self.swingPlanner_.getSwingPosition(leg)

    def getSwingFootVelocityAtTime(self, leg, t):
        self.swingPlanner_.computeSwingTrajectoryAtTime(leg, t)
        return self.swingPlanner_.getSwingVelocity(leg)

    def getFootholdLocationAtTime(self, leg, t) -> np.array:
        return self.footholdPlanner_.getFootholdLocation(leg, t)







    

    