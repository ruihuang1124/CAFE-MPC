from cmath import pi
from turtle import mode
from gait_schedule import GaitSchedule
from gait_schedule import Trot
from gait_schedule import Bound
from swing_trajectory_plan import SwingTrajectoryPlanner


gaitSchedule = GaitSchedule()
swingPlanner = SwingTrajectoryPlanner(gaitSchedule)
print("Before changing gaitSchedule")
print(swingPlanner.gaitSchedule_.periodicGait_)

gaitSchedule.setPeriodicGait(Trot)
print("After changing gaitSchedule")
print(swingPlanner.gaitSchedule_.periodicGait_)

print("Trot schedule")
print("Mode Sequence: ", Trot.modeSequence)
print("SwitchingTimes: ", Trot.switchingTimes)

print("Bound schedule")
print("Mode sequence: ", Bound.modeSequence)
print("SwitchingTimes: ", Bound.switchingTimes)

modeSchedule= gaitSchedule.getModeSchedule()
print("Oveall mode sequence before change", modeSchedule.getModeSequence())
gaitSchedule.buildModeSchedule_(1.5)
modeSchedule= gaitSchedule.getModeSchedule()
print(len(modeSchedule.getModeSequence()) - len(modeSchedule.getSwitchingTimes()))
print("Oveall mode sequence", modeSchedule.getModeSequence())
