from dataclasses import dataclass
import numpy as np
import quad_mode_definition as quad_mode
from utils import (approx__leq, approx_geq)

@dataclass()
class PeriodicGait:
    modeSequence = quad_mode.stringSeq2modeNumSeq(["Stance"])
    switchingTimes = np.array([0.0, 1.0])

Stance = PeriodicGait()
modeSeqStr = ["Stance"]
Stance.modeSequence = quad_mode.stringSeq2modeNumSeq(modeSeqStr)
Stance.switchingTimes = np.array([0.0, 0.05])

Trot = PeriodicGait()
modeSeqStr = ["FL-HR", "FR-HL"]
Trot.modeSequence = quad_mode.stringSeq2modeNumSeq(modeSeqStr)
Trot.switchingTimes = np.array([0.0, 0.25, 0.50])

FlyTrot = PeriodicGait()
modeSeqStr = ["FL-HR", "Fly", "FR-HL", "Fly"]
FlyTrot.modeSequence = quad_mode.stringSeq2modeNumSeq(modeSeqStr)
FlyTrot.switchingTimes = np.array([0.0, 0.15, 0.25, 0.4, 0.5])

Bound = PeriodicGait()
modeSeqStr = ["HR-HL", "Fly", "FR-FL", "Fly"]
Bound.modeSequence = quad_mode.stringSeq2modeNumSeq(modeSeqStr)
Bound.switchingTimes = np.array([0.0, 0.1, 0.2, 0.3, 0.4])

Pace = PeriodicGait()
modeSeqStr = ["FL-HL", "FR-HR"]
Pace.modeSequence = quad_mode.stringSeq2modeNumSeq(modeSeqStr)
Pace.switchingTimes = np.array([0.0, 0.25, 0.50])

FlyPace = PeriodicGait()
modeSeqStr = ["FL-HL", "Fly", "FR-HR", "Fly"]
FlyPace.modeSequence = quad_mode.stringSeq2modeNumSeq(modeSeqStr)
FlyPace.switchingTimes = np.array([0.0, 0.15, 0.25, 0.4, 0.5])

Pronk = PeriodicGait()
modeSeqStr = ["Stance", "Fly"]
Pronk.modeSequence = quad_mode.stringSeq2modeNumSeq(modeSeqStr)
Pronk.switchingTimes = np.array([0.0, 0.1, 0.3])

class ModeSchedule:
    def __init__(self, modeSequence=[], switchingTimes=[]):        
        self.modeSequence_ = modeSequence      # list of size N
        self.switchingTimes_ = switchingTimes     # list of size N + 1       
    
    def getNumberModes(self):
        if (len(self.modeSequence_) > 0) & (len(self.switchingTimes_) != len(self.modeSequence_) + 1):
            raise Exception("Size of swichingTimes has to be 1 greater than that of mode sequence")            
        return len(self.modeSequence_)   

    def getModeIndexAtTime(self, time):
        for i in range(self.getNumberModes()):
            if approx_geq(time, self.switchingTimes_[i]) & approx__leq(time, self.switchingTimes_[i+1]):
                return i
        if i >= self.getNumberModes():
            raise Exception("Requested time out of switching times")

    def getModeSequence(self):
        return self.modeSequence_

    def getSwitchingTimes(self):
        return self.switchingTimes_


class GaitSchedule:
    def __init__(self):         
        # Assuming always starting with a stance gait    
        self.initialGait_ = Stance                  
        self.periodicGait_ = None
        self.modeSchedule_ = None
        self.endGait_ = None

        # Schedule of each independent leg
        self.legContactStatus_ = [[],[],[],[]]  # List of 4 lists, each representing a contact status for one foot
        self.legSwitchingTimes_ = [[],[],[],[]]    # Switching times for each leg    
    
    def setPeriodicGait(self, periodicGait):
        self.periodicGait_ = periodicGait
    
    def setEndGait(self, endGait = None):
        self.endGait_ = endGait
    
    def buildSchedule(self, finalTime) -> None:
        self.buildModeSchedule_(finalTime)
        self.buildLegContactSchedule_()
    
    def buildModeSchedule_(self, finalTime) -> None: 
        modeSequence = self.initialGait_.modeSequence
        switchingTimes = self.initialGait_.switchingTimes

        periodicGaitMode = self.periodicGait_.modeSequence
        periodicGaitTimes = self.periodicGait_.switchingTimes
        numPeriodicModes = len(periodicGaitMode)

        while switchingTimes[-1] < finalTime:
            endTime = switchingTimes[-1]
            for i in range(numPeriodicModes):
                modeSequence = np.append(modeSequence, periodicGaitMode[i])
                switchTime = min(endTime+periodicGaitTimes[i+1], finalTime)
                switchingTimes = np.append(switchingTimes, switchTime)
                if switchTime >= finalTime:
                    break
        
        if self.endGait_ != None:
            modeSequence = np.hstack((modeSequence, self.endGait_.modeSequence))
            switchingTimes = np.hstack((switchingTimes, self.endGait_.switchingTimes[1:] + finalTime))        
        self.modeSchedule_ = ModeSchedule(modeSequence, switchingTimes)
    
    def buildLegContactSchedule_(self) -> None:
        # Initialize schedule of each indepedent leg
        init_mode = self.modeSchedule_.getModeSequence()[0]
        init_switchingTime = self.modeSchedule_.getSwitchingTimes()[0]
        init_contact = quad_mode.modeNumber2stanceLegs(init_mode)
        num_modes = self.modeSchedule_.getNumberModes()

        for l in range(4):
            self.legContactStatus_[l].append(init_contact[l])
            self.legSwitchingTimes_[l].append(init_switchingTime)
        
        for i in range(1,num_modes):
            mode = self.modeSchedule_.getModeSequence()[i]
            contact = quad_mode.modeNumber2stanceLegs(mode)
            switchingTime = self.modeSchedule_.getSwitchingTimes()[i]
            for l in range(4):
                if contact[l] != self.legContactStatus_[l][-1]:
                    self.legContactStatus_[l].append(contact[l])
                    self.legSwitchingTimes_[l].append(switchingTime)
            
        # For each leg, make sure it ends up with the same time
        finalTime = self.modeSchedule_.getSwitchingTimes()[-1]
        for l in range(4):
            self.legSwitchingTimes_[l].append(finalTime)
                

    def getModeSchedule(self):
        if self.modeSchedule_ is None:
            print("Warning: empty modeSchedule. Make sure the function buildModeSchedule is called")
            return ModeSchedule(self.initialGait_.modeSequence, self.initialGait_.switchingTimes)
        return self.modeSchedule_

    def getLegContactSchedule(self):
        if len(self.legContactStatus_[0]) == 0:
            print("Warning: empty legContactSchedule. Make sure the function buildLegContactSchedule is called")
        return self.legContactStatus_, self.legSwitchingTimes_
    
    def getLegModeIndexAtTime(self, leg, t):
        for i in range(len(self.legContactStatus_[leg])):
            if approx_geq(t, self.legSwitchingTimes_[leg][i]) & approx__leq(t, self.legSwitchingTimes_[leg][i+1]):                
                return i
        print("Warning: No leg contact mode is found")
        return 0

    def getContactFlagsAtTime(self, time):
        modeIndex = self.modeSchedule_.getModeIndexAtTime(time)
        modeSequence = self.modeSchedule_.getModeSequence()
        mode = modeSequence[modeIndex]

        contactFlag = quad_mode.modeNumber2stanceLegs(mode)
        return contactFlag



  





