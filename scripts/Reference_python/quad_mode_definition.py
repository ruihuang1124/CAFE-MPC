from statistics import mode
import string
from turtle import pd
import numpy as np

# QUAD_MODES = ["Fly", # Leg definition -> [FR, FL, HR, HL]
#          "HL", 
#          "HR", "HR-HL",
#          "FL", "FL-HL", "FL-HR", "FL-HR-HL",
#          "FR", "FR-HL", "FR-HR", "FR-HR-HL", "FR-FL", "FR-FL-HL", "FR-FL-HR", 
#          "Stance"]

QUAD_MODES = {
    "Fly": [0, 0, 0, 0],
    "FL" : [1, 0, 0, 0],
    "FR" : [0, 1, 0, 0],
    "HL" : [0, 0, 1, 0],
    "HR" : [0, 0, 0, 1],
    "FR-FL": [1, 1, 0, 0],
    "FR-HR": [0, 1, 0, 1],
    "FR-HL": [0, 1, 1, 0],
    "FL-HL": [1, 0, 1, 0],
    "FL-HR": [1, 0, 0, 1],
    "HR-HL": [0, 0, 1, 1],
    "FL-HR-HL": [1, 0, 1, 1],
    "FR-HR-HL": [0, 1, 1, 1],
    "FR-FL-HL": [1, 1, 1, 0],
    "FR-FL-HR": [1, 1, 0, 1],
    "Stance": [1, 1, 1, 1]
}

def string2stanceLegs(mode_name:string):
    for key, value in QUAD_MODES.items():
        if mode_name == key:
            return np.array(value)
    print("Wrong requested mode name")

def stanceLegs2string(contact:np.array) -> string:
    for key, value in QUAD_MODES.items():
        if np.array_equal(np.array(value), contact):
            return key
    print("Wrong requested contact")
            

def modeNumber2stanceLegs(mode_num:int) -> np.array:    
    numstr = np.binary_repr(mode_num, width=4)
    contact = np.array(list(numstr), dtype=int)
    return np.flip(contact)

def stanceLegs2modeNumer(contact:np.array) -> int:
    mode_num = contact.dot(2**np.arange(contact.size)[::-1])
    return mode_num

def modeNumber2string(mode_num:int) -> string:
    contact = modeNumber2stanceLegs(mode)
    return stanceLegs2string(contact)

def string2modeNumber(mode_name:string) -> int:
    contact = string2stanceLegs(mode_name)
    mode_num = stanceLegs2modeNumer(contact)
    return mode_num

def stringSeq2modeNumSeq(string_seq : list):
    modeNumseq = []
    for modeString in string_seq:
        modeNumseq.append(string2modeNumber(modeString))
    
    return np.asarray(modeNumseq)