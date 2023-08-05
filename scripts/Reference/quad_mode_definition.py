import string
import numpy as np

QUAD_MODES = ["Fly", # Leg definition -> [FR, FL, HR, HL]
         "HL", 
         "HR", "HR-HL",
         "FL", "FL-HL", "FL-HR", "FL-HR-HL",
         "FR", "FR-HL", "FR-HR", "FR-HR-HL", "FR-FL", "FR-FL-HL", "FR-FL-HR", 
         "Stance"]

def modeNumber2stanceLegs(mode_num:int) -> np.array:
    numstr = np.binary_repr(mode_num, width=4)
    contact = np.array(list(numstr), dtype=int)
    return np.flip(contact)

def modeNumber2string(mode_num:int) -> string:
    for count, item in enumerate(QUAD_MODES):
        if mode_num == count:
            return item    

def stanceLegs2modeNumer(contact:np.array) -> int:
    model_num = contact.dot(2**np.arange(contact.size)[::-1])
    return model_num

def stanceLegs2string(contact:np.array) -> string:
    contact_flipped = np.flip(contact)
    mode_num = stanceLegs2modeNumer(contact_flipped)
    return modeNumber2string(mode_num)

def string2modeNumber(mode_name:string) -> int:
    for count, item in enumerate(QUAD_MODES):
        if mode_name == item:
            return count

def stringSeq2modeNumSeq(string_seq : list):
    modeNumseq = []
    for modeString in string_seq:
        modeNumseq.append(string2modeNumber(modeString))
    
    return np.asarray(modeNumseq)