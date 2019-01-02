# dsw7@sfu.ca
# function for finding closest surface coordinate to either metal or bridging SD
# see https://github.com/dsw7/BridgingInteractions/tree/master/scalene-triangle/coordinate-to-surface

from numpy import linalg

def get_closest(SF, SD, MT):
    """
    I casted this as a function to make this modular for actual study.
    Input
    -----
    SF : iterable of surface coordinates obtained from PyMOL
    SD : a single R3 set of coordinates of form array([ 2.2, -2. ,  0. ])
    MT : a single R3 set of coordinates of form array([2. , 1.5, 0. ])
    
    Output
    ------
    The coordinate pair that yields the minimum distance in the form:
        (SD or MT, SF)
    """
    
    # SD / SF distances
    SD_SF = {}
    for xyz in SF:
        distance = linalg.norm(SD - xyz[3:6]) # modified this with slice
        SD_SF[distance] = (SD, xyz)
        
    # MT / SF distances
    MT_SF = {}
    for xyz in SF:
        distance = linalg.norm(MT - xyz[3:6]) # modified this with slice
        MT_SF[distance] = (MT, xyz)
     
    # get shortest distances
    min_SD_SF = min(SD_SF)
    min_MT_SF = min(MT_SF)
    
    # return the pair of coords yielding the shortest distance
    if min_SD_SF < min_MT_SF:
        return SD_SF.get(min_SD_SF)
    elif min_SD_SF > min_MT_SF:
        return MT_SF.get(min_MT_SF)
    else:
        return MT_SF.get(min_MT_SF)
        
