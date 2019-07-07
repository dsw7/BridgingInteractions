import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mc

size_f = 5  # set figure sizes in inches
save_f = True


# draw circle in polar coordinate system and map to Cartesian cooordinates
# also prepare some dummy MT and SD coordinates
# ===========================================================

radius = 5.0
theta  = np.arange(0, 2 * np.pi, 2 * np.pi / 20)

# these are basically our SF coordinates
x = radius * np.cos(theta)
y = radius * np.sin(theta)
z = np.zeros(x.shape[0])  # generate a third z axis

coordinates_SF = np.column_stack((x, y, z))  # zip() over ndarrays  

MT = np.array([2.2, -2, 0])
SD = np.array([2,  1.5, 0])          
            
 
           
# visualize the raw feature space
# ===========================================================

fig, ax = plt.subplots()

fig.set_figheight(size_f)
fig.set_figwidth(size_f)
ax.set_aspect('equal')

ax.set_xticks([])
ax.set_yticks([])

# SF coordinates
ax.plot(x, y, 'bo', alpha=0.5)
ax.plot(x, y, c='k', lw=0.75)

# SD / MT coordinates
dim_red = np.column_stack((MT, SD))[0:2]
ax.plot(*dim_red, 'ro')
ax.plot(*dim_red, c='r', lw=0.75)
ax.text(*MT[0:2] + 0.4, 'MT', ha='center')
ax.text(*SD[0:2] + 0.4, 'SD', ha='center')

if save_f: 
    plt.savefig('scalene_step1.png', dpi=1000, bbox_inches='tight')

plt.show()



# we now "web" our feature space with SD / SF vectors
# ===========================================================

lines_SD_SF = [[SD[0:2], item[0:2]] for item in coordinates_SF]
lc = mc.LineCollection(lines_SD_SF, lw=1, color='k', alpha=0.25)

fig, ax = plt.subplots()

fig.set_figheight(size_f)
fig.set_figwidth(size_f)
ax.set_aspect('equal')

ax.set_xticks([])
ax.set_yticks([])

# SF coordinates
ax.plot(x, y, 'bo', alpha=0.5)
ax.plot(x, y, c='k', lw=0.75)

# SD / MT coordinates
dim_red = np.column_stack((MT, SD))[0:2]
ax.plot(*dim_red, 'ro')
ax.plot(*dim_red, c='r', lw=0.75)
ax.text(*MT[0:2] + 0.4, 'MT', ha='center')
ax.text(*SD[0:2] + 0.4, 'SD', ha='center')

# add SF / SD lines
ax.add_collection(lc)

if save_f: 
    plt.savefig('scalene_step2.png', dpi=1000, bbox_inches='tight')

plt.show()



# we now "web" our feature space with MT / SF vectors
# ===========================================================

lines_SD_SF = [[SD[0:2], item[0:2]] for item in coordinates_SF]
lines_MT_SF = [[MT[0:2], item[0:2]] for item in coordinates_SF]
lc1 = mc.LineCollection(lines_SD_SF, lw=1, color='k', alpha=0.25)
lc2 = mc.LineCollection(lines_MT_SF, lw=1, color='k', alpha=0.25)

fig, ax = plt.subplots()

fig.set_figheight(size_f)
fig.set_figwidth(size_f)
ax.set_aspect('equal')

ax.set_xticks([])
ax.set_yticks([])

# SF coordinates
ax.plot(x, y, 'bo', alpha=0.5)
ax.plot(x, y, c='k', lw=0.75)

# SD / MT coordinates
dim_red = np.column_stack((MT, SD))[0:2]
ax.plot(*dim_red, 'ro')
ax.plot(*dim_red, c='r', lw=0.75)
ax.text(*MT[0:2] + 0.4, 'MT', ha='center')
ax.text(*SD[0:2] + 0.4, 'SD', ha='center')

ax.add_collection(lc1)  # SD / SF lines
ax.add_collection(lc2)  # MT / SF lines

if save_f: 
    plt.savefig('scalene_step3.png', dpi=1000, bbox_inches='tight')

plt.show()



# now we need to get distances from MT and SD to SF
# ===========================================================

def get_closest(c_SF, c_SD, c_MT):
    """
    I casted this as a function to make this modular for actual study.
    Input
    -----
    c_SF : iterable of surface coordinates obtained from PyMOL
    c_SD : a single R3 set of coordinates of form array([ 2.2, -2. ,  0. ])
    c_MT : a single R3 set of coordinates of form array([2. , 1.5, 0. ])
    
    Output
    ------
    The coordinate pair that yields the minimum distance in the form:
        (SD or MT, SF)
    """
    
    # SD / SF distances
    SD_SF = {}
    for xyz in c_SF:
        distance = np.linalg.norm(c_SD - xyz)
        SD_SF[distance] = (c_SD, xyz)
        
    # MT / SF distances
    MT_SF = {}
    for xyz in c_SF:
        distance = np.linalg.norm(c_MT - xyz)
        MT_SF[distance] = (c_MT, xyz)
     
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
        
pair = get_closest(coordinates_SF, SD, MT)
SF = pair[1]
    
    
# the final plot
# ===========================================================

lines_SD_SF = [[SD[0:2], item[0:2]] for item in coordinates_SF]
lines_MT_SF = [[MT[0:2], item[0:2]] for item in coordinates_SF]
lc1 = mc.LineCollection(lines_SD_SF, lw=0.5, color='k', alpha=0.25)
lc2 = mc.LineCollection(lines_MT_SF, lw=0.5, color='k', alpha=0.25)

fig, ax = plt.subplots()

fig.set_figheight(size_f)
fig.set_figwidth(size_f)
ax.set_aspect('equal')

ax.set_xticks([])
ax.set_yticks([])

# SF coordinates
ax.plot(x, y, 'ko', alpha=0.25)
ax.plot(x, y, c='k', lw=0.75, alpha=0.25)

# SD / MT coordinates
dim_red = np.column_stack((MT, SD))[0:2]
ax.plot(*dim_red, 'ko')
ax.plot(*dim_red, c='k', lw=0.75)
ax.text(*MT[0:2] + 0.4, 'MT', ha='center')
ax.text(*SD[0:2] + 0.4, 'SD', ha='center')


ax.add_collection(lc1)  # all SD / SF lines
ax.add_collection(lc2)  # all MT / SF lines

# the resultant SF coordinate
ax.plot(*SF[0:2], 'ko', lw=0.75)
ax.text(*SF[0:2] + 0.4, 'SF', ha='center')

lc3 = mc.LineCollection([[SF[0:2], MT[0:2]]], lw=1, color='k')
lc4 = mc.LineCollection([[SF[0:2], SD[0:2]]], lw=1, color='k')

ax.add_collection(lc3)  # all MT / SF lines
ax.add_collection(lc4)  # all MT / SF lines

if save_f: 
    plt.savefig('scalene_step4.png', dpi=1000, bbox_inches='tight')

plt.show()
      

