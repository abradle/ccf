import numpy as np
import math


def get_dimensions(line):
    #Get lengths of grid, assume angles are 90,90,90
    line = line.split()
    x_dim = line[0]
    y_dim = line[1]
    z_dim = line[2]
    return (x_dim, y_dim, z_dim)


def get_num_points(line):
    #Get number of grid points
    line = line.split()
    x_points = line[0]
    y_points = line[1]
    z_points = line[2]
    return (x_points, y_points, z_points)


def get_offsets(line):
    #Get max and min offsets for xyz
    line = line.split()
    x_min = line[1]
    x_max = line[2]
    y_min = line[3]
    y_max = line[4]
    z_min = line[5]
    z_max = line[6]
    return((x_min,x_max),(y_min,y_max),(z_min,z_max))


def parse_grid(grd_file):
    '''Returns a dictionary with all coordinates of all non-zero grid points. 
The majority of the grd file is 0, therefore these are not kept to reduce the memory usage'''
    header = grd_file[:5]
    dimensions = get_dimensions(grd_file[2])
    grid_points = get_num_points(grd_file[3])
    grid_space = float(dimensions[0])/float(grid_points[0])
    offsets = get_offsets(grd_file[4])
    x_min = float(offsets[0][0])
    x = x_min
    x_max = float(offsets[0][1])
    y_min = float(offsets[1][0])
    y = y_min
    y_max = float(offsets[1][1])
    z_min = float(offsets[2][0])
    z = z_min
    z_max = float(offsets[2][1])
    access_dict = {}
    centre_points = []
    scores = []
    for line in grd_file[5:]:
        line = line.strip()
        line = float(line)
        if x > x_max:
            x = x_min
            y = y + 1
        if y > y_max:
            y = y_min
            z = z + 1
        if z > z_max:
            print "Should be EOF"
        if line != 0:
            access_dict[(float(x*grid_space),float(y*grid_space),float(z*grid_space))] = line
        x = x + 1
    return access_dict, grid_space


def get_prop(x, y, z, grid):
    '''Looks up propensity of a coordinate. Grid points with 0 propensity were not stored,
so return 0 if KeyError'''
    try:
        return float(grid[(x,y,z)])
    except KeyError:
        return 0


def interpolate(grid_space, grid_dict, x, y, z):
    """Interpoloate a point based on grid spacings"""
    # Get the coords
    x = float(x)
    y = float(y)
    z = float(z)
    # Get the x,y and z positions either side
    x0 = (math.floor(x / grid_space)) * grid_space
    x1 = (math.floor(x / grid_space) + 1) * grid_space
    y0 = (math.floor(y / grid_space)) * grid_space
    y1 = (math.floor(y / grid_space) + 1) * grid_space
    z0 = (math.floor(z / grid_space)) * grid_space
    z1 = (math.floor(z / grid_space) + 1) * grid_space
    # Get the distances for each of these points
    xd = (x - x0) / (x1 - x0)
    yd = (y - y0) / (y1 - y0)
    zd = (z - z0) / (z1 - z0)
    # Now interpolate
    c00 = get_prop(x0, y0, z0, grid_dict) * (1 - xd) + get_prop(x1, y0, z0, grid_dict) * xd
    c10 = get_prop(x0, y1, z0, grid_dict) * (1 - xd) + get_prop(x1, y1, z0, grid_dict) * xd
    c01 = get_prop(x0, y0, z1, grid_dict) * (1 - xd) + get_prop(x1, y0, z1, grid_dict) * xd
    c11 = get_prop(x0, y1, z1, grid_dict) * (1 - xd) + get_prop(x1, y1, z1, grid_dict) * xd
    # Now combine these
    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd
    # Now combine these again
    c = c0 * (1 - zd) + c1 * zd
    return c


def parse_grid_file(grid_file):
    with open(grid_file, 'r') as grd_file:
        grd_file = grd_file.readlines()
        grid_dict, grid_space = parse_grid(grd_file)
        return grid_space, grid_dict

'''
Pseudo Chemoinformatics toolkit code below to give you an idea
#'''
#
#mol = MolFromMol2File(filepath)
#
#interaction_smarts = ['Blah','blah2','etc']
#
#for smarts in interaction_smarts:
#    match = GetMatch(mol,smarts)
#    match_propensities = []
#    for atom in match:
#        coordinates = atom.coordinates
#        propensity = interpolate(grid_space, coordinates)
#        match_propensies.append(propensity)
#    # Calculate the geometric mean of the atoms
#    interaction_score = np.exp(np.sum(np.log(match_propensities)))/len(match)

'''
If any of the atoms have a propensity of 0, it suggests a clash.
Using a geometric mean gives the interaction a score of 0 if one atom has a score of 0
'''




