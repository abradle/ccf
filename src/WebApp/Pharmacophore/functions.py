from math import sqrt,pow


def eucl_dist(mol_1, mol_2):
    """Function to find the Euclidean distance between two molecules
    Takes 2 objects with attribute x_com, y_com and z_com
    Returns a Euclidean distance"""
    return sqrt(pow((mol_1.x_com-mol_2.x_com),2)+pow((mol_1.y_com-mol_2.y_com),2)+pow((mol_1.z_com-mol_2.z_com),2))
