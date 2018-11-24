import pickle
from os import makedirs, path
from tc_python import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

"""
This program is to map the single-phase homogeneity range of an extended phase in a quaternary system and to make a 3D diagram. 

To run this program, one needs
(1) Thermo-Calc Software package, with TC-Python API
(2) A thermodynamic database containing more than 4 components
(3) Knowledge and experiences of Thermo-Calc will be a plus

This program is to map the single-phase homogeneity range of an extended phase in a quaternary
system.

For this purpose, single point calculations on many grids are utilized. By contrast, default ternary phase diagram mapping does 
not work very well, since only edges are determined on each section/diagram, which means that "holes" will be left.
"""

def plot_3d(list_of_x, list_of_y, list_of_z, xlabel, ylabel, zlabel, title):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

# to draw the compositional tetrahedron for a quaternary system
    xv = np.array([[0.0, 0.0, 0.0],[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]])
    ax.plot([xv[0, 0], xv[1, 0]], [xv[0, 1], xv[1, 1]], [xv[0, 2], xv[1, 2]], 'r')
    ax.plot([xv[0, 0], xv[2, 0]], [xv[0, 1], xv[2, 1]], [xv[0, 2], xv[2, 2]], 'r')
    ax.plot([xv[0, 0], xv[3, 0]], [xv[0, 1], xv[3, 1]], [xv[0, 2], xv[3, 2]], 'r')
    ax.plot([xv[1, 0], xv[2, 0]], [xv[1, 1], xv[2, 1]], [xv[1, 2], xv[2, 2]], 'r')
    ax.plot([xv[1, 0], xv[3, 0]], [xv[1, 1], xv[3, 1]], [xv[1, 2], xv[3, 2]], 'r')
    ax.plot([xv[2, 0], xv[3, 0]], [xv[2, 1], xv[3, 1]], [xv[2, 2], xv[3, 2]], 'r')
    ax.grid(True)
    ax.axis('equal')

# title
    ax.text2D(0.15, 0.95, 'Single-phase homogeneity in quaternary', size=16, transform=ax.transAxes)
# label the tetrahedron
    ax.text(0.0, 0.0, 0.0, "{}".format(ele_a), size=16, color='gold')
    ax.text(1.0, 0.0, 0.0, "{}".format(ele_b), size=16, color='gold')
    ax.text(0.0, 1.0, 0.0, "{}".format(ele_c), size=16, color='gold')
    ax.text(0.0, 0.0, 1.0, "{}".format(ele_d), size=16, color='gold')

# to plot the scattered composition points on the single-phase boundary
    ax.scatter(list_of_x, list_of_y, list_of_z, c='blue', marker='+')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.axes.set_xlim(0.0, 1.0)
    ax.axes.set_ylim(0.0, 1.0)
    ax.axes.set_zlim(0.0, 1.0)

    plt.show()

# for calculating an isothermal section, which tries to capure the target phase and its compositions
def calculate_an_iss(DB, A, B, C, D, x_C):
    with TCPython():
        calculation = (SetUp().select_database_and_elements(DB, [A, B, C, D])
            .get_system()
                       .with_single_equilibrium_calculation()
                       .set_condition("T", 673)
                       .set_condition("P", 100000)
                       .set_condition("N", 1)
                       .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(A), 0.1)
                       .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(B), 0.1)
                       .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(C), x_C)
                       )

        max_A = 1 - x_C
        x_A = 0.0
        while x_A <= max_A:
            max_B = 1 - x_C - x_A
            x_B = 0.0
            while x_B <= max_B:
                try:
                    results = (calculation
                               .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(A), x_A)
                               .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(B), x_B)
                               .calculate()
                               )
                    Qp_dgm = results.get_value_of('DGM({})'.format(Qphase))
                    Qp_a = results.get_value_of('X({},{})'.format(Qphase, ele_a))
                    Qp_b = results.get_value_of('X({},{})'.format(Qphase, ele_b))
                    Qp_c = results.get_value_of('X({},{})'.format(Qphase, ele_c))
                    Qp_d = results.get_value_of('X({},{})'.format(Qphase, ele_d))
                    print(x_A, x_B, x_C)
                    if Qp_dgm > -0.001:
                        list_dgm.append(round(Qp_dgm, 4))
                        list_a.append(round(Qp_a, 4))
                        list_b.append(round(Qp_b, 4))
                        list_c.append(round(Qp_c, 4))
                        list_d.append(round(Qp_d, 4))
                except Exception as e:
                    print("Exception at composition {}, {}, {}: {}".format(x_A, x_B, x_C, str(e)))
                x_B += 0.01
            x_A += 0.01

# made a folder if does not exist
def make_subfolder(fd, kw):
    global  subfolder
    subfolder = fd + kw + '/'
    if not path.exists(subfolder):
        makedirs(subfolder)
        print('create folder ' + subfolder)
    else:
        pass

# save the data to a file, which can be directly used for plotting afterwards
def save_file(list_tuple):
    with open(subfolder + 'Qp_griddatapoints.dat', 'wb') as fp:
        pickle.dump(list_tuple, fp)
        print("{} data are saved to a .dat file".format(len(list_tuple[0])))

# read the data file if one chooses to make a plot using the saved data.
def read_file():
    with open(subfolder + 'Qp_griddatapoints.dat', 'rb') as fp:
         contents = pickle.load(fp)
    return contents

def main():
    global list_a, list_b, list_c, list_d, list_dgm
    list_a, list_b, list_c, list_d, list_dgm = [], [], [], [], []
    global myDB, Qphase, ele_a, ele_b, ele_c, ele_d

# specify the database, the phase name, the system (i.e. the four elements)
# note that the following inputs are fictitious. Use the real ones
    myDB = "TCAL6"
    Qphase = "C14_laves"
    ele_a = "Al"
    ele_b = "Cu"
    ele_c = "Mg"
    ele_d = "Zn"

    keyword = "Eta"
    folder = 'E:/myfolder/'
    make_subfolder(folder, keyword)

# decide whether use the existing data or to perform new calculations
    read_saved_data = False  # manual switch controled by the user
    if read_saved_data:
        data = read_file()
        list_a, list_b, list_c, list_d, list_dgm = data[0], data[1], data[2], data[3], data[4]
        print("There are {} sets of data.".format(len(list_a)))
    else:
# in order to map the whole homogeneity range of the target phase, several isothermal sections
# may have to be calculated.
        calculate_an_iss(myDB, ele_a, ele_b, ele_c, ele_d, 0.30)
        save_file((list_a, list_b, list_c, list_d, list_dgm))
        calculate_an_iss(myDB, ele_a, ele_b, ele_c, ele_d, 0.35)
        save_file((list_a, list_b, list_c, list_d, list_dgm))
        calculate_an_iss(myDB, ele_b, ele_c, ele_a, ele_d, 0.46)
        save_file((list_a, list_b, list_c, list_d, list_dgm))
        calculate_an_iss(myDB, ele_a, ele_c, ele_b, ele_d, 0.25)
        save_file((list_a, list_b, list_c, list_d, list_dgm))
        calculate_an_iss(myDB, ele_a, ele_b, ele_d, ele_c, 0.70)
        save_file((list_a, list_b, list_c, list_d, list_dgm))
    print(list_dgm, '\n', list_a, '\n', list_b, '\n', list_c, '\n', list_d)

    plot_3d(list_a, list_b, list_d, 'x({},{})'.format(Qphase,ele_a), 'x({},{})'.format(Qphase,ele_b),
            'x({},{})'.format(Qphase,ele_d), "Composition of {}".format(Qphase))

main()
