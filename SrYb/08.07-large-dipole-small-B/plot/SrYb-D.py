# -*- coding: utf-8 -*-
"""
Plot of a spectrum of two charged dipoles.
"""
import numpy as np
import matplotlib.pyplot as plt
from os.path import expanduser
import sys 
sys.path.insert(0, expanduser('~')+'/ions/lib')
from collect_data import *
from label_lines import *
from latex import *

# =============================================================================
# Plotting functions
# =============================================================================

def show_dipole_effect_as_a_fucntion_of_B(dataset_0D, dataset_50D, lvl, start, 
                                          fname='test.eps'):
    domain = [] # rotational constants
    energy_changes = [] # difference between the energy at d=50 and d=0

    omega_rho = get_omega_rho(dataset_0D[0])    
    omega_z = get_omega_z(dataset_0D[0])
    omega_1 = omega_z*np.sqrt(3)
    omega_1_2pi = omega_1/2.0/np.pi
    MHzTOkHz = 1.0e3

    # data_sets is a list of pairs [[descriptors, spectrum]]
    for data0, data50 in zip(dataset_0D, dataset_50D):
        # sanity test to assure that both d=0 and d=50 are computed at same B
        B0 = get_B(data0)
        B50 = get_B(data50)
        if (B0 - B50) > 10.0:
            print("Different values of B!")
            exit(0)
        
        domain += [ B0 ]
        # difference in the enery of ith state as a function of B
        energy_changes += [ [(data50[1][i] - data0[1][i])*omega_1_2pi*MHzTOkHz 
                             for i in range(start, lvl)] ]

    # present result
    for i in range(start, lvl):
        plt.plot(np.log(domain), [ sp[i-start] for sp in energy_changes ], 
                 label=f"{i}", color = f"C{i%10}")

#    plt.title(f"$\omega_\\rho$={omega_rho} MHz, $\omega_z$={omega_z} MHz,"\
#              "truncation: "+get_truncation_string(dataset_0D))
    plt.xlabel("Rotational constant, log(B/MHz)")
    plt.ylabel("$\\frac{E_i(D=0) - E_i(D)}{2\pi\hbar}$ (kHz)")

#    plt.text(0.3, -1.05, "("+fig_name+")")
#    plt.text(6.3, -0.0001, "(c)")
    
    labelLines(plt.gca().get_lines(),zorder=2.5)
    plt.savefig(fname, dpi=300, orientation='portrait', 
                format='eps', transparent=True,
                bbox_inches='tight')
    return 0

# =============================================================================
# Main
# =============================================================================
def main():
    home = expanduser("~")
    path = home+"/ions/SrYb/08.07-large-dipole-small-B/"
    Bs=["503.6513", "400.0", "200.0", "100.0", "50.0", "25.0", "12.0", "5.0"] 
    #MHz
    filenames0 = [ "D0.0_B"+B+".out" for B in Bs ] 
    filenames50 = [ "D50.0_B"+B+".out" for B in Bs ]
    ##
    ##
    excitations_dataset_0D = get_excitations(filenames0, path)
    excitations_dataset_0D.sort(key=get_B) 
    #
    excitations_dataset_50D = get_excitations(filenames50, path)
    excitations_dataset_50D.sort(key=get_B)
    ##
    ##
    dataset_0D = get_dataset(filenames0, path)
    dataset_0D.sort(key=get_B)
    #
    dataset_50D = get_dataset(filenames50, path)
    dataset_50D.sort(key=get_B)
    ##
    ##
    latex_fonts()
    ##
    ##
    show_dipole_effect_as_a_fucntion_of_B(dataset_0D, dataset_50D, 6, 0, 
                                          fname='figr2a.eps')
#    show_dipole_effect_as_a_fucntion_of_B(excitations_dataset_0D, 
#                                          excitations_dataset_50D, 10, 0, 
#                                          fname='figr2b.eps')

    return 0

if __name__ == '__main__':
    main()