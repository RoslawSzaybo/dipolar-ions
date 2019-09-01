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
def get_truncation_string(dataset):
#    I skip here n_1 as it is the variable of the model
    truncation = dataset[0][0]['basis_truncation']
    n1 = truncation['n1']
    n3 = truncation['n3']
    j1 = truncation['j1']
    out_string = "$n_1 = $"+f" {n1}"
    out_string +=", $n_{3/5}=$"+f" {n3}"
    out_string +=", $j_{1/2} = $"+f" {j1}"
    return out_string

def show_spectrum(dataset, lvl=10, start = 0, fname='test.eps'):
    domain = [] # dipole moments
    spectra = [] # system energy
    
    omega_rho = get_omega_rho(dataset[0])
    omega_z = get_omega_z(dataset[0])
    omega_1 = omega_z*np.sqrt(3)
    omega_1_2pi = omega_1/2/np.pi
    MHzTOkHz = 1e3

    for data in dataset:
        domain += [ get_dipole(data) ]
        spectra += [ [data[1][i]*omega_1_2pi*MHzTOkHz 
                      for i in range(start,lvl)] ]

    # present result
    for i in range(start,lvl):
        plt.plot(domain, [ sp[i-start] for sp in spectra ], label=f"{i}" )

#    plt.legend()

#    plt.title("Energy of the $i$th level of a SrYb$^+$-alike system\n"\
#              f"$\omega_\\rho$={omega_rho} MHz, $\omega_z$={omega_z} MHz,"\
#              "truncation: "+get_truncation_string(dataset))
#    plt.title(f"$\omega_\\rho={omega_rho}$ MHz, $\omega_z={omega_z}$ MHz "\
#              "("+get_truncation_string(dataset)+")")
    
#    plt.text(6.35, 256, "(a)")
    
    plt.xlabel("Dipole moment, $d$ (D)")
    plt.ylabel("$E_i/2\pi\hbar$ (kHz)")
    
    labelLines(plt.gca().get_lines(),zorder=2.5)
    plt.savefig(fname, dpi=300, orientation='portrait', 
                format='eps', transparent=True,
                bbox_inches='tight')
    return 0

def show_one_energy_level_change_together(dataset, lvl=10, start = 0, 
                                          descriptors=["","test.eps",""]):
    tit = descriptors[0]
    fname = descriptors[1]
    fig_name = descriptors[2]
    domain = [] # dipole moments
    spectra = [] # system energy

    omega_rho = get_omega_rho(dataset[0])
    omega_z = get_omega_z(dataset[0])
    omega_1 = omega_z*np.sqrt(3)
    omega_1_2pi = omega_1/2.0/np.pi
    MHzTOHz = 1.0e6

    energy0 = dataset[0][1][0:lvl]
    for data in dataset:
        domain += [ get_dipole(data) ]
        spectra += [ [(data[1][i] - energy0[i])*omega_1_2pi*MHzTOHz 
                      for i in range(start, lvl)] ]
#    print(f"len of domain is {len(domain)}")
#    print(f"len of spectra is {len(spectra)}")

    # present result
    for i in range(start, lvl):
        plt.plot(domain, [ sp[i-start] for sp in spectra ], label=f"{i}",
                 color = f"C{i%10}")


#    plt.title(tit+f"$\omega_\\rho$={omega_rho} MHz, $\omega_z$={omega_z} MHz,"\
#              "truncation: "+get_truncation_string(dataset))
    plt.xlabel("Dipole moment, $d$ (D)")
    plt.ylabel("$\\frac{E_i(D=0) - E_i(D)}{2\pi\hbar}$ (Hz)")

    if fig_name == 'b':
        plt.text(0.3, -1.05, "("+fig_name+")")
    elif fig_name == 'c':
        plt.text(0.3, -0.75, "("+fig_name+")")
    
    labelLines(plt.gca().get_lines(),zorder=2.5)
    plt.savefig(fname, dpi=300, orientation='portrait', 
                format='eps', transparent=True,
                bbox_inches='tight')
    return 0

def show_change_of_both(dataset, excitations, lvl=10, start = 0, 
                        descriptors=["","test.eps",""]):
    tit = descriptors[0]
    fname = descriptors[1]
    fig_name = descriptors[2]
    domain = [] # dipole moments
    spectra = [] # system energy
    
    domain_excited = [] # dipole moments
    spectra_excited = [] # system energy

    omega_rho = get_omega_rho(dataset[0])
    omega_z = get_omega_z(dataset[0])
    omega_1 = omega_z*np.sqrt(3)
    omega_1_2pi = omega_1/2.0/np.pi
    MHzTOHz = 1.0e6

    energy0 = dataset[0][1][0:lvl]
    for data in dataset:
        domain += [ get_dipole(data) ]
        spectra += [ [(data[1][i] - energy0[i])*omega_1_2pi*MHzTOHz 
                      for i in range(start, lvl)] ]
        
    energy0 = excitations[0][1][0:lvl]
    for data in excitations:
        domain_excited += [ get_dipole(data) ]
        spectra_excited += [ [(data[1][i] - energy0[i])*omega_1_2pi*MHzTOHz 
                      for i in range(start, lvl)] ]
#    print(f"len of domain is {len(domain)}")
#    print(f"len of spectra is {len(spectra)}")

    # present result
    for i in range(start, lvl):
        plt.plot(domain, [ sp[i-start] for sp in spectra ], label=f"{i}",
                 color = f"C{i%10}")
        
    for i in range(start, lvl):
        plt.plot(domain_excited, [ sp[i-start] for sp in spectra_excited ], 
                 label=f"e{i}", color = f"C{i%10}")


    plt.title(tit+f"$\omega_\\rho$={omega_rho} MHz, $\omega_z$={omega_z} MHz,"\
              "truncation: "+get_truncation_string(dataset))
    plt.xlabel("Dipole moment, $d$ (D)")
    plt.ylabel("$\\frac{E_i(D=0) - E_i(D)}{2\pi\hbar}$ (Hz)")

    plt.text(0.3, -1.05, "("+fig_name+")")
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
    path = home+"/ions/SrYb/07.23-dipole/SrYb_"
    Ds=["0.0", "0.95", "1.0", "1.05", "2.0", "3.0", "4.0", "4.745", "5.0", 
        "6.0", "7.0"]
    filenames = [ "D"+D+".out" for D in Ds ]
    #
    excitations_dataset = get_excitations(filenames, path)
    excitations_dataset.sort(key=get_dipole)
    #
    dataset = get_dataset(filenames, path)
    dataset.sort(key=get_dipole)
    #
    latex_fonts()
    
    title_figb="Change in the $i$th energy level of a SrYb$^+$-alike system\n"
    descriptors_b = [title_figb, "figr1b.eps", "b"]
    
    title_figc="Change in the $i$th (excited $j$) energy level of a SrYb$^+$-alike system\n"
    descriptors_c = [title_figc, "figr1c.eps", "c"]
    
    title_figd="Change in the $i$th (excited $j$) energy level of a SrYb$^+$-alike system\n"\
    +"both ground and excited\n"
    descriptors_d = [title_figd, "figr1d.eps", "d"]

    # fix the position of the file name like (b) or (c) or alike thing

    show_spectrum(dataset, 10, fname='figr1a.eps')
#    show_one_energy_level_change_together(dataset, 10, 0, descriptors_b)
#    show_one_energy_level_change_together(excitations_dataset, 7, 0, descriptors_c)
#    show_change_of_both(dataset, excitations_dataset, 5, 0, descriptors_d)
#    show_spectrum(excitations_dataset, 17, 11, fname='figr1e.eps')

    return 0

if __name__ == '__main__':
    main()