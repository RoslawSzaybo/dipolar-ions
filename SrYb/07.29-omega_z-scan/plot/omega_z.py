# -*- coding: utf-8 -*-
"""
Plot of a spectrum of two charged dipoles.
"""
import matplotlib.pyplot as plt
from os.path import expanduser
import numpy as np
import sys 
sys.path.insert(0, expanduser('~')+'/ions/lib')
from collect_data import *
from label_lines import *
from latex import *

# =============================================================================
# Plotting functions
# =============================================================================
def get_truncation_string(dataset):
    truncation = dataset[0][0]['basis_truncation']
    n1 = truncation['n1']
    n3 = truncation['n3']
    j1 = truncation['j1']
    out_string = f"$n_1=$ {n1}"+", $n_{3/5}$ ="+f" {n3}, "\
    "$j_{1/2}=$"+f" {j1}"
    return out_string

# presents  energies of the first `lvl` excited states
def spectrum(dataset, lvl=10, fname='test.eps'):
    domain = [] # omega_z
    spectra = [] # energy
    
    MHzTOkHz=1e3
    
    for data in dataset:
        omega_z = get_omega_z(data)
        omega_1 = omega_z * np.sqrt(3)
        omega_1_2pi = omega_1 /2/np.pi
        domain += [ omega_z ]
        # program outpus in omega_1
        spectra += [ [data[1][i]*omega_1_2pi*MHzTOkHz 
                      for i in range(lvl)] ]
    
    # pic
    for i in range(lvl):
        plt.plot(domain, [ sp[i] for sp in spectra ], label=f"{i}" )
# =============================================================================
#     omega_rho = get_omega_rho(dataset[0])
#     plt.title("SrYb$^+$ spectrum as a function of $\omega_z$\n"\
#               f"$\omega_\\rho$= {omega_rho}MHz, "\
#               "truncation: "+get_truncation_string(dataset))
# =============================================================================
    plt.xlabel("$\omega_z$ (MHz)")
    plt.ylabel("$E/2\pi\hbar$ (kHz)")
    plt.legend(labelspacing=-0.08, loc=2)
    
    plt.text(0.18, 260, "(a)")
    
# =============================================================================
#     labelLines(lines, reverse=True, zorder=2.5)
# =============================================================================
    plt.savefig(fname, dpi=300, orientation='portrait', 
                format='eps', transparent=True,
                bbox_inches='tight')
    return 0

# blow-up of the interesting avoided-crossing
def zoom_in(dataset, lvl=10, start=0, fname='test.eps'):
    domain = [] # omega_z
    spectra = [] # energy
    
    MHzTOkHz=1e3
    
    for data in dataset:
        omega_z = get_omega_z(data)
        omega_1 = omega_z * np.sqrt(3)
        omega_1_2pi = omega_1/2.0/np.pi
        domain += [ omega_z ]
        # program outpus in omega_1
        spectra += [ [data[1][i+start]*omega_1_2pi*MHzTOkHz 
                      for i in range(lvl-start)] ]
    
    # pic
    for i in range(lvl-start):
        plt.plot(domain, [ sp[i] for sp in spectra ], label=f"{i+start}",
                 color = f'C{i+start}' )
# =============================================================================
#     omega_rho = get_omega_rho(dataset[0])
#     plt.title("SrYb$^+$ spectrum as a function of $\omega_z$\n"\
#               f"$\omega_\\rho$= {omega_rho}MHz, "\
#               "truncation: "+get_truncation_string(dataset))
# =============================================================================
    plt.xlabel("$\omega_z$ (MHz)")
    plt.ylabel("$E/2\pi\hbar$ (kHz)")
    
    plt.ylim(464,465.5)
    plt.xlim(0.1599, 0.1615)

    
#    labelLines(plt.gca().get_lines(), reverse=True, zorder=2.5)
    plt.legend()
    
    plt.savefig(fname, dpi=300, orientation='portrait', 
                format='eps', transparent=True,
                bbox_inches='tight')
    return 0

# =============================================================================
# Main
# =============================================================================
def main():
    home = expanduser("~")
    path = home+"/ions/SrYb/07.29-omega_z-scan/"
    omegas = ["0.01", "0.02", "0.04", "0.06", "0.08", "0.1", "0.105", "0.110", "0.112", 
        "0.114", "0.116", "0.118", "0.12", "0.122", "0.124", "0.126", "0.128", 
        "0.130", "0.132", "0.134", "0.136", "0.138", "0.14", "0.142", "0.144", 
        "0.146", "0.148", "0.150", "0.152", "0.154", "0.156", "0.158", "0.159",
        "0.1595", "0.16", "0.1605", "0.161", "0.1615", 
        "0.162", "0.164", "0.168", "0.170", "0.172", "0.174", "0.176", 
        "0.178", "0.18", "0.182", "0.184", "0.186", "0.188", "0.190", "0.192", 
        "0.194", "0.196", "0.198", "0.2"]
#    float_omegas = [float(o) for o in omegas]
#    float_omegas.sort()
#    omegas = [str(o) for o in float_omegas]
    filenames = ['SrYb_w_z-'+str(o)+'.out' for o in omegas]
    latex_fonts()
    dataset = get_dataset(filenames, path)
#    spectrum(dataset, 10, fname="fig2a.eps")
    zoom_in(dataset, 8, 5, fname="zoom-in.eps")
    
    return 0
    
if __name__ == '__main__':
    main()