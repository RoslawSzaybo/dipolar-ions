# -*- coding: utf-8 -*-
"""
Plot of a spectrum of two charged dipoles.
"""
import matplotlib.pyplot as plt
from os.path import expanduser
import sys
sys.path.insert(0, expanduser('~')+'/ions/lib')
from label_lines import *
import numpy as np

# latex font 
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'text.usetex': True})

# =============================================================================
# Universal set of functions which serve to read the program oputput and 
# prepare it for plotting
# =============================================================================
def get_value(file, name):
    for line in file:
        if line[0] == '#':
            content = line.split()
            if content[1] == name+":":
                file.seek(0)
                return float(content[2])

def get_basis_truncation(file):
    truncation = {}
    for line in file:
        if line[0] == '#':
            content = line.split()
            if content[1] == 'Basis' and content[2] == 'truncation:':
                file.seek(0)
                # content[3] = "|n1,n3,n5;j1,j2>"
                data = content[3][1:-1]
                data = data.split(',')
                truncation['n1'] =  int(data[0])
                truncation['n3'] =  int(data[1])
                subdata = data[2].split(';')
                truncation['n5'] =  int(subdata[0])
                truncation['j1'] = int(subdata[1])
                truncation['j2'] = int(data[3])
                return truncation

def get_descriptors(f, name):
    pack = {}
    pack['mass'] = get_value(f, 'mass')
    pack['charge'] = get_value(f, 'charge')
    pack['dipole'] = get_value(f, 'dipole')
    pack['B'] = get_value(f, 'B')
    pack['omega_rho'] = get_value(f, 'omega_rho')
    pack['omega_z'] = float(name[12:17])  
    pack['basis_truncation'] = get_basis_truncation(f)
    
    return pack

def get_spectrum(f):
    for line in f:
        if line[0] == '#' or line [0] == ' ' :
            content = line.split()
            if content[1] == '100' and content[2] == 'smallest' and\
            content[3] == 'eigenvalues':
                eigenvalues = f.readline()
                f.seek(0)
                spectrum = [float(i) for i in eigenvalues.split()]
                return spectrum
                    
def get_dataset(filenames, path):
    dataset = []
    
    for name in filenames:
        f = open(path+name, "r")
        if f.mode != "r":
            print("Problem with readipoleding "+name+".\n")
            continue
        
        descriptors = get_descriptors(f, name)            
        spectrum = get_spectrum(f)
        f.close()
        dataset += [[descriptors,spectrum]]

    # integrity tests
    # checks if all spectrums have the same number of elements
    spectrum_length = len(dataset[0][1])
    for d in dataset:
        if len(d[1]) != spectrum_length:
            print("Error!")
            print("Not all spectra have the same length.")
            print("The different element is:")
            print(d[0])
            exit( 0 )
        
    return dataset

def get_dipole(data):
    return data[0]['dipole']

def get_omega_z(data):
    return data[0]['omega_z']

def get_omega_rho(data):
    return data[0]['omega_rho']

def get_n1_truncations(data):
    return data[0]['basis_truncation']['n1']

# =============================================================================
# Plotting functions
# =============================================================================
def get_truncation_string(dataset):
    truncation = dataset[0][0]['basis_truncation']
    n1 = truncation['n1']
    n3 = truncation['n3']
    j1 = truncation['j1']
    out_string = f"$n_1$= {n1}"+", $n_{3/5}$ = "+f"{n3}, "\
    "$j_{1/2}$"+f"= {j1}"
    return out_string

# presents  energies of the first `lvl` excited states
def show_one_energy_level_change_together(dataset, last=10, first=0, 
                                          fname='test.eps'):
    domain = [] # omega_z
    spectra = [] # energy
    
    MHzTOkHz=1e3
    
    omega_rho = get_omega_rho(dataset[0])
    for data in dataset:
        omega_z = get_omega_z(data)
        domain += [ omega_z ]
        omega_1 = omega_z * np.sqrt(3)
        omega_1_2pi = omega_1/2/np.pi
        # program outpus in omega_1
        spectra += [ [data[1][i+first]*omega_1_2pi*MHzTOkHz  
                      for i in range(last-first)] ]
    
    # pic
    for i in range(last-first):
        plt.plot(domain, [ sp[i] for sp in spectra ], label=f"{i+first}",
                 color = f'C{i+first}')
#    plt.title("SrYb$^+$ spectrum as a function of $\omega_z$\n"\
#              f"$\omega_\\rho$= {omega_rho}MHz, "\
#              "truncation: "+get_truncation_string(dataset))
    plt.xlabel("$\omega_z$ (MHz)")
    plt.ylabel("$E/2\pi\hbar$ (kHz)")
    plt.legend(labelspacing=-0.1, loc=2)
    
    plt.text(0.155, 440, "(b)")
#    labelLines(plt.gca().get_lines(),zorder=2.5)
    plt.savefig(fname, dpi=300, orientation='portrait', 
                format='eps', transparent=True,
                bbox_inches='tight')
    return 0

# =============================================================================
# Main
# =============================================================================
def main():
    home = expanduser("~")
    path = home+"/ions/SrYb/05.30-omega_z/"
    omegas = ["0.120", "0.122", "0.124", "0.126", "0.128", "0.130", "0.132", 
              "0.134", "0.136", "0.138", "0.140", "0.142", "0.146", "0.148", 
              "0.150", "0.152", "0.154", "0.156", "0.158"]
    filenames = ['SrYb-omega_z'+str(o)+'.out' for o in omegas]
    dataset = get_dataset(filenames, path)
    dataset.sort(key=get_dipole)
    show_one_energy_level_change_together(dataset, 9, 6, fname="fig2b.eps")
    
    return 0
    
if __name__ == '__main__':
    main()