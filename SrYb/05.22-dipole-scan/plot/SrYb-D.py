# -*- coding: utf-8 -*-
"""
Plot of a spectrum of two charged dipoles.
"""
import numpy as np
import matplotlib.pyplot as plt
from os.path import expanduser
import sys 
sys.path.insert(0, expanduser('~')+'/ions/lib')
from label_lines import *

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

def get_descriptors(f):
    pack = {}
    pack['mass'] = get_value(f, 'mass')
    pack['charge'] = get_value(f, 'charge')
    pack['dipole'] = get_value(f, 'dipole')
    pack['B'] = get_value(f, 'B')
    pack['omega_rho'] = get_value(f, 'omega_rho')
    pack['omega_z'] = get_value(f, 'omega_z')
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
        
        descriptors = get_descriptors(f)            
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
#    I skip here n_1 as it is the variable of the model
    truncation = dataset[0][0]['basis_truncation']
    n1 = truncation['n1']
    n3 = truncation['n3']
    j1 = truncation['j1']
    out_string = "$n_1 = $"+f" {n1}"
    out_string +=", $n_{3/5}=$"+f" {n3}"
    out_string +=", $j_{1/2} = $"+f" {j1}"
    return out_string

def show_one_energy_level_change_together(dataset, lvl=10, start = 0, 
                                          fname="test.eps"):
    domain = [] # dipole moments
    spectra = [] # system energy

    omega_rho = get_omega_rho(dataset[0])
    omega_z = get_omega_z(dataset[0])
    omega_1 = omega_z*np.sqrt(3)
    omega_1_2pi = omega_1/2/np.pi
    MHzTOkHz = 1e3

    energy0 = dataset[0][1][0:lvl]
    for data in dataset:
        domain += [ get_dipole(data) ]
        spectra += [ [(data[1][i] - energy0[i])*omega_1_2pi*MHzTOkHz 
                      for i in range(start, lvl)] ]

    # present result
    for i in range(start, lvl):
        plt.plot(domain, [ sp[i-start] for sp in spectra ], label=f"{i}" )

#
#    plt.title("Change in the $i$th energy level of a SrYb$^+$-alike system\n"\
#              f"$\omega_\\rho$={omega_rho} MHz, $\omega_z$={omega_z} MHz,"\
#              "truncation: "+get_truncation_string(dataset))
    plt.xlabel("Dipole moment, $d$ (D)")
    plt.ylabel("$(E_i(D=0) - E_i(D))/2\pi\hbar$ (kHz)")

    plt.text(6.3, 8, "(b)")
    
    labelLines(plt.gca().get_lines(),zorder=2.5)
    plt.savefig(fname, dpi=300, orientation='portrait', 
                format='eps', transparent=True,
                bbox_inches='tight')
    return 0

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
    
    plt.text(0.0, 480, "(a)")
    
    plt.xlabel("Dipole moment, $d$ (D)")
    plt.ylabel("$E_i/2\pi\hbar$ (kHz)")
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
    path = home+"/ions/SrYb/05.22-dipole-scan/"
    Ds=["0.0", "1.0", "2.0", "3.0", "4.0", "4.745", "5.0", "6.0", "7.0"]
    filenames = [ "D"+D+".out" for D in Ds ]
    dataset = get_dataset(filenames, path)
    dataset.sort(key=get_dipole)
#    show_one_energy_level(dataset,5)
#    show_one_energy_level_change_together(dataset, 5, fname="fig1b.eps")
    show_spectrum(dataset, 11, fname='fig1a.eps')

    return 0

if __name__ == '__main__':
    main()