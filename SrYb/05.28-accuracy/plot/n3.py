# -*- coding: utf-8 -*-
"""
Plot of a spectrum of two charged dipoles.
"""
import matplotlib.pyplot as plt
from os.path import expanduser
import numpy as np
import sys
sys.path.insert(0, expanduser('~')+'/ions/lib')
from label_lines import *

# latex font 
plt.rcParams.update({'font.size': 22})
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

def get_n3_truncations(data):
    return data[0]['basis_truncation']['n3']

# =============================================================================
# Plotting functions
# =============================================================================
def get_truncation_string(dataset):
#    I skip here n_3/5 as it is the variable of the model
    truncation = dataset[0][0]['basis_truncation']
    n1 = truncation['n1']
    j1 = truncation['j1']
    out_string = f"$n_1 = {n1}$"+", $j_{1/2} = $"+f"${j1}$"
    return out_string

def get_parameters_string(dataset):
    omega_z = get_omega_z(dataset[0])
    omega_rho = get_omega_rho(dataset[0])
    out_string = f"$\omega_z = {omega_z}$ MHz, "+\
    f"$\omega_\\rho = {omega_rho}$ MHz"
    return out_string


def change_of_energy_levels(dataset, lvl=10):
    ns = []
    spectra = []
    
    energy0 = dataset[0][1][0:lvl]
    for data in dataset:
        ns += [ get_n3_truncations(data) ]
        omegaz = get_omega_z(data)
        omega1 = omegaz*np.sqrt(3)
        omega1_2pi = omega1/2/np.pi
        MHzTOkHz = 1e3
        spectra += [ [(data[1][i] - energy0[i])*omega1_2pi*MHzTOkHz
                      for i in range(lvl)] ]
    
    # present result
    for i in range(lvl):
        plt.plot(ns, [ sp[i] for sp in spectra ], label=f"{i}" )
    plt.title("Energy levels of SrYb, relative changes\n"+
              get_parameters_string(dataset)+
              ", "+get_truncation_string(dataset))
    plt.xlabel("$n_3$ truncation")
    plt.ylabel("$\\frac{E_{i}(5) - E_i(n_3)}{2\pi\hbar}$ (kHz)")
    labelLines(plt.gca().get_lines(),zorder=2.5)
    return 0


# =============================================================================
# Main
# =============================================================================
def main():    
    home = expanduser("~")
    path = home+'/ions/SrYb/05.28-accuracy/'
    ns = np.arange(4)+5
    filenamesA = ['SrYb_n1_34_n35_'+str(n)+'_j1.out' for n in ns]
    dataset = get_dataset(filenamesA, path)
    dataset.sort(key=get_dipole)
    change_of_energy_levels(dataset, 10)

    return 0
    
if __name__ == '__main__':
    main()