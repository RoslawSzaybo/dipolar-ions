# -*- coding: utf-8 -*-
"""
Plot of a spectrum of two charged dipoles.
"""
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
from label_lines import *
import numpy as np

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
    truncation = dataset[0][0]['basis_truncation']
    n1 = truncation['n1']
    n3 = truncation['n3']
    j1 = truncation['j1']
    out_string = f"$n_1$={n1}"+", $n_{3/5}$ = "+f"{n3}, "\
    "$j_{1/2}$"+f"={j1}"
    return out_string

def show_spectrum(dataset, m=100, n=0):    
    for data in dataset:
        omega = get_omega_z(data)
        omega_1 = np.sqrt(3)*omega
        # collect spectrum, avoid reaching for elements that are not accesible
        if m > len(data[1]):
            m = len(data[1])
        if n<0 or n >= m:
            print("Incorrect n - No of the smallest eigenvalue")
            n = m-1
        spectrum = [ d*omega_1 for d in data[1][n:m] ] # in MHz
        # present result
        plt.scatter([omega]*(m-n), spectrum, color='k')
        plt.title(f"Spectrum: eigenvalues {n} through {m-1}")
    return 0

# Presents how energy of the lvl-th energy level changes with one of the system
# parameters
def show_energy_level(dataset, lvl=0):  
    domain = [] # omega_z
    spectrum = [] # energy

    for data in dataset:
        omega_z = get_omega_z(data)
        domain += [ omega_z ]
        # the selected energy
        # omega_1 = omega_z * sqrt(3)
        # program outpus in omega_1
        spectrum += [ data[1][lvl]*omega_z*np.sqrt(3) ] # in MHz
    # present result
    plt.scatter(domain, spectrum, color='k')
    plt.title(f"Energy level {lvl}")
    return 0

# presents  energies of the first `lvl` excited states
def show_one_energy_level_change_together(dataset, lvl=10):
    domain = [] # omega_z
    spectra = [] # energy
    omega_rho = get_omega_rho(dataset[0])
    for data in dataset:
        omega_z = get_omega_z(data)
        domain += [ omega_z ]
        # omega_1 = omega_z * sqrt(3)
        # program outpus in omega_1
        spectra += [ [data[1][i]*omega_z*np.sqrt(3)  for i in range(lvl)] ]
    
    # pic
    for i in range(lvl):
        plt.plot(domain, [ sp[i] for sp in spectra ], label=f"{i}" )
    plt.title("MgH$^+$ spectrum as a function of $\omega_z$\n"\
              f"$\omega_\\rho$= {omega_rho}MHz, "\
              "truncation: "+get_truncation_string(dataset))
    plt.xlabel("$\omega_z$ (MHz)")
    plt.ylabel("$E$ (MHz)")
    labelLines(plt.gca().get_lines(),zorder=2.5)
    return 0

# =============================================================================
# Main
# =============================================================================
def main():
    path = "/home/pawel/ions/MgH/05.08-omega_z-scan/"
    omegas = ["0.01", "0.02", "0.04", "0.16", "0.18", "0.2"]
    filenames = ['omega_z-'+str(o)+'.out' for o in omegas]
    dataset = get_dataset(filenames, path)
    dataset.sort(key=get_dipole)
#    show_spectrum(dataset)
#    show_energy_level(dataset, 1)
    show_one_energy_level_change_together(dataset, 8)
    
    return 0
    
if __name__ == '__main__':
    main()