# -*- coding: utf-8 -*-
"""
Plot of a spectrum of two charged dipoles.
"""
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import numpy as np
from matplotlib.ticker import MaxNLocator


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
    out_string = f"$n_3$/$n_5$ = {n1}, $j_1$/$j_2$ = {j1}"
    return out_string

def show_spectrum(dataset, m=100, n=0):    
    for data in dataset:
        n1 = get_n3_truncations(data)
        # collect spectrum, avoid reaching for elements that are not accesible
        if m > len(data[1]):
            m = len(data[1])
        if n<0 or n >= m:
            print("Incorrect n - No of the smallest eigenvalue")
            n = m-1
        spectrum = data[1][n:m]
        # present result
        plt.scatter([n1]*(m-n), spectrum, 
                    color='k')
        plt.title(f"Spectrum: eigenvalues {n} through {m-1}")
    return 0

def show_one_energy_level(dataset, lvl=0):  
    ns = []
    spectrum = []
    if lvl > len(dataset[0][1]) or lvl < 0:
        print("Incorrect energy level.\n")
        exit(0)

    for data in dataset:
        ns += [ get_n3_truncations(data) ]
        # the selected energy
        spectrum += [ data[1][lvl] ]
    
    # present result
    plt.scatter(ns, spectrum, color='k')
    plt.title(f"Energy level {lvl}")
    return 0

def show_one_energy_level_change(dataset, lvl=0):
    ns = []
    spectrum = []
    
    energy0 = dataset[0][1][lvl]
    for data in dataset:
        ns += [ get_n3_truncations(data) ]
        # the selected energy
        spectrum += [ data[1][lvl] - energy0 ]
    
    # present result
    plt.scatter(ns, spectrum, color='k')
    plt.title(f"Energy level No {lvl}, relative change")
    plt.xlabel("$n_3$/$n_5$ truncation ")
    plt.ylabel("$E_{lvl} - E_0$ ($\hbar \omega_0$)")
    return 0

def show_one_energy_level_change_together(dataset, lvl=10):
    ns = []
    spectra = []
    
    energy0 = dataset[0][1][0:lvl]
    for data in dataset:
        ns += [ get_n3_truncations(data) ]
        spectra += [ [data[1][i] - energy0[i] for i in range(lvl)] ]
    
    # present result
    for i in range(lvl):
        plt.plot(ns, [ sp[i] for sp in spectra ], label=f"{i}" )
    plt.title("Energy levels, relative changes\n"+
              get_truncation_string(dataset))
    plt.legend()
    plt.xlabel("$n_3$ truncation")
    plt.ylabel("$E_{lvl} - E_0$ ($\hbar \omega_0$)")
    return 0


# =============================================================================
# Main
# =============================================================================
def main():
    path = "/home/pawel/ions/SrYb/accuracy/"
    ns = np.arange(3)+3
    filenames = ['SrYb_n1_26_n35_'+str(n)+'_j2.out' for n in ns]
    dataset = get_dataset(filenames, path)
    dataset.sort(key=get_dipole)
#    show_spectrum(dataset)
#    show_one_energy_level(dataset, 0)
#    show_one_energy_level_change(dataset, 4)
    show_one_energy_level_change_together(dataset, 10)

    return 0
    
if __name__ == '__main__':
    main()