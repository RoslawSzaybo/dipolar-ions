# -*- coding: utf-8 -*-
"""
Plot of a spectrum of two charged dipoles.
"""
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join

def get_value(file, name):
    for line in file:
        if line[0] == '#':
            content = line.split()
            if content[1] == name+":":
                file.seek(0)
                return content[2]

def get_basis_truncation(file):
    truncation = {}
    for line in file:
        if line[0] == '#':
            content = line.split()
            if content[1] == 'Basis' and content[2] == 'truncation:':
                file.seek(0)
                # content[3] = "|n1,n3,n5,j1,j2>"
                data = content[3][1:-1]
                data = data.split(',')
                truncation['n1'] =  int(data[0])
                truncation['n3'] =  int(data[1])
                truncation['n5'] =  int(data[2])
                truncation['j1'] = int(data[3])
                truncation['j2'] = int(data[4])
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
            print("Problem with reading "+name+".\n")
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
    return float(data[0]['dipole'])

def get_omega_z(data):
    return float(data[0]['omega_z'])

"""
`m` is a number of the largest eigenvalue to display.
`n` is a number of the smallest eigenvalue to display.
"""
def show_dipole(dataset, m=100, n=0):    
    for data in dataset:
        dipole = get_dipole(data)
        if m > len(data[1]):
            m = len(data[1])
        if n<0 or n >= m:
            print("Incorrect n - No of the smallest eigenvalue")
            n = m-1
        spectrum = data[1][n:m]
        plt.scatter([dipole]*(m-n), spectrum, 
                    color='k')
        plt.title(f"Spectrum: eigenvalues {n} through {m-1}")
    return 0

"""
m - How many energy levels to display
"""
def show_dipole_lines(dataset, m=10):    
    dipoles = []
    all_spectra = []

    for data in dataset:
        dipoles += [get_dipole(data)]
        all_spectra += [data[1]] 
        
    # transpose all_spectra
    # each row has equal number of elements
    all_spectra = list(map(list, zip(*all_spectra)))
    
    states_no = len(all_spectra)
    m = states_no if m > states_no else m
    
    for i in range(m):
        plt.plot(dipoles, all_spectra[i])
    plt.title(f"First {m} energy levels")
    plt.xlabel("electric dipole moment (D)")
    plt.ylabel("Energy of the level ($\hbar \omega_1$)")
    return 0

"""
m - How many energy levels to display
"""
def show_omega_z_lines(dataset, m=10):    
    omegas = []
    all_spectra = []

    for data in dataset:
        omegas += [get_omega_z(data)]
        all_spectra += [data[1]] 
        
    # transpose all_spectra
    # each row has equal number of elements
    all_spectra = list(map(list, zip(*all_spectra)))
    
    states_no = len(all_spectra)
    m = states_no if m > states_no else m
    
    for i in range(m):
        plt.plot(omegas, all_spectra[i])
    plt.xlim([100,700])
    plt.ylim([5.0,10.0])
    plt.title(f"First {m} energy levels")
    plt.xlabel("$\omega_z$ (MHz)")
    plt.ylabel("Energy of the level ($\hbar \omega_1$)")
    return 0

def main():
    path = "/home/pawel/dipolar-ions/results/04.22-dipole-scan/"
    path = "/home/pwojcik/ions/results/04.22-dipole-scan/"
    path = "/home/pwojcik/ions/results/04.26-omega_z-scan/"
    #filenames = [ "test.txt" ]
    
    filenames = [f for f in listdir(path) if isfile(join(path, f))]
    dataset = get_dataset(filenames, path)
    
    #show_dipole(dataset)
    #dataset.sort(key=get_dipole)
    #show_dipole_lines(dataset)
    
    dataset.sort(key=get_omega_z)
    show_omega_z_lines(dataset)
    
    return 0
    
if __name__ == '__main__':
    main()