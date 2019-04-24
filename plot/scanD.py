# -*- coding: utf-8 -*-
"""
Plot of a spectrum of two charged dipoles.
"""
import numpy as np
#import matplotlib.pyplot as plt

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
                spectrum = eigenvalues.split()
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

        
        
    return dataset
        
    
def show(dataset):
    for d in dataset:
        print(d[0])
    return 0

def main():
    path = "/home/pawel/dipolar-ions/results/"
    filenames = ["D0.0.txt", "D0.4.txt", "D0.8.txt", "D1.6.txt", "D3.txt",
                 "D6.0.txt"]
    filenames = [ "test.txt" ]
    dataset = get_dataset(filenames, path)
    show(dataset)
    
    return 0
    
if __name__ == '__main__':
    main()