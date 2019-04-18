# -*- coding: utf-8 -*-
"""
Plot of a spectrum of two charged dipoles.
"""
import numpy as np
import matplotlib.pyplot as plt
    

def get_dataset(filenames, path):
    dataset = []
    
    for name in filenames:
        f = open(path+name, "r")
        if f.mode != "r":
            print("Problem with reading "+name+".\n")
            continue
        
        
    return dataset
        
    
def show(dataset):
    return 0

def main():
    path = "/home/pwojcik/ions/results/"
    filenames = ["D0.0.txt", "D0.4.txt", "D0.8.txt", "D1.6.txt", "D3.txt",
                 "D6.0.txt"]
    dataset = get_dataset(filenames, path)
    show(dataset)
    
    return 0
    
if __name__ == '__main__':
    main()