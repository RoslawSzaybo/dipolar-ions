from os.path import expanduser
import matplotlib.pyplot as plt
import sys 
sys.path.insert(0, expanduser('~')+'/ions/lib')
from label_lines import *

# latex font 
plt.rcParams.update({'font.size': 22})
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'text.usetex': True})

# =============================================================================
# Data readout
# =============================================================================
def read_time(line):
    t = float(line.split()[2])
    return t

def get_state(line):
    # the position of numbers in each line of the output is fixed
    n1 = str( int( line[32:36] ) )
    n3 = str( int( line[37:39] ) )
    n5 = str( int( line[40:42] ) )
    j1 = str( int( line[43:44] ) )
    m1 = str( int( line[45:48] ) )
    j2 = str( int( line[49:50] ) )
    m2 = str( int( line[51:54] ) )
    # construc the final string
    state = n1+n3+n5+j1+m1+j2+m2
    return state

def get_amplitude(line):
    # the position of numbers in each line of the output is fixed
    re = float(line[9:19])
    im = float(line[20:30])
    return [re, im]

# =============================================================================
# Extracts all the data from file `filename` and converts it into a dictionary. 
# Keys correspond to versors. Values are lists of pairs (also lists). First 
# element of each pair is time, and the second element is the amplitude sotred 
# as a pair [re, im]. 
# Keys are strings formed is the following manner 
# |n1,n3,n5;j1,m1;j2,m2> corresponds to 'n1n3n5j1m1j2m2'
# =============================================================================
def read_file(filename):
    file = open(filename, 'r')
    if file.mode != 'r':
        print("Error!")
        print(f"The progam cannot read file {filename}")
        exit(0)
    
    states = {}
    
    line = file.readline()
    while ( line != '' ):
        # omit comments
        if (line[0] == '#'):
            line = file.readline()
            continue        
        # look for data entires
        elif (line[0] == 't'):
            t = read_time(line)
            line = file.readline()
            while( line.split()[0] != '...' ):
                state = get_state(line)
                amp = get_amplitude(line)
                if state in states:
                    states[state] += [[t, amp]]
                else:
                    states[state] = [[t,amp]]
                line = file.readline()

        # otherwise print default error message        
        else:
            print("There is a line which is difficult to understand it is:")
            print(f"{line}")
            
        line = file.readline()
    
    file.close()
    return states

# extracts the three states which appeared as the first ones
def get_main_oscillations(states,n):
    top3 = list(states)[0:n]
    limited = {}
    for state in top3:
        limited[state] = states[state]
    return limited

# =============================================================================
# Plotting
# =============================================================================
def get_amp2(cnumber):
    amp2 = cnumber[0]**2 + cnumber[1]**2
    return amp2
    

def fig3a(states, fname = "fig3a.eps", pltname="(a)"):
    if (fname=="fig3a.eps"):
        #interesting moments of the evolution
        plt.axvline(x=507, color="red")
        plt.axvline(x=1250, color="red")
        plt.axvline(x=1730, color="red")
    # evolution amplitudes
    labels = {}
    labels["0001000"] = "$|\\alpha\\rangle$"
    labels["0000000"] = "$|\\beta\\rangle$"
    labels["0001010"] = "$|\gamma\\rangle$"
    labels["0002000"] = "$|\delta\\rangle$"
    labels["0000010"] = "$|\epsilon\\rangle$"
    for state, fringes in states.items():
        amp2s = [ get_amp2(point[1]) for point in fringes]
        t = [point[0] for point in fringes]
        if state in labels:
            plt.plot(t, amp2s, label = labels[state])
        else:
            plt.plot(t, amp2s, label="$|$"+state+"$\\rangle$")
    
    plt.legend(labelspacing=0.2, 
               loc=5,
               handletextpad=0.3,
               handlelength=0.7)
#    labelLines(plt.gca().get_lines(), zorder=2.5)
    plt.xlabel("Propagation time, $t$ ($\mu$s)")
    
    if ( pltname != "noname" ):
        plt.text(0, 0.28, pltname)
    
    plt.savefig(fname, dpi=300, orientation='portrait', 
                format='eps', transparent=True,
                bbox_inches='tight')


def oscillations(states, fname = "test.eps"):
    labels = {}
    labels["0001000"] = "$|\\alpha\\rangle$"
    labels["0000000"] = "$|\\beta\\rangle$"
    labels["0001010"] = "$|\gamma\\rangle$"
    labels["0002000"] = "$|\delta\\rangle$"
    labels["0000010"] = "$|\epsilon\\rangle$"
    for state, fringes in states.items():
        amp2s = [ get_amp2(point[1]) for point in fringes]
        t = [point[0] for point in fringes]
        if state in labels:
            plt.plot(t, amp2s)
        else:
            plt.plot(t, amp2s, label="$|$"+state+"$\\rangle$")
    
    plt.legend(labelspacing=0.2, 
               loc=5,
               handletextpad=0.3,
               handlelength=0.7)
#    labelLines(plt.gca().get_lines(), zorder=2.5)
    plt.xlabel("Propagation time, $t$ ($\mu$s)")
    plt.text(0, 0.25, "(b)")
    
    plt.savefig(fname, dpi=300, orientation='portrait', 
                format='eps', transparent=True,
                bbox_inches='tight')

def main():
    home = expanduser('~')
    longfname = home + '/ions/quench/dat/j1.dat'
    longstates = read_file(longfname)
    longstates2 = get_main_oscillations(longstates, 5)
    fig3a(longstates2, fname = "fig3c.eps", pltname="noname")

    
#    filename = home + '/ions/quench/dat/j1-2500.dat'
#    states = read_file(filename)
##    oscillations(states, fname="fig3b.eps")
#    states2 = get_main_oscillations(states,5)
#    fig3a(states2)
    
    return 0
    
if __name__ == '__main__':
    main()