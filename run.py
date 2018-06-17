#----------------------------------------------------------------------------------------------------
#       run.py
#   capable of building, running, and managing all modules in the project
#   requires pexpect python library to run, may not be inlcuded in your base python distribution
#   also requires python 3.6
#----------------------------------------------------------------------------------------------------


import subprocess
import os
import re
import sys
import imp
from resourceFunctions import * #see resourceFunctions.py for the functions that do the heavy lifting

def menu(): #simple menu function. Prints options and then returns the value inputted by the user
    print("Menu:")
    print("     1. Simulate")
    print("     2. Create New Hamiltonian")
    print("     3. Edit hamiltonian")
    print("     4. Install Hamiltonian")
    print("     5. Uninstall Hamiltonian")
    print("     6. Quit")
    
    access = int(input("Make a numerical selection from the above list: "))
    return access



if sys.version_info < (3, 6): #Checks to ensure the the python version is at least 3.6
    sys.tracebacklimit = 0
    version = ".".join(map(str, sys.version_info[:3]))
    print('This code requires Python 3.6. You are currently running version ' + version)
    raise Exception('Please update your python version and try again')

try: #Checks if the pexpect module is installed
    imp.find_module('pexpect')
except:
    sys.tracebacklimit = 0
    print("This software requires the Python module pexpect. It can most easily be installed using pip")
    sys.exit(3)


import pexpect

buildCanRun = False #boolean that is used to check if compilations worked properly


# Main loop - runs continuisouly until user quits

while True:
    access = menu()

    if access == 1: #Monte Carlo simulation
        path = None
        simType = input("Do a Monte-Carlo simulation? (Y/N): ").upper()
        monteC = True
        if simType == 'N':
            monteC = False
            path = "Simulate/BasicSim" #The basic simulation will probably never be used anymore, but it has been left just in case
        else:
            path = "Simulate/Monte-Carlo"
        hFile = open("hamiltonian/hamiltonians.txt", "r")
        tempContents = hFile.readlines()
        hFile.close()
        found = False
        choice = 0
        while not found:
            print("Installed hamiltonians: ") #Lists out the hamiltonians with a number in front of them
            for i, line in enumerate(tempContents):
                line = line.strip()
                print(str(i) + ": " + line)
            choice = int(input("\nPlease make a numerical choice from the list above: "))
            if choice < 0 or choice >= len(tempContents): #Takes in a number corresponding to a hamiltonian
                print("Invalid selection ")
            else:
                found = True
            


        setHam(tempContents[choice].strip(), monteC)#rewrites the appropriate lines in the Monte Carlo simulation, declared in resourceFunctions.py
        print("Building simulation....")
        proc = subprocess.Popen(['make'],
                                stdout=subprocess.PIPE, cwd=path)#attempts to build the simulation

        
        exit_code = proc.wait()#gets the exit code of the compilation

        if exit_code == 0: #a 0 indicates sucess
            simCanRun = True
            print("Simulate module built sucessfully")
        else:
            print("Simulate module failed to build properly:")
            print(proc.communicate()[0].decode('utf-8'))#outputs the errors if it was unable to compile
            print("terminating...")
            break
        if simCanRun:
            if monteC:#Runs the simulation of compilation was successful
                proc2 = pexpect.spawn('./Simulate/Monte-Carlo/MonteCarlo')
                proc2.interact()
            else:
                proc2 = pexpect.spawn('./Simulate/BasicSim/Simulate')
                proc2.interact()

        proc = subprocess.Popen(['make', 'clean'],
                stdout=subprocess.PIPE, cwd=path)

        

        
    elif access == 2:#Calls a function that can create a new hamiltonian
        newHam()    #declared in resourceFunctions.py
    
    elif access == 3: #Allows the user to edit hamiltonians without opening the header file
        print("Currently installed hamiltonians: ")
        rFile = open("hamiltonian/hamiltonians.txt", "r")
        tempContents = rFile.readlines()
        rFile.close()
        choice = 0
        while True:
            for i, line in enumerate(tempContents):
                print(str(i) + ": "  + line.strip())
            choice = int(input("\nWhich hamiltonian would you like to edit? Make a numberical selection: "))
            if choice < 0 or choice > len(tempContents) - 1:
                print("Invalid selection ")
            else:
                break
        
        
        editCalcCode(tempContents[choice].strip())#defined in resourceFunctions.py

    elif access == 4:#installs a hamiltonian that is declared in a header file
        filename = input("Please enter the name of the header file: ")
        path = "hamiltonian/" + filename
        
        try:
            open(path, "r")
            installHam(filename)

        except:
            print("Hamiltonian not found.")
            print("Make sure the spelling of the filename is correct and the file is located in the hamiltonian folder")

    elif access == 5:#uninstalls a currently installed hamiltonian
        print("Currently installed hamiltonians: ")
        rFile = open("hamiltonian/hamiltonians.txt", "r")
        tempContents = rFile.readlines()
        rFile.close()
        choice = 0
        while True:
            for i, line in enumerate(tempContents):
                print(str(i) + ": "  + line.strip())
            choice = int(input("\nWhich hamiltonian would you like to edit? Make a numberical selection: "))
            if choice < 0 or choice > len(tempContents) - 1:
                print("Invalid selection ")
            else:
                break
        
        
        uninstallHam(tempContents[choice].strip())#defined in resourceFunctions.py

    
    elif access == 6:
        break

    
    else:
        print("Invalid option")
