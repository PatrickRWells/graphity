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
from resourceFunctions import *

if sys.version_info < (3, 6):
    sys.tracebacklimit = 0
    version = ".".join(map(str, sys.version_info[:3]))
    print('This code requires Python 3.6. You are currently running version ' + version)
    raise Exception('Please update your python version and try again')

try:
    imp.find_module('pexpect')
except:
    sys.tracebacklimit = 0
    print("This software requires the Python module pexpect. It can most easily be installed using pip")
    sys.exit(3)


import pexpect

buildCanRun = False


# Main loop - runs continuisouly until user quits

while True:
    access = menu()

    if access == 1:
        path = None
        simType = input("Do a Monte-Carlo simulation? (Y/N): ").upper()
        monteC = True
        if simType == 'N':
            monteC = False
            path = "Simulate/BasicSim"
        else:
            path = "Simulate/Monte-Carlo"
        hFile = open("hamiltonian/hamiltonians.txt", "r")
        found = False
        choice = ""
        while not found:
            print("Installed hamiltonians: ")
            print(hFile.read())
            choice = input("Please enter the hamiltonian you wish to simulate: ")
            hFile.seek(0)
            for line in hFile:
                if line.strip() == choice:
                    found = True
                    hFile.close()
                    break
            if found:
                break
            print("Hamiltonian not found")
            hFile.seek(0)

        setHam(choice, monteC)
        print("Building simulation....")
        proc = subprocess.Popen(['make'],
                                stdout=subprocess.PIPE, cwd=path)

        
        exit_code = proc.wait()

        if exit_code == 0:
            simCanRun = True
            print("Simulate module built sucessfully")
        else:
            print("Simulate module failed to build properly:")
            print(proc.communicate()[0].decode('utf-8'))
            print("terminating")
            break
        if simCanRun:
            if monteC:
                proc2 = pexpect.spawn('./Simulate/Monte-Carlo/MonteCarlo')
                proc2.interact()
            else:
                proc2 = pexpect.spawn('./Simulate/BasicSim/Simulate')
                proc2.interact()

        proc = subprocess.Popen(['make', 'clean'],
                stdout=subprocess.PIPE, cwd=path)

        

        
    elif access == 2:
        newHam()
    
    elif access == 3:
        print("Currently installed hamiltonians: ")
        rFile = open("hamiltonian/hamiltonians.txt", "r")
        print(rFile.read())
        choice = input("Which hamiltonian would you like to edit? ")
        editCalcCode(choice)

    elif access == 4:
        filename = input("Please enter the name of the header file: ")
        path = "hamiltonian/" + filename
        
        try:

            open(path, "r")
            installHam(filename)

        except:
            print("Hamiltonian not found.")
            print("Make sure the spelling of the filename is correct and the file is located in the hamiltonian folder")

    elif access == 5:
        print("Currently installed hamiltonians: ")
        rFile = open("hamiltonian/hamiltonians.txt", "r")
        print(rFile.read())
        choice = input("Which hamiltonian would you like to uninstall? ")
        uninstallHam(choice)
    
    
    elif access == 6:
        break

    
    else:
        print("Invalid option")
