#----------------------------------------------------------------------------------------------------
#       run.py
#   capable of building, running, and managing all modules in the project
#   requires pexpect python library to run, may not be inlcuded in your base python distribution
#   also requires python 3.6
#----------------------------------------------------------------------------------------------------


import subprocess
import pexpect
import os
import re

def menu(): #simple menu function
    print("Menu:")
    print("     0. Generate Graphs")
    print("     1. Simulate")
    print("     2. Install Hamiltonian")
    print("     3. Uninstall Hamiltonian")
    print("     4. Clean installation")
    print("     5. Quit")


    access = int(input("Make a numerical selection from the above list: "))
    return access


# Installs a hamiltonian for use in simulation
# Hamiltonains must be in the 'hamiltonian' subfolder
# Edits hamiltonian.txt with name and hamiltonians.h with an appropriate include statement
# Requires input of a filename (c header file) where hamiltonian is implemented

def installHam(name):
    wFile = open("hamiltonian/hamiltonians.txt", "a+")
    wFile.seek(0)
    hamName = ""
    for char in name:
        if char == '.':
            break
        hamName += char
    
    installed = False
    
    for line in wFile:
        if line.strip() == hamName:
            print("Hamiltonian has already been installed")
            wFile.close()
            installed = True
            break
    
    
    if installed == False:
        wFile.write(hamName + '\n')
        wFile.close()
        hFile = open("hamiltonian/hamiltonians.h", 'r')
        lookup = "#endif /* hamiltonians_h */"
        linenum = 0
        for num, line in enumerate(hFile, 0):
            if line.strip() == lookup:
                linenum = num
                break
        hFile.seek(0)
        contents = hFile.readlines()
        hFile.close()
        contents.insert(linenum, '#include "' + name + '"' + '\n')
        eFile = open("hamiltonian/hamiltonians.h", "w")
        contents = "".join(contents)
        eFile.write(contents)
        eFile.close

def uninstallHam(name):
    linenum = 0
    found = False
    rFile = open("hamiltonian/hamiltonians.txt", "r")
    for num, line in enumerate(rFile, 0):
        if line.strip() == name:
            linenum = num
            found = True
            break
    if not found:
        print("Error: hamiltonian not found")
        return
    rFile.seek(0)
    contents = rFile.readlines()
    rFile.close()
    contents.pop(linenum)
    oFile = open("hamiltonian/hamiltonians.txt", "w")
    contents = "".join(contents)
    oFile.write(contents)
    oFile.close

    name += ".h"
    lookup = '#include "' + name + '"'
    hFile = open("hamiltonian/hamiltonians.h", "r")
    for num, line in enumerate(hFile, 0):
        if line.strip() == lookup:
            linenum = num
            break
    hFile.seek(0)
    contents = hFile.readlines()
    hFile.close()
    contents.pop(linenum)
    oFile = open("hamiltonian/hamiltonians.h", "w")
    contents = "".join(contents)
    oFile.write(contents)
    oFile.close






# Modifies Simulation.cpp to set it to run the correct hamiltonian
# Can be found on line that reads 'simFunction = ......."
# As an input, takes the name of a hamiltonian as it reads in hamiltonian.txt (assumes associated header file is there)


def setHam(name):
    name += "Ham" #Naming convention
    simFile = open("Simulate/Simulate.cpp" , 'r')
    lookup = "simFunction = "
    linenum = 0
    for num, line in enumerate(simFile, 0):
        if line.strip().startswith(lookup):
            linenum = num
            break
    simFile.seek(0)
    contents = simFile.readlines()
    simFile.close()
    contents[linenum] = '\t' + lookup + name + ";\n"
    oFile = open("Simulate/Simulate.cpp", "w")
    contents = "".join(contents)
    oFile.write(contents)
    oFile.close
    

buildCanRun = False


# Main loop - runs continuisouly until user quits

while True:
    access = menu()

    if access == 0:
        
        
        if buildCanRun == False:
            print("Checking if buildGraphs needs to be built...")
            proc = subprocess.Popen(['make'],
                            stdout=subprocess.PIPE, cwd='buildGraphs')
            exit_code = proc.wait()
            if exit_code == 0:
                buildCanRun = True;
                print("buildGraphs built sucessfully")
            else:
                print("buildGraphs failed:")
                print(proc.communicate()[0].decode('utf-8'))
                print("terminating")
                break

        proc2 = pexpect.spawn('./buildGraphs/buildGraphs')
        proc2.interact()


    elif access == 1:
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

        setHam(choice)
        print("Building simulation....")
        proc = subprocess.Popen(['make'],
                stdout=subprocess.PIPE, cwd='Simulate')

        
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
            proc2 = pexpect.spawn('./Simulate/Simulate')
            proc2.interact()
            
        proc = subprocess.Popen(['make', 'clean'],
                stdout=subprocess.PIPE, cwd='Simulate')

        

        
        
    elif access == 2:
        filename = input("Please enter the name of the header file: ")
        path = "hamiltonian/" + filename
        print(filename)

        try:

            open(path, "r")
            installHam(filename)

        except:
            print("Hamiltonian not found.")
            print("Make sure the spelling of the filename is correct and the file is located in the hamiltonian folder")

    elif access == 3:
        print("Currently installed hamiltonians: ")
        rFile = open("hamiltonian/hamiltonians.txt", "r")
        print(rFile.read())
        choice = input("Which hamiltonian would you like to uninstall? ")
        uninstallHam(choice)
    
    

    elif access == 4:

        proc = subprocess.Popen(['make', 'clean'],
                            stdout=subprocess.PIPE, cwd='buildGraphs')
        proc.communicate()[0].decode('utf-8')

        proc2 = subprocess.Popen(['make', 'clean'],
                    stdout=subprocess.PIPE, cwd='Simulate')

    elif access == 5:
        break

    
    else:
        print("Invalid option")
