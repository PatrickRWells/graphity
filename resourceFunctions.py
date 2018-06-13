#----------------------------------------------------------------------------------------------------
#       resourceFunctions.py
#   contains resources for the run.py script
#----------------------------------------------------------------------------------------------------

import pexpect
import subprocess
import os
import re
import sys
import imp
from time import sleep

def menu(): #simple menu function
    print("Menu:")
    print("     1. Simulate")
    print("     2. Create New Hamiltonian")
    print("     3. Install Hamiltonian")
    print("     4. Uninstall Hamiltonian")
    print("     5. Clean installation")
    print("     6. Quit")


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
            sleep(2)
            break
    
    
    if not installed:
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
        print("Hamiltonian installed sucessfully")
        sleep(2)


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
    print("Hamiltonian uninstalled")
    sleep(2)



# Modifies simulation source file to set it to run the correct hamiltonian
# Can be found on line that reads 'simFunction = ......."
# For the Monte-Carlo simulation, also edits the partial hamiltonian line
# As an input, takes the name of a hamiltonian as it reads in hamiltonian.txt (assumes associated header file is there)


def setHam(name, monte):
    partialName = name + "Partial"
    name += "Ham" #Naming convention
    simFile = None
    if not monte:
        simFile = open("Simulate/BasicSim/Simulate.cpp" , 'r')
    else:
        simFile = open("Simulate/Monte-Carlo/MonteCarlo.cpp" , 'r')
    lookup = "simFunction ="
    lookup2 = "simPartial = "
    linenum = 0
    for num, line in enumerate(simFile, 0):
        if line.strip().startswith(lookup):
            linenum = num
            break

    simFile.seek(0)
    contents = simFile.readlines()
    simFile.close()
    contents[linenum] = '\t' + lookup + ' ' + name + ";\n"
    if monte:
        contents[linenum + 1] = '\t' + lookup2 + partialName + ";\n"
        oFile = open("Simulate/Monte-Carlo/MonteCarlo.cpp", "w")
        contents = "".join(contents)
        oFile.write(contents)
        oFile.close
    else:
        oFile = open("Simulate/BasicSim/Simulate.cpp", "w")
        contents = "".join(contents)
        oFile.write(contents)
        oFile.close

def newHam():
    print("Welcome to the Hamiltonian creation wizard:")
    name = str(input("Input a name for the new hamiltonian: "))
    if name.upper() == "TEMPLATE":
        print("ERROR: " + name + " is not a valid hamiltonian name")
        sleep(2)
        return
    filename = name + ".h"
    path = "hamiltonian/" + filename
    if os.path.isfile(path):
        print("Error: Hamiltonian with this name already exists")
        sleep(2)
        return

    writeCode = False
    if input("Would you like to write the calculation code now? (y/n) ").upper() == 'Y':
        writeCode = True

    if writeCode:
        writeCalcCode()
        calcFile = open(".tempcode", "r")
        calcContents = calcFile.readlines()
        calcFile.close()
        os.remove(".tempcode")
        lookup = "void calculate"
        while True:
            line = calcContents.pop(0)
            if line.startswith(lookup):
                break
        while True:
            line = calcContents.pop()
            if "//End Calculation" in line:
                break

    templateFile = open("hamiltonian/.template", 'r')
    lookup = "Template"
    tempContents = templateFile.readlines()
    templateFile.close()
    oFile = open(path, "w")
    for line in tempContents:
        line = line.replace("Template", name)
        line = line.replace("template", name)
        oFile.write(line)
        if line.startswith("void " + name + "::calculate"):
            if writeCode:
                oFile.write("\n")
                for item in calcContents:
                    oFile.write(item.strip('\n'))
                    oFile.write("\n")
            else:
                oFile.write("\n\tINPUT CALCULATION CODE HERE\n")


    oFile.close()
    print("File " + name + ".h created in the hamiltonian folder")
    if writeCode:
        if input("Would you like to install the hamiltonian? (y/n)").upper() == 'Y':
            installHam(filename)

    sleep(2)





def writeCalcCode():
    filename = ".tempcode"
    oFile = open(filename, "w")
    oFile.write("//The c++ code written in this file will be put into the calculation function of the hamiltinian\n")
    oFile.write("//The graph the hamiltonian is running in is accessed through the variable host\n")
    oFile.write("//The values of the complete calculation should be assigned to _result and _partial respectively\n")
    oFile.write("//The source term can be accessed via the sourceT variable.\n")
    oFile.write("//See the documentation for available functions\n\n")
    oFile.write("//IMPORTANT: Please SAVE and CLOSE the file when you are done editing it. You may add your own comments, but please do not edit the ones already present\n\n")
    oFile.write("void calculate(hGraph host) { \n\t\n\t\n} //End Calculation")
    oFile.close()
    proc = pexpect.spawn('open -a TextEdit .tempcode')
    proc.interact()
    input("When you are finished, save the file and press enter")



