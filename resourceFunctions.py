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


# Installs a hamiltonian for use in simulation
# Hamiltonains must be in the 'hamiltonian' subfolder
# Edits hamiltonian.txt with name and hamiltonians.h with an appropriate include statement
# Requires input of a filename (c header file) where hamiltonian is implemented

def installHam(name): #takes in the name of a file containing a new hamiltonian and installs it
    wFile = open("hamiltonian/hamiltonians.txt", "a+")
    wFile.seek(0)
    hamName = ""
    for char in name:
        if char == '.':
            break
        hamName += char #gets the hamiltonian name from the name of the header file
    
    installed = False
    
    for line in wFile:#Chekcs if a hamiltonian by that name is already installed
        if line.strip() == hamName:
            print("Hamiltonian has already been installed")
            wFile.close()
            installed = True
            sleep(2)
            break
    
    
    if not installed:
        wFile.write(hamName + '\n')#adds the hamiltonian to the end of the hamiltonians.txt file
        wFile.close()
        hFile = open("hamiltonian/hamiltonians.h", 'r')#opens the hamiltonians header file
        lookup = "#endif /* hamiltonians_h */"
        linenum = 0
        for num, line in enumerate(hFile, 0):
            if line.strip() == lookup: #finds the end of the file
                linenum = num
                break
        hFile.seek(0)
        contents = hFile.readlines()
        hFile.close()
        contents.insert(linenum, '#include "' + name + '"' + '\n')#puts the include statement for the new hamiltonian at the end of the file
        eFile = open("hamiltonian/hamiltonians.h", "w")#overwrites the old file
        contents = "".join(contents)
        eFile.write(contents)
        eFile.close
        print("Hamiltonian installed sucessfully")
        sleep(2)


def uninstallHam(name): #uninstalls a hamiltonian
    print(name)
    linenum = 0
    found = False
    rFile = open("hamiltonian/hamiltonians.txt", "r")
    for num, line in enumerate(rFile, 0): #finds the hamiltonian in the hamiltonians.txt file
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
    contents.pop(linenum) #removes the hamiltonian
    oFile = open("hamiltonian/hamiltonians.txt", "w")#overwrites the hamiltonians.txt file
    contents = "".join(contents)
    oFile.write(contents)
    oFile.close

    name += ".h"
    lookup = '#include "' + name + '"'
    hFile = open("hamiltonian/hamiltonians.h", "r")#Findes the hamiltonian header file include statement in hamiltonians.h
    for num, line in enumerate(hFile, 0):
        if line.strip() == lookup:
            linenum = num
            break
    hFile.seek(0)
    contents = hFile.readlines()
    hFile.close()
    contents.pop(linenum)#remobes the hamiltonian header file include function
    oFile = open("hamiltonian/hamiltonians.h", "w")#overwrites the hamiltonian header file
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
        simFile = open("Simulate/BasicSim/Simulate.cpp" , 'r')#Opens the code for the simulation
    else:
        simFile = open("Simulate/Monte-Carlo/MonteCarlo.cpp" , 'r')
    lookup = "simFunction ="
    lookup2 = "simPartial = "
    linenum = 0
    for num, line in enumerate(simFile, 0): #finds the lines that set the hamiltonian function to simulate on
        if line.strip().startswith(lookup):
            linenum = num
            break

    simFile.seek(0)
    contents = simFile.readlines()
    simFile.close()
    contents[linenum] = '\t' + lookup + ' ' + name + ";\n"
    if monte:
        contents[linenum + 1] = '\t' + lookup2 + partialName + ";\n" #Replaces with the new hamiltonian and overwrites the old file
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
    if name.upper() == "TEMPLATE":#Template is not a valid hamiltonian name
        print("ERROR: " + name + " is not a valid hamiltonian name")
        sleep(2)
        return
    filename = name + ".h"
    path = "hamiltonian/" + filename
    if os.path.isfile(path):
        print("Error: Hamiltonian with this name already exists")#Ensures hamiltonian doesn't already exist
        sleep(2)
        return

    writeCode = False
    if input("Would you like to write the calculation code now? (y/n) ").upper() == 'Y':
        writeCode = True

    if writeCode:
        writeCalcCode()#see appropriate file
        calcFile = open(".tempcode", "r")#code is saved as a dotfile so it does not appear in the finder
        calcContents = calcFile.readlines()
        calcFile.close()
        os.remove(".tempcode")
        lookup = "void calculate"
        while True: #The next two loops extrat the actual calculation code from the other pieces of the dotfile
            line = calcContents.pop(0)
            if line.startswith(lookup):
                break
        while True:
            line = calcContents.pop()
            if "//End Calculation" in line:
                break

    templateFile = open("hamiltonian/.template", 'r') #There is a dotfile named .template that contians a template hamiltonian
    lookup = "Template"
    tempContents = templateFile.readlines()
    templateFile.close()
    oFile = open(path, "w")
    for line in tempContents:
        line = line.replace("Template", name)#overwrites all instances of the word "Template" with the name of the new hamiltonian
        line = line.replace("template", name)
        oFile.write(line)
        if line.startswith("void " + name + "::calculate"): #inserts the calculation code where it is needed
            if writeCode:
                oFile.write("\n")
                for item in calcContents:
                    oFile.write(item.strip('\n'))
                    oFile.write("\n")
            else:
                oFile.write("\n\tINPUT CALCULATION CODE HERE\n")#if the user did not write the code right away, inserts a clear message (also prevents compilation)


    oFile.close()
    print("File " + name + ".h created in the hamiltonian folder")
    if input("Would you like to install the hamiltonian? (y/n)").upper() == 'Y':
        installHam(filename)

    sleep(2)





def writeCalcCode():
    filename = ".tempcode"
    oFile = open(filename, "w")#creates a filename calld temp code
    calcCodeHeader(oFile)      #Writes some information into the header
    oFile.write("void calculate(hGraph host) { \n\t\n\t\n} //End Calculation")
    oFile.close()
    proc = pexpect.spawn('open -a TextEdit .tempcode') #opens the file in textedit
    proc.interact() #Waits until the user closes textedit
    input("When you are finished, save the file and press enter")

def editCalcCode(name): #Allows the user to edit the calculation code of an already-installed hamiltonian
    filename = name + ".h"
    path = "hamiltonian/" + filename
    iFile = open(path, "r")
    calcContents = iFile.readlines()
    iFile.close()
    while True:#The next several lines extract the current calculation code from the header file
        line = calcContents.pop(0)
        if line.strip().startswith("void " + name + "::calculate"):
            break

    calcContents.reverse()


    found = False
    while True:
        line = calcContents.pop(0)
        if line.strip().startswith("void " + name + "Ham"):
            found = True
        if found and "}" in line:
            break

    calcContents.reverse()
    calcFile = open(".tempcode", "w")#Creates a new dotfile with the current calculation code
    calcCodeHeader(calcFile)
    calcFile.write("void calculate(hGraph host) {\n")
    for line in calcContents:
        calcFile.write(line)
    calcFile.write("}")
    calcFile.close()

    proc = pexpect.spawn('open -a TextEdit .tempcode')#Opens the dotfile in textedit to be edited by the user
    proc.interact()
    input("When you are finished, save the file and press enter")
    newCalcFile = open(".tempcode", "r")
    data = newCalcFile.readlines()
    newCalcFile.close()
    os.remove(".tempcode")

    while True: #The next to loops grab the calculation code from the dotfile
        line = data.pop(0)
        if line.strip().startswith("void calculate"):
            break
    while True:
        line = data.pop()
        if line.strip().startswith("}"):
            break


    templateFile = open("hamiltonian/.template", "r") #Creates a new file based on the template with the hamiltonian name and the new calculation code
    found = False
    i = 0
    for line in templateFile:
        line = line.replace("Template", name)
        line = line.replace("template", name)
        if not found:
            data.insert(i, line)
            i += 1
        else:
            data.append(line)

        if line.strip().startswith("void " + name + "::calculate"):
            found = True

    finalOut = open(path, "w")#overwrites the old header file with the new data
    for line in data:
        finalOut.write(line)
    print("Hamiltonian edited sucessfully")
    sleep(2)



def calcCodeHeader(oFile):
    oFile.write("//The c++ code written in this file will be put into the calculation function of the hamiltonian\n")
    oFile.write("//The graph the hamiltonian is running in is accessed through the variable host\n")
    oFile.write("//The values of the complete calculation should be assigned to _result and _partial respectively\n")
    oFile.write("//The source term can be accessed via the sourceT variable.\n")
    oFile.write("//See the documentation for available functions\n\n")
    oFile.write("//IMPORTANT: Please SAVE and CLOSE the file when you are done editing it. You may add your own comments, but please do not edit the ones already present\n\n")



