#Can be used to buildGraphs, simulate, and clean installation. Requires pexpect, which may not be in your base version of Python

import subprocess
import pexpect
import os
import re

def menu():
    print("Menu:")
    print("     0. Generate Graphs")
    print("     1. Simulate")
    print("     2. Install Hamiltonian")
    print("     3. Uninstall Hamiltonian")
    print("     4. Clean installation")
    print("     5. Quit")


    access = int(input("Make a numerical selection from the above list: "))
    return access


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

def setHam(name):
    name += "Ham"
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
        print("This functionality is not yet supported")

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
