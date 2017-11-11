#Can be used to buildGraphs and clean installation. Requires pexpect, which may not be in your base version of Python

import subprocess
import pexpect
import os
import re

def menu():
    print("Menu:")
    print("     0. Generate Graphs")
    print("     1. Clean installation")
    print("     2. Calculate Hamiltonians")
    print("     3. Install Hamiltonian")
    print("     4. Quit")


    access = int(input("Make a selection from the above list: "))
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

        proc = subprocess.Popen(['make', 'clean'],
                            stdout=subprocess.PIPE, cwd='buildGraphs')
        proc.communicate()[0].decode('utf-8')

    elif access == 2:
        break;

    elif access == 3:
        filename = input("Please enter the name of the header file: ")
        path = "hamiltonian/" + filename
        print(filename)

        try:

            open(path, "r")
            installHam(filename)

        except:
            print("Hamiltonian not found.")
            print("Make sure the spelling of the filename is correct and the file is located in the hamiltonian folder")



    elif access == 4:
        break
    
    else:
        print("Invalid option")


                



