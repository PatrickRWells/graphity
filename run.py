#Can be used to buildGraphs and clean installation. Requires pexpect, which may not be in your base version of Python

import subprocess
import pexpect
import os

def menu():
    print("Menu:")
    print("     0. Generate Graphs")
    print("     1. Clean installation")
    print("     2. Quit")

    access = int(input("Make a selection from the above list: "))
    return access

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
        break
    
    else:
        print("Invalid option")


                



