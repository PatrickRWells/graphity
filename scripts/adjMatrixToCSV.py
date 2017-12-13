#converts text file with space seperated adjacency matrix into csv file that can be read into simulation


while True:
    filename = input("Please enter the name of the graph file: ")
    try:
        open(filename, "r")
        break
    except:
        print("File not found.")

size = input("Please enter number of nodes in the graphs: ")
intSize = int(size)

outFile = input("Please ener an ouptut file name: ")
iFile = open(filename, "r")
oFile = open(outFile, 'w')
oFile.write(size + '\n')

contents = ""

for num, line in enumerate(iFile,0):
    contents = contents + line.strip('\n') + ' '
    if ( (num + 1) % intSize) == 0:
        contents = contents.replace(' ', ',')
        contents = contents + "\n"
        oFile.write(contents)
        contents = ""

iFile.close()
oFile.close()

