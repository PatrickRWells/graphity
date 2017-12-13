# 7.1 Write a program that prompts for a file name, then opens that file and reads 
# through the file, and print the contents of the file in upper case. Use the file 
# words.txt to produce the output below.
# You can download the sample data at http://www.pythonlearn.com/code/words.txt

# Use words.txt as the file name

fh = None
while fh is None:
    fname = raw_input("Enter file name: ")
    if len(fname) < 1 : fname = "graph10.txt"
    try:
        fh = open(fname)
        break
    except: 
        print 'File name is invalid. Please enter a correct file name.'
        continue


outFile = open("graph10out.txt","w");

# print len(line)
# print range(9)
# newline = ""
# for i in range(len(line)):
# 	print(i)	
# 	newline+=line[i]+" "
# 
# print newline
# print len(newline)

for line in fh:
	line = line.rstrip()
#  	print line	
#  	print len(line)
	newline = ""
	if len(line) == 10:
		for i in range(len(line)):
			newline+=line[i]+" "
		newline = newline.rstrip()
# 		print newline
# 		print len(newline)
	if len(newline) > 2:
		outFile.write(newline+"\n")
	
fh.close
outFile.close()