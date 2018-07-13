import os

#in Python, the listdir command  returns files and directories (like typeing in "dir").
listOfDirectoriesAndFiles = os.listdir(".")


#Below is going to become a list of directories only in the next loop.
directoryList = []
for elem in listOfDirectoriesAndFiles:
    if os.path.isdir(elem) == True:
        directoryList.append(elem)

#This loop goes into each directories, runs the specified command, and comes back.
for directory in directoryList:
    print("\nChanging directory to"+directory)
    os.chdir(directory)
    listOfFilesInDirectory=os.listdir(".")\
    
    #Loops through each of the files in the directory and runs any file that begins with 'test_'
    for name in listOfFilesInDirectory:
        if "test_" in name:
            print('\n'+ name)
            os.system("python " + name)
            
    os.chdir("..")