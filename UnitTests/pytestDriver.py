import os
import shutil

#in Python, the listdir command  returns files and directories (like typeing in "dir").
listOfDirectoriesAndFiles = os.listdir(".")

#Below is going to become a list of directories only in the next loop.
directoryList = []
for elem in listOfDirectoriesAndFiles:
    if os.path.isdir(elem) == True:
        directoryList.append(elem)

#This loop goes into each directories, runs the specified command, and comes back.
for directory in directoryList:
    print("Changing directory to {}".format(directory))
    os.chdir(directory)
    try:
        shutil.rmtree("__pycache__")
    except:
        pass
    os.system("pytest")
    os.chdir("..")
