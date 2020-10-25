# Code that creates color-magnitude diagrams using data from Modules for Experiments in Stellar Astrophysics

import numpy as np
import matplotlib.pyplot as plt
import statistics
import os

# filename of color magnitude diagram
CMD_name = 'CMD.png'

# gather all directories
cwd = os.getcwd()
parent_dir = cwd+"/test"
dirs = os.listdir(parent_dir)

abs_V = []
abs_I = []

# color index V-I
abs_V_minus_I = []

# count used files
used_files = 0

for dir in dirs:
	if os.path.exists(parent_dir+"/"+dir+"/history.data"):
		current_file = parent_dir+"/"+dir+"/history.data"
		#print(current_file + " exists.")
		used_files += 1
	else:
		#print(current_file + " does NOT exist. \n")
		continue

	data = np.loadtxt(current_file, skiprows=6)
	#print("Loaded " + current_file) 

	temp_V = np.log10(data[:, 32])
	temp_I = np.log10(data[:, 33])
	temp_V_minus_I = temp_V - temp_I
	
	abs_V.append(statistics.fmean(temp_V))
	abs_I.append(statistics.fmean(temp_I))
	
	abs_V_minus_I.append(statistics.fmean(temp_V_minus_I))
	#print("Added absolute magnitude data to array for plotting")
	#print("Completed " + current_file + "\n")


plt.plot(abs_V_minus_I, abs_V, 'k.', markersize=3)
plt.xlabel(r'$M_{V} - M_{I}$ (mag)')
plt.ylabel(r'$M_{V}$ (mag)')

plt.gca().invert_yaxis()

#plt.show()

plt.savefig(CMD_name)
print("Color magnitude diagram saved as " + CMD_name + " in " + cwd)
print("Total directories: " + str(len(dirs)))
print("Total used files: " + str(used_files))

plt.clf()


