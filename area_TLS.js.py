###Intro

# Hi!

# Characteristics
	# written for Jython (implementation of Python 2.7 in ImageJ, version 2.9.0/1.53t)
	# contact for help: solenedeutsch@gmail.com

# Aim of the script
	# Measure TLS area, T- and B-cells area inside, return data as a csv file
	# Retrieve thresholding data as well for T- and B-cells channels

# /!\ Getting started /!\
	# opening the images: "Split channels" must be unchecked
	# in ImageJ, measurement results must include "Area" (Analyze -> Set Measurement...)
	# before starting, check the directory for the results file (the script will ask you where you want to save it)
	# running this script on one image is pretty fast, but it asks for user input for each image
		# depending on the size of your dataset, this could take some time
	# results are saved in a csv file called Results_area.csv (can be modified)
		# /!\ currently, every time you run this script, it will erase the content of the previous Results_area.csv and create a new one
		# if you want to save the previous results, change the name of the file before running the script again
		# or alternatively, don't forget to transfer the data to an Excel file before running the script again

# Have fun!

### Imports

from ij import IJ, WindowManager, Prefs
from ij.measure import ResultsTable
from ij.plugin.frame import RoiManager
from ij.plugin import ChannelSplitter, ZProjector
from ij.process import AutoThresholder
from ij.gui import WaitForUserDialog
import csv
import os

###Initializing

rm = RoiManager.getRoiManager()
DATA = []
field_names = ['Name', 'Area TLS', 'Area CD3', 'Threshold CD3', 'Area B220', 'Threshold B220']

count = int(IJ.getNumber("How many images to analyze ?", 1))

os.chdir("my_directory")

### Analysis

for _ in xrange(count):
	# Initializaing for the current image
	
	new_dict = {}
	
	whole = IJ.getImage()
	title = whole.getTitle()
	
	channels = ChannelSplitter.split(whole)
	
	Immune_channel = channels[2]
	T_channel = channels[3]
	B_channel = channels[0]
	
	# Borders of the TLS
	
	immune = ZProjector.run(Immune_channel,"sum")
	immune.show()
	IJ.setAutoThreshold(immune, "Otsu dark")
	Prefs.blackBackground = True
	IJ.run(immune, "Convert to Mask", "")
	IJ.run(immune, "Analyze Particles...", "size=800.00-Infinity circularity=0.00-0.50 show=Overlay clear add")
	ra = rm.getRoisAsArray()
	for ROI in ra:
		Immune_channel.setRoi(ROI)
	Immune_channel.show()
	check = WaitForUserDialog("Time for a manual check", "The ROIs in the ROI Manager are going to be fused together to create the TLS ROI. \nIf some of those ROIs are inapropriate, delete them. \nIf a part of the TLS was not detected, please draw the corresponding ROI manually (and sorry for that). \nClick Ok when finished.")
	check.show()
	Immune_channel.close()
	ra = rm.getRoisAsArray()
	if len(ra) >1:
		indexes = range(len(ra))
		rm.setSelectedIndexes(indexes)
		rm.runCommand(immune,"OR")
		roi = immune.getRoi()
		rm.addRoi(roi)
		rm.setSelectedIndexes(range(len(ra)))
		rm.runCommand(immune,"Delete")
		
	# Measurement for T-cells
	
	t_cells = ZProjector.run(T_channel,"max")
	IJ.setAutoThreshold(t_cells, "Shanbhag dark")
	ImProc = t_cells.getProcessor()
	lower_T = ImProc.getMinThreshold()
	upper_T = ImProc.getMaxThreshold()
	IJ.run(t_cells, "Convert to Mask", "")
	IJ.run(t_cells, "Create Selection", "")
	roi = t_cells.getRoi()
	rm.addRoi(roi)
	
	rm.setSelectedIndexes([0,1])
	rm.runCommand(t_cells,"AND")
	roi = t_cells.getRoi()
	rm.addRoi(roi)
	rm.setSelectedIndexes([1])
	rm.runCommand(t_cells,"Delete")
	
	# Measurement for B-cells
	
	b_cells = ZProjector.run(B_channel,"max")
	IJ.setAutoThreshold(b_cells, "Shanbhag dark")
	ImProc = b_cells.getProcessor()
	lower_B = ImProc.getMinThreshold()
	upper_B = ImProc.getMaxThreshold()
	IJ.run(b_cells, "Convert to Mask", "")
	IJ.run(b_cells, "Create Selection", "")
	roi = b_cells.getRoi()
	rm.addRoi(roi)
	
	rm.setSelectedIndexes([0,2])
	rm.runCommand(b_cells,"AND")
	roi = b_cells.getRoi()
	rm.addRoi(roi)
	rm.setSelectedIndexes([2])
	rm.runCommand(b_cells,"Delete")
	
	rm.setSelectedIndexes([0,1,2])
	rm.runCommand(immune,"Measure")
	
	immune.close()
	whole.close()
	
	# Save results
	
	table = ResultsTable.getResultsTable()
	areas = [area for area in table.getColumn(0)]
	
	new_dict['Name'] = title[37:]
	new_dict['Area TLS'] = areas[0]
	new_dict['Area CD3'] = areas[1]
	new_dict['Area B220'] = areas[2]
	new_dict['Threshold CD3'] = str(lower_T)+'/'+str(upper_T)
	new_dict['Threshold B220'] = str(lower_B)+'/'+str(upper_B)
	
	DATA.append(new_dict)

### Saving whole results in csv file

csvfile = open('Results_area.csv', 'w')
writer = csv.DictWriter(csvfile, fieldnames = field_names)
writer.writeheader()
writer.writerows(DATA)
csvfile.close()
	