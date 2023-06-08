###Intro

# Hi!

# Characteristics
	# written for Jython (implementation of Python 2.7 in ImageJ, version 2.9.0/1.53t)
	# contact for help: solenedeutsch@gmail.com

# Aim of the script
	# Measure the size of a tumor

# /!\ Getting started /!\
	# opening the images: "Split channels" must be unchecked
	# in ImageJ, measurement results must include "Area" (Analyze -> Set Measurement...)
	# before starting, check the directory for the results file
	# results are saved in a text file called results_area_tumor.txt (can be modified)
		# /!\ currently, every time you run this script, it will erase the content of the previous results file and create a new one
		# if you want to save the previous results, change the name of the file before running the script again
	# check which channel you want to use (start counting from 0 instead of 1) and change that in the initializing for the current image part

# Have fun!

### Imports

from ij import IJ, WindowManager, Prefs
from ij.measure import ResultsTable
from ij.plugin.frame import RoiManager
from ij.plugin import ChannelSplitter
from ij.process import AutoThresholder
from ij.gui import WaitForUserDialog
import os

###Initializing

rm = RoiManager.getRoiManager()

DATA = {}

count = int(IJ.getNumber("How many images to analyze ?", 1))

os.chdir("my_directory")

### Analysis

for _ in xrange(count):
	# Initializing for the current image

	whole = IJ.getImage()
	title = whole.getTitle()

	channels = ChannelSplitter.split(whole)

	tumor_channel = channels[0]

	# Borders of the tumor

	tumor = tumor_channel.duplicate()
	tumor.show()
	IJ.setAutoThreshold(tumor, "Otsu dark")
	Prefs.blackBackground = True
	IJ.run(immune, "Convert to Mask", "")
	check = WaitForUserDialog("Time for a manual check", "This is for you to check the threshold, you can delete that part once you're confident in the thresholding method.")
	check.show()
	IJ.run(immune, "Analyze Particles...", "size=300-Infinity circularity=0.00-1.00 show=Overlay clear add")
	check = WaitForUserDialog("Time for a manual check", "Check the areas that the script detected (useful for debugging).")
	check.show()
	#Calculating area
	ra = rm.getRoisAsArray()
	if len(ra) >1:
		#Combining teh ROIs if several have been detected
		indexes = range(len(ra))
		rm.setSelectedIndexes(indexes)
		rm.runCommand(tumor,"OR")
		roi = immune.getRoi()
		rm.addRoi(roi)
		rm.setSelectedIndexes(range(len(ra)))
		rm.runCommand(tumor,"Delete")
		rm.setSelectedIndexes([0])
		rm.runCommand(tumor,"Measure")
		table = ResultsTable.getResultsTable()
		area = table.getColumn(0)[0]
		DATA[title] = area
	elif len(ra) == 1:
		rm.setSelectedIndexes([0])
		rm.runCommand(tumor,"Measure")
		table = ResultsTable.getResultsTable()
		area = table.getColumn(0)[0]
		DATA[title] = area
	else:
		DATA[title] = 0

	tumor.close()
	whole.close()

### Saving whole results in csv file

results = open('Results_area_tumor.txt', 'w')
results.write(str(DATA))
results.close()