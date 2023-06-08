### Intro

# Hi!

# Characteristics:
	# script written for Jython (Python 2.7 implemented in ImageJ, current version 1.53t)
	# contact for help: solenedeutsch@gmail.com

# Aim of the script
	# Running a colocalization analysis between two channels and retrieving the result as a text file
		# the data can then be analyzed using the colocalization_stat.py python script

# /!\ Getting started /!\
	# opening the images : "Split channels" box must be checked
	# this script is currently only suited for z-stacks
	# this script can currently handle only one marker to analyze in the colocalization with the vessels
	# if there are a lot of images to be analyzed, run the script on the HIVE to avoid memory issues
	# in ImageJ, measurement has to include "Stack position" and "Min Max intensity" (Analyze -> Set Measurement...)
	# a first try to ensure the thresholding works correctly is necessary for every new dataset
	# before starting : check the naming part in run_coloc() and modify according to your data
	# the script will take a bit of time to run. Go have Fika in the meantime.
	# if you are working on your computer while the script is running be aware that:
		# it might slow down your computer (if it's really annoying, use the HIVE)
		# some images might pop up at the beginning: do not close them
	# change the directory in which you want to save the result files in line 48
	# the results are written in a file, you can change the name of the file in the indicated line at the end of the script
		# /!\ if you keep a file name that is already taken by another one of your files, this other file will be erased

# Code details for curious people
	# xrange() is used instead of range() for memory optimization (Python 2.7) when possible
	# reversed() is used every time there is a deletion based on indexes, so as to not disturb the loop

# Have fun !

### Imports

from ij import IJ, WindowManager, Prefs
from ij.plugin import ImageCalculator, ZProjector
from ij.measure import ResultsTable
from ij.plugin.frame import RoiManager
from ij.gui import WaitForUserDialog
import os

### Global variables

DATA = {}

os.chdir("my_directory")

### Part 1: colocalization analysis

def select_vessel(OriImage):
	"""
	Takes the z-stack containing the vessel staining and returns a z projection to isolate the vessels
	list of relevant variables :
		OriImage = original image containing the vessel staining
		Thresholded = same image to be thresholded
		Filtered = thresholded image on which a median filter is applied
		FinalIg = image displaying only the vessels
	"""
	# Delimiting borders around the vessels
	Thresholded = OriImage.duplicate()
	IJ.setAutoThreshold(Thresholded, "Moments dark")
	IJ.run(Thresholded, "Make Binary", "method=Moments background=Dark calculate")
	Filtered = Thresholded.duplicate()
	IJ.run(Filtered, "Median...", "radius=2")
	IJ.run(Filtered, "Make Binary", "method=Otsu background=Dark calculate")
	ic = ImageCalculator()
	FinalIg = ic.run("Subtract create stack", Thresholded, Filtered)
	# 'Cleaning up'
	for _ in range(7):
		IJ.run(FinalIg, "Despeckle", "stack")
	IJ.run(FinalIg, "Options...", "iterations=2 count=1 do=Dilate stack")
	IJ.run(FinalIg, "Options...", "iterations=10 count=1 do=Close stack")
	# z-projection
	FinalIg = ZProjector.run(FinalIg,"max")
	return FinalIg

def define_ROI(vessels, OriImage, control = False):
	'''
	Takes the z-projection and returns a list of z-stacks, each containing a single vessel.
	Both parameters are images, returns a list of images.
	'''
	# Isolate individual vessels
	rm = RoiManager.getRoiManager()
	IJ.run(vessels, "Analyze Particles...", "size=30-Infinity show=Nothing display clear add stack")
	table = ResultsTable.getResultsTable()
	n_roi = rm.getCount() # number of ROIs (=number of vessels detected)
	n_slices = int(table.getValue("Slice", n_roi-1)) # number of slices per z-stack
	list_vessels = []
	ra = rm.getRoisAsArray()
	# Apply ROI and clear outside to keep only the vessel part
	for i in range(n_roi):
		single_vessel = OriImage.duplicate()
		single_vessel.setRoi(ra[i])
		IJ.setBackgroundColor(0,0,0)
		IJ.run(single_vessel, "Clear Outside", "stack")
		IJ.run(single_vessel, "Select None", "stack")
		list_vessels.append(single_vessel)
		if control:
			single_vessel.show()
			myWait = WaitForUserDialog("Time to look at the data", "Check if the detected shape corresponds to an actual vessel.\nIf it is not the case, you might want to change the thresholding method.")
			myWait.show()
			single_vessel.changes = False
			single_vessel.close()
	return list_vessels

def blur_image(current_image):
	'''
	Adds a Gaussian filter to try and remove a bit of noise.
	Parameter is an image, returns an image.
	'''
	duplicate_image = current_image.duplicate()
	IJ.run(duplicate_image, "Gaussian Blur...", "sigma=2 stack")
	return duplicate_image

def coloc(current_vessel, current_chan1):
	'''
	Runs the colocalization analysis on the the two channels current_vessel and current_chan1.
	Results of the analysis are displayed in both a colocalization results window, and the log window.
	List of relevant variables :
		OriImage = original image containing the vessel staining
		Thresholded = same image to be thresholded using the triangle method
		Filtered = thresholded image on which a median filter is applied
		FinalIg = image displaying only the vessels
	'''
	# Here we refine the ROI by making a mask -> more precise ROI for each slice of the z-stack
	# Delimiting borders around the vessel
	Thresholded = current_vessel.duplicate()
	IJ.setAutoThreshold(Thresholded, "Triangle dark")
	IJ.run(Thresholded, "Make Binary", "method=Otsu background=Dark calculate")
	Filtered = Thresholded.duplicate()
	IJ.run(Filtered, "Median...", "radius=2 stack")
	IJ.run(Filtered, "Make Binary", "method=Otsu background=Dark calculate")
	ic = ImageCalculator()
	FinalIg = ic.run("Subtract create stack", Thresholded, Filtered)
	IJ.run(FinalIg, "Options...", "iterations=5 count=1 do=Close stack")
	for _ in range(5):
		IJ.run(FinalIg, "Despeckle", "stack")
	IJ.run(FinalIg, "Options...", "iterations=5 count=1 do=Dilate stack")
	Prefs.blackBackground = True
	IJ.run(FinalIg, "Convert to Mask", "method=Otsu background=Light calculate")
	IJ.run(FinalIg, "Invert LUT", "")
	FinalIg.show()
	ROI = FinalIg.getTitle()
	# Applying Gaussian filters to remove some noise
	current_vessel = blur_image(current_vessel)
	current_vessel.show()
	title_vessel = current_vessel.getTitle()
	current_chan1.show()
	title_chan1 = current_chan1.getTitle()
	# Running the colocalization analysis
	IJ.run("Coloc 2", "channel_1=["+title_chan1+"] channel_2=["+title_vessel+"] roi_or_mask=["+ROI+"] threshold_regression=Costes manders'_correlation psf=3 costes_randomisations=10")
	# Closing useless windows
	FinalIg.changes = False
	FinalIg.close()

def retrieve_data(isolate_vessels, image_name, vessel_number):
	'''
	Stores Manders tM1 and tM2 values (from the colocalization analysis of image_name) in the data dictionary, then adds the mean pixel intensity as well
	key -> image_name (short for readability)
	value -> list containing tM1 value and tM2 value
	'''
	# Creates an entry in the dictionary
	DATA[image_name] = []
	# Retrieves text from the log window
	logText = IJ.getLog()
	list_log = logText.split("\n")
	# Takes specifically Manders' M1 and M2 values
	for i in range(len(list_log)):
		if "Manders" in list_log[i] and "tM" in list_log[i]:
			for k in range(len(list_log[i])):
				if list_log[i][k] == ",":
					value = list_log[i][k+2:]
					# Adds values to the correct dictionary entry
					DATA[image_name].append(value)
	# Closes the log window to clear it for the next image
	IJ.selectWindow("Log")
	IJ.run("Close")

def run_colocalization():
	'''
	Runs the entire colocalization process for a batch of images.
	Colocalization data is stored in the DATA dictionnary and in a text file.
	Additionally: returns a list containing the names of the individual vessels, and another containing all the images names (format = that of the t-cell channel).
	'''
	# Number of time to run the colocalization (= number of images to analyze)
	count = int(IJ.getNumber("How many images to analyze ?", 1))
	control = bool(IJ.getString("Do you want to control the detection of vessel shapes? (True or False)", "False"))
	format = IJ.getString("Will you use excel or python for downstream analysis (to change the output format)?", "python")

	# Parameters of the experiment
	experiment_name = IJ.getString("Enter name of lif file", "")
	stainings = IJ.getString("Enter stained molecules, in the order and case of the image titles", "")
	panels = IJ.getString("Enter panel description for the full experiment (ex: P1-2 P2-2 for 2 panels with two images in each).\nIf there is only one image to analyze, the name can be like P1-5 for panel 1, image number 5.", "")
	magnification = IJ.getString("Enter magnification", "x63")
	vessel_channel = IJ.getString("Enter number of channel for vessel staining", "2")
	coloc_channel = IJ.getString("Enter number of channel for the other staining to be used in the colocalization analysis", "1")
	variable_part = "" # to be modified depending on the experiment
	panel_description = panels.split(" ")
	# Images to be analyzed (whole batch or only one image for testing)
	if count == 1:
		image_list = [panels]
	else:
		image_list = []
		# Getting the short version of the images names for the dictionary's keys
		for images in panel_description:
			for image_number in range(int(images[4:])):
				image_list.append(images[0:3]+"_"+str(int(image_number)+1))
	# Adding a list containing every vessel code name
	total_vessels = []

	# Colocalization of each images
	for i in range(count):
		# Construction of image titles
		current_panel = image_list[i][0:3]
		current_panel_image = image_list[i][4:]
		current_dict_key = image_list[i]
		current_image_name = experiment_name+".lif - "+variable_part+current_panel+" "+stainings+" "+magnification+" "+current_panel_image
		# Channels to be analyzed
		curr_vessel = current_image_name+" - C="+vessel_channel
		current_chan1 = current_image_name+" - C="+coloc_channel
		# Making a blurred version of the first channel
		IJ.selectWindow(current_chan1)
		active_chan1 = IJ.getImage()
		blur_current_chan1 = blur_image(active_chan1)
		# Selecting the vessels in the image
		IJ.selectWindow(curr_vessel)
		curr_vessel = IJ.getImage()
		isolate_vessels = select_vessel(curr_vessel)
		current_vessel_list = define_ROI(isolate_vessels, curr_vessel, control)
		# Apply colocalization analysis to each vessels
		for k in xrange(len(current_vessel_list)):
			current_dict_key_vessel = current_dict_key + "_"+ str(k) # image code name
			current_vessel_list[k].setTitle(current_dict_key_vessel)
			total_vessels.append(current_vessel_list[k].getTitle())
			coloc(current_vessel_list[k], blur_current_chan1) # performing colocalization analysis
			retrieve_data(isolate_vessels, current_dict_key_vessel, k) # storing values in the dictionary
		# Closing useless windows
		blur_current_chan1.changes = False
		blur_current_chan1.close()
		# Sort of a timer
		print(str(count-i-1)+" images remaining")
	# Printing data to check the results
	print(DATA)
	print("finished colocalization")
	#Saving the data in the proper format
	if format == "python":
		file = open("result_file_colocalization.txt", "w") #change the name here
		file.write(str(DATA))
		file.close()
	else:
		field_names = ['Name', 'tM1', 'tM2']
		excel_data = []
		for key in DATA.keys():
			new_dict = {}
			new_dict['Name'] = key
			new_dict['tM1'] = DATA[key][0]
			new_dict['tM2'] = DATA[key][1]
			excel_data.append(new_dict)
		csvfile = open('excel_results_colocalization.csv', 'w') #change the name here
		writer = csv.DictWriter(csvfile, fieldnames = field_names)
		writer.writeheader()
		writer.writerows(excel_data)
		csvfile.close()
	return total_vessels

run_colocalization()