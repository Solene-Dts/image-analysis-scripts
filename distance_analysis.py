### Intro

# Hi!

# Characteristics:
    # date: 16/03/23
    # /!\ Python version 3.8, matplotlib version 3.7 /!\
    # contact for help: solenedeutsch@gmail.com

# Aim of the script
    # Plotting relevant graphs for T-cell distance analysis

# Getting started
    # Check that you're working in the directory where your text files are (os.getcwd() in the shell)
    # If the directory is wrong use os.chdir("my_right_directory") in the shell
    # This scripts calculate statistical significance but IT DOES NOT PLOT IT AUTOMATICALLY
        # According to the results of the satistical tests, change what you want to show on the plot in the indicated lines
    # Choose which graph you want to plot 
        # Write in the command shell: name_of_the_plotting_function()

# Have fun!

### Imports
import matplotlib.pyplot as plt #plotting graphs
import numpy as np #useful maths
import random as rd #generating random numbers
import os #handle working directories
import ast #data formatting
import tkinter as tk #interactive windows
import tkinter.filedialog as fd #retrieving files
import scikit_posthocs as sp #stats
import scipy.stats as st #more general stats

### Global variables
# Please feel free to choose better colors...
PALETTE = [np.array([(3/255, 3/255, 17/255)]), np.array([(62/255, 15/255, 114/255)]), np.array([(122/255, 34/255, 129/255)]), np.array([(167/255, 49/255, 125/255)]), np.array([(215/255, 69/255, 107/255)]), np.array([(242/255, 99/255, 92/255)]), np.array([(254/255, 168/255, 115/255)])]

###Preliminary functions

def select_files():
    '''
    Asks user for input. Return a list with correct file paths.
    '''
    # Retrieving data from text file
    merging = int(input('How many tumors do you want to combine ? (0 or more than 1) ')) #0 (no combination), or any other number (combination)
    root = tk.Tk()
    #Merging files if necessary
    if bool(merging):
        #Making user open files
        files = list(fd.askopenfilenames(title="Choose result files to combine"))
        root.destroy()
    else:
        #Making user open files
        files = [fd.askopenfilename(title="Choose result file")]
        #Avoiding tk window freeze
        root.withdraw()
    return files

def collect_data_from_file(file_title):
    '''
    Opens a text file and retrieves the data as a dictionary.
    '''
    file = open(file_title, 'r')
    content = file.readlines()
    file.close()
    data = ast.literal_eval(content[0])
    return data

def merge_dict_tot(ref, list_dict):
    '''
    Combine several dictionnaries into ref.
    list_dict: list of dictionaries to merge.
    '''
    for element in list_dict:
        ref = {**ref, **element}
    return ref

def split_dict_by_info(ref):
    '''
    Takes a dictionary and split it based on the kind of information it contains.
    Returns two dictionaries: one with vessel name and size, one with the distance data.
    '''
    #Two new dictionaries
    vessel_info = {}
    distance_info = {}
    #Splitting
    for vessel in ref.keys():
        vessel_info[vessel] = ref[vessel][:2] #Taking the first two elements of the list
        distance_info[vessel] = ref[vessel][2:] #Taking the remaining elements
    return vessel_info, distance_info

def merge_dict_info(vessel_info, distance_info):
    '''
    Merge two dictionaries based on the information they contain
    '''
    #New dictionary
    whole = {}
    #Merging
    for vessel in vessel_info.keys():
        whole[vessel] = vessel_info[vessel]
        if vessel in distance_info.keys():
            #Taking care of single numbers vs list of numbers format
            if type(distance_info[vessel]) != list:
                whole[vessel].append(distance_info[vessel])
            else:
                whole[vessel] += distance_info[vessel]
    return whole

def remove_doubles(dico, position):
    '''
    Loops through every vessel of every image in the dico dictionnary.
    key : vessel code name
    value : distance of the T-cell to that vessel
    If a T-cell appears more than once in an image, keeps only the vessel that is closest.
    Like this, every T-cell is only counted once in the analysis.
    Returns a dictionary, same key and values, with only the relevant information.
    '''
    #Creating new dictionary
    new = {}
    #Control: how many T cell are we counting before correction?
    count = 0
    for value in list(dico.values()):
            count+=len(value)
    print("Before correction : "+ str(count))
    #Looping on every vessel
    images_done = []
    for vessel in dico.keys():
        #Getting the common part of the name of the vessels from the same image
        current_image = vessel[:position]
        #Analyze grouped for every image (which contains the same detected T cells)
        if current_image not in images_done:
            images_done.append(current_image)
            n_tcells = len(dico[vessel])
            #Taking into account every vessel from that image
            to_analyze = [key for key in dico.keys() if current_image in key]
            #Applying correction to every T cell
            for t_cell in reversed(range(n_tcells)):
                previous = ""
                positions = []
                too_far = []
                for vessel in to_analyze:
                    #First time encountering a close vessel
                    if dico[vessel][t_cell] != "too far" and previous == "":
                        previous = dico[vessel][t_cell]
                        positions.append((vessel, dico[vessel][t_cell]))
                    #Encountering other close vessels
                    elif dico[vessel][t_cell] != "too far" and previous != "":
                        positions.append((vessel, dico[vessel][t_cell]))
                    #Encountering a vessel far away
                    else:
                        too_far.append(vessel) #Name of the vessel
                #If there are only "too far" vessels for that T cell
                # Randomly assigning that T cell to only one vessel
                if previous == "":
                    # Random assignment
                    choosen = rd.choice(too_far)
                    positions = [(choosen, "too far")]
                # If there are several vessels close to that T cell
                # Assigning that T cell to the closest vessel
                if len(positions) > 1:
                    min_dist = positions[0][1]
                    min_vessel = positions[0][0]
                    for i in range(1, len(positions)):
                        if positions[i][1] < min_dist:
                            min_dist = positions[i][1]
                            min_vessel = positions[i][0]
                    final_vessel = min_vessel
                    final_value = min_dist
                else:
                    final_vessel = positions[0][0]
                    final_value = positions[0][1]
                #New dict with the corrected positions
                if final_vessel in new.keys():
                    #Case where the key already exists
                    new[final_vessel].append(final_value)
                else:
                    #Case where the key doesn't already exist
                    new[final_vessel] = [final_value]
    #Control: how many T cells are we counting after correction?
    count = 0
    for value in list(new.values()):
            count+=len(value)
    print("After correction: "+str(count))
    return new

def prep_data(processed_dict, average = False):
    '''
    Makes a suitable dictionnary to build a histogram out of.
    key -> tags on x axis
    values -> height of the bar for each category
    Computes the average distance of T cells from each vessel.
    '''
    #Creating new dictionary
    histo = {}
    #Converting distances to proper tags for the histogram
    distances = [0, 1, 2, 3, 4, 5, "too far"]
    #Counting the number of T cells per distance category
    for vessel in processed_dict.keys():
        histo[vessel] = []
        for i in range(len(distances)):
            count = 0
            for value in processed_dict[vessel]:
                if value == distances[i]:
                    count+=1
            histo[vessel].append(count)
    if average:
        for vessel in histo.keys():
            histo[vessel] = np.mean(histo[vessel])
    return histo

def better_prep_data(processed_dict):
    '''
    Makes a suitable dictionnary to make a histogram out of.
    key -> vessel name
    value -> average distance of the T cells to that vessel
    Computes the average T cell distance from a vessel.
    '''
    #Creating new dictionary
    new = {}
    for vessel in processed_dict.keys():
        #Computing the average distance from each vessel
        cell_distances = processed_dict[vessel]
        for t in range(len(cell_distances)):
            if cell_distances[t] == "too far":
                cell_distances[t] = 6
        #Handling null data
        if len(processed_dict[vessel]) > 0:
            #Computing average distance
            new[vessel] = np.mean(cell_distances)
    return new
    
def lists_are_easier(histo_dict):
    '''
    Makes a list of values out of a dictionary suited for histogram plotting.
    '''
    results = []
    for i in range(9): #area + intensity + 7 distance categories
        results.append([])
        for vessel in histo_dict.keys():
            if len(histo_dict[vessel]) > i:
                results[i].append(histo_dict[vessel][i])
            else :
                results[i].append(0)
    return results

def merge_lists(ref, new):
    '''
    Handling null data.
    '''
    for i in range(len(new)):
        if new[i] != "nan":
            ref[i].append(new[i])
    return ref

def prep_histo(list_data):
    '''
    Makes a suitable dictionnary to build a histogram out of.
    key -> tags on x axis
    values -> height of the bar for each category
    Computes the average number of T cells per distance category.
    '''
    #Creating new dictionary
    histo = []
    #Converting distances to proper tags for the histogram
    distances = [0, 1, 2, 3, 4, 5, "too far"]
    #Total number of counted T cells for normalization
    n_tcells = len(list_data)
    #Counting the number of T cell per distance category
    for i in range(len(distances)):
        count = 0
        for value in list_data:
            if value == distances[i]:
                count+=1
        #Normalizing over total number of T cells
        histo.append(count/n_tcells)
    return histo

def barplot_annotate_brackets(num1, num2, data, center, height, yerr=None, dh=.05, barh=.05, fs=None, maxasterix=None):
    """ 
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05

        while data < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                break

        if len(text) == 0:
            text = 'n. s.'

    lx, ly = center[num1], height[num1]
    rx, ry = center[num2], height[num2]

    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh)

    plt.plot(barx, bary, c='black')

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs)

### Plotting graphs

def histo_average_distance():
    '''
    Plots a histogram showing the average T cell distance from each vessel category, based on the intensity of a marker's staining in those vessels.
    '''
    #Collecting files
    files = select_files()
    whole_data_list = [[] for _ in range(5)]
    #Data processing
    for file in files:
        #Retrieve proper data
        raw_data = collect_data_from_file(file)
        vessel_info, vessel_distance = split_dict_by_info(raw_data)
        #Correct T cell count
        position = [pos for pos, char in enumerate(list(vessel_info.keys())[0]) if char == "_"][-1]
        processed_data = remove_doubles(vessel_distance, position)
        #Get the average distance to each vessel
        graph_data = better_prep_data(processed_data)
        #Get proper data format for plotting
        final_data = merge_dict_info(vessel_info, graph_data)
        #Split into the different vessel categories
        histo_data = [[] for i in range(5)]
        for vessel in final_data.keys():
            if final_data[vessel][1] <= 25 and len(final_data[vessel]) >2:
                histo_data[0].append(final_data[vessel][-1])
            elif 25 < final_data[vessel][1] <= 50 and len(final_data[vessel]) >2:
                histo_data[1].append(final_data[vessel][-1])
            elif 50 < final_data[vessel][1] <= 75 and len(final_data[vessel]) >2:
                histo_data[2].append(final_data[vessel][-1])
            elif 75 < final_data[vessel][1] <= 100 and len(final_data[vessel]) >2:
                histo_data[3].append(final_data[vessel][-1])
            elif 100 < final_data[vessel][1] and len(final_data[vessel]) >2:
                histo_data[4].append(final_data[vessel][-1])
        #Takes the average distance for each vessel category
        for i in range(5):
            if len(histo_data[i])>0:
                histo_data[i] = np.mean(histo_data[i])
            else:
                #Handling null data
                histo_data[i] = "nan"
        whole_data_list = merge_lists(whole_data_list, histo_data)
    #Get proper distance value
    for element in whole_data_list:
        for i in range(len(element)):
            element[i] = element[i]*10
    
    #Calculate statistics
    for i in range(5):
        #Normal distribution test
        print(st.kstest(whole_data_list[i], "norm"))
    #Homogeneity of variances test
    print(st.bartlett(whole_data_list[0], whole_data_list[1], whole_data_list[2], whole_data_list[3], whole_data_list[4]))
    #No normal distribution: one-way analysis of variances
    print(st.kruskal(whole_data_list[0], whole_data_list[1], whole_data_list[2], whole_data_list[3], whole_data_list[4]))
    #If precedent test is significant: post-hoc analysis
    print(sp.posthoc_dunn(whole_data_list, p_adjust = "bonferroni"))
    
    #Calculate error bars to plot
    histo_std = [np.std(list_values) for list_values in whole_data_list]
    
    #Gets the height of the bars for each category
    histo_height = []
    for i in range(5):
        histo_height.append(np.mean(whole_data_list[i]))
    
    #Plot parameters
    bars = ["[0, 25]", "]25, 50]", "]50, 75]", "]75, 100]", "]100, 255]"]
    bars_num = [1,2,3,4,5]
    
    #Plotting the graph
    fig, ax = plt.subplots(1, 1, figsize = (15,15))
    
    #Setting font size for axis labelling
    plt.rc("xtick", labelsize = 20)
    plt.rc("ytick", labelsize = 20)
    
    #Plotting the histogram
    ax.bar(bars_num, histo_height, color = PALETTE[3])
    
    #Adding error bars
    ax.errorbar(bars_num, histo_height, yerr = histo_std, fmt = "none", color =PALETTE[1])
    
    #Adding Significance indications TO BE MODIFIED ACCORDING TO STATISTICS RESULTS
    #barplot_annotate_brackets(0, 4, 1, bars_num, histo_height, dh=0.2, fs=20)
    #barplot_annotate_brackets(0, 2, 0.22, bars_num, histo_height, dh=0.1, fs=20)
    
    #Setting axis labels
    ax.set_xticks(bars_num)
    ax.set_xticklabels(bars)
    #Hiding right and top spines
    ax.spines[["right", "top"]].set_visible(False)
    
    #Setting axis titles
    ax.set_ylabel("Average T-cell distance to the vessel (µm)", labelpad = 10, fontsize = 20)
    ax.set_xlabel("vWF intensity in the vessel", rotation = 0, labelpad = 10, fontsize = 20)
    
    #Setting plot title
    ax.set_title('T-cell distance to tumor blood vessels', fontsize=25, fontweight='bold')
    
    #Showing plot
    plt.show()

def histo_distribution():
    '''
    Plot the distribution of T cells around two types of vessels based on the intensity of a marker's staining.
    '''
    #Choosing files
    files = select_files()
    #Computing data
    whole_data_pos = []
    whole_data_neg = []
    histo_data_pos = {}
    histo_data_neg = {}
    for i in range(len(files)):
        #Retrieve proper data format
        file = files[i]
        whole_data_pos.append([])
        whole_data_neg.append([])
        raw_data = collect_data_from_file(file)
        #T cell count correction
        vessel_info, vessel_distance = split_dict_by_info(raw_data)
        position = [pos for pos, char in enumerate(list(vessel_info.keys())[0]) if char == "_"][-1]
        processed_data = remove_doubles(vessel_distance, position)
        final_data = merge_dict_info(vessel_info, processed_data)
        #Split data into the two vessel categories
        for vessel in final_data.keys():
            if final_data[vessel][1] <= 50 and len(final_data[vessel]) >2:
                if len(final_data[vessel]) >3:
                    whole_data_neg[i] += final_data[vessel][2:]
                else:
                    whole_data_neg[i].append(final_data[vessel][2])
            elif 75 < final_data[vessel][1] and len(final_data[vessel]) >2:
                if len(final_data[vessel]) >3:
                    whole_data_pos[i] += final_data[vessel][2:]
                else:
                    whole_data_pos[i].append(final_data[vessel][2])
        #Proper data format for histogram plotting
        if len(whole_data_pos[i])>0:
            whole_data_pos[i] = prep_histo(whole_data_pos[i])
        else:
            #Handling null data
            whole_data_pos[i] = "nan"
        whole_data_neg[i] = prep_histo(whole_data_neg[i])
    for i in reversed(range(len(whole_data_pos))):
        if whole_data_pos[i] == 'nan':
            del whole_data_pos[i]
    
    #Plot parameters: axis labels
    bars = ["Inside", "0-10", "10-20", "20-30", "30-40", "40-50", "+50"]
    bars_num = [1,2,3,4,5,6,7]
    
    #Calculate error bars
    heights_neg = []
    heights_pos = []
    std_neg = []
    std_pos = []
    for i in range(len(bars)):
        heights_neg.append(np.mean([whole_data_neg[k][i] for k in range(len(whole_data_neg))]))
        std_neg.append(np.std([whole_data_neg[k][i] for k in range(len(whole_data_neg))]))
        heights_pos.append(np.mean([whole_data_pos[k][i] for k in range(len(whole_data_pos))]))
        std_pos.append(np.std([whole_data_pos[k][i] for k in range(len(whole_data_pos))]))
    
    #Calculate statistical difference between two distributions
    print(st.ks_2samp(heights_neg, heights_pos))
    
    #Plotting the graph
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (15,15))
    
    #Setting axis label font
    plt.rc("xtick", labelsize = 20)
    plt.rc("ytick", labelsize = 20)
    
    #Plotting the histogram (first category)
    ax1.bar(bars_num, heights_neg, color = (0.5,0.1,0.5,0.6), label = "vWF intensity below 50")
        #Plotting error bars
    ax1.errorbar(bars_num, heights_neg, yerr = std_neg, fmt = "none", color =PALETTE[1])
        #Settings of the plot
    ax1.set_xticks(bars_num)
    ax1.set_xticklabels(bars)
    ax1.spines[["right", "top"]].set_visible(False)
    ax1.set_ylim([0, 0.4])
    ax1.legend(fontsize = 20)
    
    #Plotting the histogram (second category)
    ax2.bar(bars_num, heights_pos, color = (0.1,0.1,0.5,0.6), label = "vWF intensity above 75")
        #Plotting error bars
    ax2.errorbar(bars_num, heights_pos, yerr = std_pos, fmt = "none", color =PALETTE[1])
        #Settings of the plot
    ax2.set_xticks(bars_num)
    ax2.set_xticklabels(bars)
    ax2.spines[["right", "top"]].set_visible(False)
    ax2.set_ylim([0, 0.4])
    ax2.legend(fontsize = 20)
    
    #Setting axis titles
    fig.supylabel("Average fraction of T cells", fontsize = 20)
    fig.supxlabel("T-cell distance to the vessel (µm)", fontsize = 20)
    
    #Setting graph title
    fig.suptitle('T-cell distribution around blood vessels', fontsize=25, fontweight='bold')

    #Show plot
    plt.show()