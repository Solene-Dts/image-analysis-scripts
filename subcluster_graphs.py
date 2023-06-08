###Introduction

# Hi!

# Characteristics:
    # written for Python 3.8 and matplotlib 3.7
    # contact for help: solenedeutsch@gmail.com

# Aim of the script:
    # Plot relevant graphs for the subcluster analysis

# Getting started:
    # Choose the graph you want to plot
    # Type the corresponding function's name in the command shell
        # Example: histo_marker_intensity("plvap", True)

# Have fun!

### Imports

import ast #Data formatting
import matplotlib.pyplot as plt #Graph plotting
import tkinter as tk #interactive windows
import tkinter.filedialog as fd #retrieving files
import numpy as np #useful maths
import scipy.stats as st #stats
import scikit_posthocs as sp #other stats

### Global variables

# Please feel free to choose better colors...
PALETTE = ["blue", "black", "green", "red", "yellow", "pink", "purple", "orange", "brown", "cyan"]
NICER_PALETTE = [np.array([(3/255, 3/255, 17/255)]), np.array([(62/255, 15/255, 114/255)]), np.array([(122/255, 34/255, 129/255)]), np.array([(167/255, 49/255, 125/255)]), np.array([(215/255, 69/255, 107/255)]), np.array([(242/255, 99/255, 92/255)]), np.array([(254/255, 168/255, 115/255)])]

### Preliminary functions

def select_files():
    '''
    Asks user for input. Return a list with correct file paths.
    '''
    # Retrieving data from text file
    merging = int(input('How many files do you want to combine ? ')) #0 (no combination), or any other number (combination)
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

def select_regions_files():
    '''
    Asks user for input. Return a list with correct file paths.
    '''
    # Retrieving data from text file
    merging = int(input('How many tumors do you want to combine ? (0 or more than 1) ')) #0 (no combination), or any other number (combination)
    root = tk.Tk()
    #Merging files if necessary
    if bool(merging):
        #Making user open files
        files_tum = list(fd.askopenfilenames(title="Choose tumor vessel files to combine"))
        files_rim = list(fd.askopenfilenames(title="Choose rim vessel files to combine"))
        files_out = list(fd.askopenfilenames(title="Choose outside vessel files to combine"))
        files_ctr = list(fd.askopenfilenames(title="Choose contralateral hemisphere vessel files to combine"))
        root.destroy()
    else:
        #Making user open files
        files_tum = [fd.askopenfilename(title="Choose tumor vessel file")]
        #Avoiding tk window freeze
        root.withdraw()
        root = tk.Tk()
        files_rim = [fd.askopenfilename(title="Choose rim vessel file")]
        #Avoiding tk window freeze
        root.withdraw()
        root = tk.Tk()
        files_out = [fd.askopenfilename(title="Choose outside vessel file")]
        #Avoiding tk window freeze
        root.withdraw()
        root = tk.Tk()
        files_ctr = [fd.askopenfilename(title="Choose contralateral hemisphere vessel file")]
        #Avoiding tk window freeze
        root.withdraw()
    return [files_tum, files_rim, files_out, files_ctr]

def collect_intensity(file_title):
    '''
    Opens a text file, retrieves the data as a dictionary and returns lists of values from the dictionary.
    '''
    #Open file
    file = open(file_title, 'r')
    content = file.readlines()
    file.close()
    #Proper data format
    data = ast.literal_eval(content[0])
    #Retrieve interesting data
    int_vwf = [data[i][-2] for i in list(data.keys())]
    int_plvap = [data[i][1] for i in list(data.keys())]
    area = [data[i][0] for i in list(data.keys())]
    area_plvap = [data[i][2] for i in list(data.keys())]
    area_vwf = [data[i][-1] for i in list(data.keys())]
    return area, int_plvap, int_vwf, area_plvap, area_vwf

def data_for_stats(files):
    '''
    Open the results files and retrieve all data for each areas in the brain
    '''
    total_stat = [[] for _ in range(5)]
    for file_title in files:
        file = open(file_title, 'r')
        content = file.readlines()
        file.close()
        data = ast.literal_eval(content[0])
        for vessel in data.keys():
            for i in range(len(data[vessel])):
                total_stat[i].append(data[vessel][i])
    return total_stat

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

### Graph plotting

def graph_both_intensities():
    '''
    Plots the data as a scatter plot : Vwf intensity vs Plvap intensity.
    '''
    files = select_files()
    fig, ax = plt.subplots(1, 1)
    #ax.set_xlim(0, 200)

    tot_vwf = []
    tot_plvap = []
    #Retrieving data for each individual tumor
    for i in range(len(files)):
        file_data = files[i]

        # Retrieving the intensity and area values separately
        area, int_plvap, int_vwf = collect_intensity(file_data)[:3]
        tot_vwf += int_vwf
        tot_plvap += int_plvap

        # Adding to the graph
        ax.scatter(int_plvap, int_vwf, s=4, color = PALETTE[3])

    # Final settings

    #Titles and axes
    fig.suptitle('Vessel phenotype in the tumor', fontsize=17, fontweight='bold')
    fig.supylabel('Vwf intensity')
    fig.supxlabel('Plvap intensity')

    #Showing graph
    plt.show()

def graph_ratioVSarea():
    '''
    Plots the data as a scatter plot: ratio Plvap/Vwf intensity vs vessel area.
    '''
    files = select_files()

    fig, ax = plt.subplots(1, 1, figsize = (10,5))

    #Retrieving data for each individual tumor
    for i in range(len(files)):
        file_data = files[i]

        # Retrieving the intensity and area values separately
        area, int_plvap, int_vwf = collect_intensity(file_data)[:3]
        ratio = []

        for k in range(len(area)):
            plvap = int_plvap[k]
            vwf = int_vwf[k]
            if plvap*vwf == 0:
                ratio.append(0)
            else:
                ratio.append(plvap/vwf)

        # Adding to the graph
        ax.scatter(area, ratio, s=4, color = "blue")
        #ax.scatter(area, int_vwf, s=4, color = PALETTE[i])

    # Final settings

    #Titles and axes
    fig.suptitle('Vessel phenotype in the tumor', fontsize=17, fontweight='bold')
    fig.supylabel('PLVAP/vWF intensity ratio')
    fig.supxlabel('Vessel area')

    #Showing graph
    plt.show()

def histo_marker_intensity(marker, area = False):
    '''
    Plots the average intensity of the marker in different areas of the brain.
    marker can be "plvap" or "vwf".
    area is a parameter (True or False) that indicates if the measurement of the intensity is normal (False)
    or if it requires the area of above-average background (True)
    '''
    #Opening files
    files = select_regions_files()
    print("Now select files from Control brains.")
    files.append(select_files())

    categories = ["TUM", "RIM", "OUT", "CLH", "CTR"]

    #Settings for the data to be plotted
    if marker == "plvap":
        index = 1
        title = "PLVAP distribution across brain areas"
        ylabel = "Percentage of the vessel with above-background signal"
    else:
        index = 2
        title = "vWF distribution across brain areas"
        ylabel = "Intensity of vWF signal inside the vessel"

    if area:
        index+=2

    #Computing data for histogram
    full_histo = [0 for _ in range(5)]
    stock = [[] for _ in range(5)]

    for i in range(len(files[0])):
        file_tum = files[0][i]
        file_rim = files[1][i]
        file_out = files[2][i]
        file_clh = files[3][i]
        if len(files[4]) > i :
            file_ctr = files[4][i]
        else :
            file_ctr = []

        current_files = [file_tum, file_rim, file_out, file_clh, file_ctr]
        histo = [0 for _ in range(5)]

        for type in range(len(current_files)):
            if len(current_files[type]) != 0:
                all_data = collect_intensity(current_files[type])
                histo[type] += sum(all_data[index])/len(all_data[0])
                stock[type].append(sum(all_data[index])/len(all_data[0]))

        for k in range(len(full_histo)):
            full_histo[k] += histo[k]

    for k in range(4):
        full_histo[k] = full_histo[k]/len(files[0])
    full_histo[4] = full_histo[4]/len(files[4])

    #Retrieve proper data for statistical analysis
    stats_tum = data_for_stats(files[0])
    stats_rim = data_for_stats(files[1])
    stats_out = data_for_stats(files[2])
    stats_clh = data_for_stats(files[3])
    stats_ctr = data_for_stats(files[4])

    #Statistical analysis
    stats = [stats_tum[index], stats_rim[index], stats_out[index], stats_clh[index], stats_ctr[index]]
    for stat in stock:
        #Normal distribution test
        print(st.kstest(stat, 'norm'))
    #Homogeneity of variances test
    print(st.bartlett(stock[0], stock[1], stock[2], stock[3], stock[4]))
    #No normal distribution: one-way variance analysis (non-parametric)
    print(st.kruskal(stock[0], stock[1], stock[2], stock[3], stock[4]))

    #Calculate error bars
    histo_std = [np.std(list_values) for list_values in stock]

    categories_num = [1,2,3,4,5]

    #Plotting the graph
    fig, ax = plt.subplots(1, 1, figsize = (15,15))

    #Setting axis label font
    plt.rc("xtick", labelsize = 20)
    plt.rc("ytick", labelsize = 20)

    #Plotting the graph
    ax.bar(categories_num, full_histo, color = NICER_PALETTE[3])
    #Adding error bars
    ax.errorbar(categories_num, full_histo, yerr = histo_std, fmt = "none", color =PALETTE[1])
    #Adding significance TO BE MODIFIED ACCORDING TO THE RESULTS OF THE TESTS
    barplot_annotate_brackets(0, 1, 0.4, categories_num, full_histo, dh = 0.25, fs = 20)
    barplot_annotate_brackets(0, 2, 0.002, categories_num, full_histo, dh = 0.3, fs = 20)
    barplot_annotate_brackets(2, 4, 1, categories_num, full_histo, dh = 0.25, fs = 20)

    #Setting axis labels
    ax.set_xticks(categories_num)
    ax.set_xticklabels(categories)
    ax.spines[["right", "top"]].set_visible(False)

    #Setting axis title
    ax.set_ylabel(ylabel, labelpad = 10, fontsize = 20)
    #Setting graph title
    ax.set_title(title, fontsize=25, fontweight='bold')

    #Showing graph
    plt.show()

def histo_intensity_all_vessels(marker, area = False):
    '''
    Plots all vessel intensity of the marker in different areas of the brain.
    marker can be "plvap" or "vwf".
    area is a parameter (True or False) that indicates if the measurement of the intensity is normal (False)
    or if it requires the area of above-average background (True)
    '''
    files = select_regions_files()

    categories = ["TUM", "RIM", "OUT", "CTR"]
    additional_categories = []

    if marker == "plvap":
        index = 1
        title = "Plvap average intensity"
    else:
        index = 2
        title = "Vwf average intensity"

    if area :
        index += 2

    full_histo = [0 for _ in range(4)]
    stock = []

    for i in range(len(files[0])):
        file_tum = files[0][i]
        file_rim = files[1][i]
        file_out = files[2][i]
        file_ctr = files[3][i]

        current_files = [file_tum, file_rim, file_out, file_ctr]
        histo = [0 for _ in range(4)]

        for type in range(len(current_files)):
            all_data = collect_intensity(current_files[type])
            histo[type] += sum(all_data[index])/len(all_data[0])
            stock.append(all_data[index])
            additional_categories.append([categories[type] for _ in range(len(all_data[index]))])

        for k in range(len(full_histo)):
            full_histo[k] += histo[k]

    for k in range(4):
        full_histo[k] = full_histo[k]/len(files[0])

    fig, ax = plt.subplots(1, 1, figsize = (15,15))

    ax.bar(categories, full_histo)
    for i in range(len(stock)):
        ax.scatter(additional_categories[i], stock[i], s = 25, c = "red")

    fig.suptitle('Average intensity of the markers in different areas of the brain', fontsize=17, fontweight='bold')
    fig.supylabel(title)

    plt.show()