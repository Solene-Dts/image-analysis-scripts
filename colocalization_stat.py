##Intro

#Date : 23/02/23
#Python version : 3.8

#This script assumes the existence of a file containing a dictionary where the keys are images names and the values are a list of relevant stats

#Expected image name format : P1*1 where P1 is the name of the panel, 1 is the number of the image taken for that panel, and * represents any character exept a blank space


#BEFORE STARTING :
    #Modify file name and path name to fit your configuration (lines 26, 28 and 35)

##Importations

import os #handles paths and working directories
import pandas as pd #handles dataframe format (tables)
import ast #special formatting
import matplotlib.pyplot as plt #handles graphs
import numpy as np #useful math

##Core script


def data_handling(file_name, number):
    #retrieve saved data

    file = open(file_name, "r")
    data = file.readlines()
    res = ast.literal_eval(data[0])
    file.close()

    #create and edit dataframe from dictionary
    df = pd.DataFrame(res.items())
    df.columns = ['Image', 'Values'] #column names
    df['Values'] = df['Values'].apply(str)
    df[['tM1', 'tM2']] = df['Values'].str.split(",", expand = True) #split the list to give the separate values

    df['tM1'] = df['tM1'].apply(list) #easiest way to edit a string
    df['tM2'] = df['tM2'].apply(list)

    for i in range(len(df['Values'])): #keeps only the numbers
                                        #there has to be a better way to do this
        df['tM1'][i][1] =''
        df['tM1'][i][0] =''
        df['tM1'][i][-1] =''
        df['tM2'][i][-1] =''
        df['tM2'][i][1] =''
        df['tM2'][i][0] =''
        df['Image'][i] = str(number) + df['Image'][i]

    df['tM1'] = df['tM1'].apply(''.join) #and back to a string
    df['tM2'] = df['tM2'].apply(''.join)

    df.drop('Values', inplace=True, axis=1)

    df['tM1'] = df['tM1'].apply(float) #and finally we get proper number format
    df['tM2'] = df['tM2'].apply(float)

    #get the panel names (yes this could have been automated)

    type = input("Do you want to plot each vessel separately? (y/n) : ")
    if type == "y":
        list_image = list(df['Image'])
    else:
        list_image = list(str.split(input("Name of panels again please (expected format is P1spaceP2) : "), " "))

    #another dataframe to stock relevant statistic data
    stats = pd.DataFrame()
    whole = pd.DataFrame()

    #fill the dataframe for each panel
    for panel in list_image:
        boolean = []
        for i in range(len(df['Image'])):
            if panel in df['Image'][i]:
                boolean.append(True)
            else :
                boolean.append(False)
        #temporary dataframe which contains the dataframe information about all the images of a certain panel
        intermediate = df.iloc[boolean]
        whole = pd.concat([whole, intermediate])
        mean_M1 = intermediate['tM1'].mean() #calculates the mean
        mean_M2 = intermediate['tM2'].mean()
        stdev_M1 = intermediate['tM1'].var()**0.5 #calculates the standard deviation
        stdev_M2 = intermediate['tM2'].var()**0.5
        #add proper data to the stats dataframe
        stats[panel] = [mean_M1, mean_M2, stdev_M1, stdev_M2]

    #adds row names
    stats.index = ['Mean tM1', 'Mean tM2', 'St-Dev tM1', 'St-Dev tM2']
    return list_image, stats, whole

def make_graph(list_image, stats):
    #plots the graph : here shows tM1 and tM2 value with standard deviation as error bars

    X_axis = np.arange(len(list_image))

    Title = input("Title of the graph : ")

    #to add for error bars : , yerr= list(stats.loc['St-Dev tM1'])

    #Plotting the histogram

    fig, (ax1) = plt.subplots(1, 1, figsize = (8,8))
    fig.subplots_adjust(left = 0.09, wspace = 0.5, hspace = 0.1)  # adjust space between axes

    ax1.bar(X_axis+0.2, list(stats.loc['Mean tM1']), 0.4, color = "green", label = 'tM1')
    ax1.bar(X_axis-0.2, list(stats.loc['Mean tM2']), 0.4, color = "blue", label = 'tM2')

    # Add thresholds

    #ax1.axhline(3.5, ls='-', color='purple', label = 'Tumor-specific tM1 threshold')

    ax1.legend()

    fig.suptitle(Title, fontsize=17, fontweight='bold')
    fig.supylabel('Marker coefficients for colocalization analysis')
    fig.supxlabel('Vessel')
    plt.show()

def define_thresholds(whole):
    print("Pay attention you only selected the control panels before")
    print("Max of tM1 is :")
    print(max(list(whole['tM1'])))
    print("Max of tM2 is :")
    print(max(list(whole['tM2'])))
    return max(list(whole['tM1'])), max(list(whole['tM2']))

def run_mode(file_list = []):
    if len(file_list) ==0:
        os.chdir("E:\ENS Lyon\IGP") #replace with your own working directory
        file_list.append(input("Enter file name : "))
    mode = int(input("Type 0 for graph mode \nType 1 for threshold mode\n-> "))
    total_per = {}
    for i in range(len(file_list)):
        print(file_list[i])
        list_image, stats, whole = data_handling(file_list[i], i)
        if mode ==0:
            make_graph(list_image, stats)
        elif mode ==1:
            define_thresholds(whole)
    if len(total_per) >0:
        return total_per

run_mode()


