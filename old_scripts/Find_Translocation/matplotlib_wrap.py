#this script is provide simple wrapper for plot functions

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from matplotlib.ticker import MaxNLocator
from os import path
from datetime import datetime


print("Let me know which data you want to use for making plot (Please enter complete path):")
data = input()
print("Choose one plot type: bar, histogram")
plot_type = input()

fig, ax = plt.subplots()

if path.isfile(data) and plot_type == "bar":
    with open(data) as d_f:
        d_f_data = np.array(d_f.readlines())
        print("File read at %s" % (datetime.now()))
        print("Which column do you want to use for bar plot (write 0 for column 1 and 0 for column 2 and so on)")
        y_column = int(input())
        print("Which column do you want to use for x-axis (write 0 for column 1 and 1 for column 2 and so on)")
        x_column = int(input())
        line_split_tab = np.array([line.split('\t') for line in d_f_data])
        line_split_tab[:, -1] = np.array([l[-1].split('\n')[0] for l in line_split_tab], dtype='object')
        #add threshold for bar values
        #line_split_tab = line_split_tab[0:1000, :]
        data_bar_plot = line_split_tab[:, y_column]
        data_bar_plot = data_bar_plot.astype(float)
        for i in range(len(data_bar_plot)):
            if data_bar_plot[i] > 100:
                data_bar_plot[i] = 100

        print("Please let me know the title of your plot")
        title = input()
        print("Plotting bar plot")
        bar_plot = plt.bar(line_split_tab[:, x_column], data_bar_plot)
        plt.xlabel(line_split_tab[:, x_column])
        plt.title(title)
        print("Saving")
        #bar_plot.xticks()
        plt.savefig(data+"_bar_plot.png")

if path.isfile(data) and plot_type == "histogram":
    with open(data) as d_f:
        d_f_data = np.array(d_f.readlines())
        print("File read at %s" % (datetime.now()))
        print("Which column do you want to use for histogram plot (write 0 for column 1 and 0 for column 2 and so on)")
        x_column = int(input())
        line_split_tab = np.array([line.split('\t') for line in d_f_data])
        line_split_tab[:, -1] = np.array([l[-1].split('\n')[0] for l in line_split_tab], dtype='object')
        # add threshold for bar values
        line_split_tab = line_split_tab[0:1000, :]
        data_plot = line_split_tab[:, x_column]
        data_plot = data_plot.astype(float)
        for i in range(len(data_plot)):
            if data_plot[i] > 100:
                data_plot[i] = 100

        print("Number of bins required in the output")
        num_bins = input()
        print("Please let me know the title of your plot")
        title = input()
        print("Please let me know the title of your plot")
        y_label = input()
        print("Plotting histogram")
        plt.hist(data_plot)
        #ax.plot(bins)
        plt.ylabel(y_label)
        plt.title(title)
        print("Saving")
        # bar_plot.xticks()
        plt.savefig(data + "_hist_plot.png")

