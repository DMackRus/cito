# This python file needs to loop through testing folders and creat ssingle csv file for each method

import os
import numpy as np
import matplotlib.pyplot as plt

def box_plot(data, fill_color, yAxisTitle, ax, labels, logyAxis = False, baseline_yLine = False):
    normalPosterColour = "#103755"
    highlightPosterColor = "#EEF30D"


    bp = ax.boxplot(data, patch_artist=True, meanprops={"marker":"s","markerfacecolor":highlightPosterColor, "markeredgecolor":highlightPosterColor}, showmeans=True, showfliers=False)
    if logyAxis:
        ax.set_yscale('log')
    black = "#1c1b1a"

    for element in ['medians']:
        plt.setp(bp[element], color=black)

    for element in ['means']:
        plt.setp(bp[element], color=highlightPosterColor)
        # element.set(color="#a808c4")

    # for bpMean in bp['means']:
    #     bpMean.set(color="#a808c4")

    # for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
    #     plt.setp(bp[element], color=black)


    dynamicColor = "#ff8400"
    baselineColor = "#0057c9"
    setIntervalColour = "#9d00c9"
    backgroundColour = "#d6d6d6"

    
    index = 0
    for patch in bp['boxes']:
        patch.set(facecolor=normalPosterColour)

    labelSize = 11

    ax.set_ylabel(yAxisTitle, fontsize=labelSize)

    if(baseline_yLine):
        ax.axhline(y=0, color=baselineColor, linewidth=1.5, alpha=0.5)

    xticks = []
    for i in range(len(labels)):
        xticks.append(i + 1)

    ax.set_xticks(xticks)
    ax.set_xticklabels(labels, fontsize=11)
        
    return bp 

def main():
    num_trajecs = 100
    methods = ["baseline", "SI2", "SI20", "adaptive_jerk", "magvel_change", "iterative_error"]
    labels = ["optTime", "costReduction", "derivsTime", "qpTime"]

    data = np.zeros((num_trajecs, len(labels) * len(methods)))



    for i in range(len(methods)):
        for j in range(num_trajecs):

            # root = "home/davidrussell/cito_ws/src/cito/src"
            file_name = "../testingData/data/" + methods[i] + "/" + str(j) + ".csv"
            tempData = np.genfromtxt(file_name, delimiter=",")

            data[j, i * len(labels) : (i + 1) * len(labels)] = tempData


    optTimes = np.zeros((num_trajecs, len(methods)))
    costReductions = np.zeros((num_trajecs, len(methods)))
    derivsTime = np.zeros((num_trajecs, len(methods)))
    qpTime = np.zeros((num_trajecs, len(methods)))

    for i in range(len(methods)):
        optTimes[:, i] = data[:, i * len(labels)]
        costReductions[:, i] = data[:, i * len(labels) + 1]
        derivsTime[:, i] = data[:, i * len(labels) + 2]
        qpTime[:, i] = data[:, i * len(labels) + 3]

    orange = "#edb83b"

    fig, axes = plt.subplots(4, 1, figsize = (18,8))
    boxPlotTitle = "opt times"
    yAxisLabel = "opt time (s)"
    orange = "#edb83b"
    bp3 = box_plot(optTimes, orange, yAxisLabel, axes[0], methods, False)

    boxPlotTitle = "cost reduction"
    yAxisLabel = "cost reduction"
    orange = "#edb83b"
    bp4 = box_plot(costReductions, orange, yAxisLabel, axes[1], methods)

    boxPlotTitle = "time derivs"
    yAxisLabel = "derivs time (s)"
    orange = "#edb83b"
    bp4 = box_plot(derivsTime, orange, yAxisLabel, axes[2], methods)

    boxPlotTitle = "time qp"
    yAxisLabel = "qp time (s)"
    orange = "#edb83b"
    bp4 = box_plot(qpTime, orange, yAxisLabel, axes[3], methods)

    meanOptTimes = np.mean(optTimes, axis=0)
    meanCostReductions = np.mean(costReductions, axis=0)
    meanderivsTime = np.mean(derivsTime, axis=0)
    meanQPTimes = np.mean(qpTime, axis=0)

    print("--------------------- methods ----------------------")
    print(methods)
    print("--------------------- mean opt times ----------------------")
    print(meanOptTimes)
    print("--------------------- mean cost reductions ----------------------")
    print(meanCostReductions)
    print("--------------------- mean derivs times ----------------------")
    print(meanderivsTime)




    fig.suptitle("sawyer" + " - derivative information", fontsize=16)
    plt.show()

    # Save data to new csv file
    np.savetxt("../testingData/data.csv", data, delimiter=",")
    


if __name__ == "__main__":
    main()
