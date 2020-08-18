import numpy as np
import matplotlib
matplotlib.use('Agg') # Plot when using ssh -X (avoid server ... error)

import matplotlib.pyplot as plt
import seaborn as sns

colors_bw = ['white', 'whitesmoke', 'lightgray', 'silver', 'darkgrey', 'gray', 'dimgrey', "black"]
colors = ["green", 'blue', 'red', "black", "maroon", "magenta", "cyan"]
linestyles = ['solid', 'dashdot', 'dashed', 'dashed', 'dashdot', 'dotted', 'solid']
linewidths = [1.75, 1.75, 2.5, 2.5, 3.25, 3.75, 2]

def plotTrend(name_to_data, image_file, xlabel, ylabel, order=None):
    if order is None:
        order = list(name_to_data)

    # get median
    plotobj = {name: {'x':list(data.keys()), 'y':[y for _,y in data.items()]} for name, data in name_to_data.items()}

    plt.figure(figsize=(13, 9))
    plt.gcf().subplots_adjust(bottom=0.27)
    #plt.style.use(u'ggplot')
    #sns.set_style("ticks")
    sns.set_style("whitegrid")
    plt.rcParams["axes.edgecolor"] = "0.15"
    plt.rcParams["axes.linewidth"]  = 1.25
    #sns.set_context("talk")
    fontsize = 26
    maxx = max([max(plotobj[t]['x']) for t in order])
    for ti,tech in enumerate(order):
        plt.plot(plotobj[tech]['x'], plotobj[tech]['y'], color=colors[ti], linestyle=linestyles[ti], linewidth=linewidths[ti], label=tech, alpha=0.8)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    step = int(min(maxx, 10))
    plt.xticks(list(range(1, maxx+1, int(maxx/step))) + [maxx] if (maxx % step == 0 or type(maxx) == int) else np.arange(1,maxx+1, maxx/float(step)), fontsize=fontsize-5)
    plt.yticks(np.arange(0,1.01,0.2), fontsize=fontsize-5)
    legendMode=1 if len(order) <= 3 else 2
    if legendMode==1:
        lgd = plt.legend(bbox_to_anchor=(0., 0.98, 1., .102), loc=2, ncol=3, mode="expand", fontsize=fontsize, borderaxespad=0.)
    elif legendMode==2:
        lgd = plt.legend(bbox_to_anchor=(0., 0.98, 1.02, .152), loc=2, ncol=3, mode="expand", fontsize=fontsize, borderaxespad=0.)
    else:
        assert False, "invalid legend mode (expect either 1 or 2)"
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.00), shadow=True, ncol=3)
    #sns_plot.set_title('APFD - '+allkonly)
    plt.tight_layout()
    ybot, ytop = plt.gca().get_ylim()
    ypad = (ytop - ybot) / 50
    #ypad = 2
    plt.gca().set_ylim(ybot - ypad, ytop + ypad)
    plt.savefig(image_file+".pdf", format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close('all')
#~ def plotTrend()

def plot_Box_Grouped(groupedData, imagefile, colors_bw, ylabel, selectGroups=None, selectData=None, \
				groupsLabelsDict=None, dataLabelsDict=None, yticks_range=range(0,101,20), fontsize=26, title=None): 
    plt.figure(figsize=(16, 8)) 
    plt.gcf().subplots_adjust(bottom=0.27)
    #plt.style.use(u'ggplot')
    #plt.style.use(u'seaborn-ticks')
    sns.set_style("ticks")
    
    if selectGroups is None:
        selectGroups = list(groupedData)
    if groupsLabelsDict is None:
        groupsLabelsDict = {g:g for g in selectGroups}
    if selectData is None:
        selectData = list(groupedData[selectGroups[0]])
    if dataLabelsDict is None:
        dataLabelsDict = {d:d for d in selectData}

    nBoxes = 0
    plotobjList = []
    labels = []
    coloring = []
    positions = []
    lastpos = -0.5
    for g_ind, group in enumerate(selectGroups):
        plotobjList += [groupedData[group][data] for data in selectData]
        labels += [' '.join([groupsLabelsDict[group], dataLabelsDict[data]]) for data in selectData]
        coloring += colors_bw[:len(selectData)]
        positions += [lastpos+0.5+i for i in range(1, len(selectData)+1)]
        lastpos = positions[-1]
    nBoxes = len(plotobjList) - g_ind
    bp = plt.boxplot(plotobjList, labels=labels, widths=0.75, positions=positions, patch_artist=True)

    medianValues = []
    for ind,box in enumerate(bp['boxes']):
        box.set(color='black')
        box.set(facecolor = coloring[ind])
    for ind,med in enumerate(bp['medians']):
        med.set(color='black', lw=4)
        medianValues.append(med.get_xydata()[1][1])
    for ind,wh in enumerate(bp['whiskers']):
        wh.set(color='black')
    for ind,wh in enumerate(bp['fliers']):
        wh.set(mew=2)

    plt.ylabel(ylabel, fontsize=fontsize)
    if nBoxes > 2:
        plt.xticks(fontsize=fontsize, rotation=30, ha='right')
    else:
        plt.xticks(fontsize=fontsize) # do not rotate x ticks
    if yticks_range is not None:
        plt.yticks(yticks_range, fontsize=fontsize)
    else:
        plt.yticks(fontsize=fontsize)
    if title is not None:
        plt.title(title, fontsize=fontsize, fontdict={"weight":"bold"})
    plt.tight_layout()
    ybot, ytop = plt.gca().get_ylim()
    ypad = (ytop - ybot) / 50
    #ypad = 2
    plt.gca().set_ylim(ybot - ypad, ytop + ypad)
    plt.savefig(imagefile+".pdf", format='pdf')
    plt.close('all')
    return medianValues
#~ def plot_Box_Grouped()

def plotBoxes(plotobj, order, imagefile, colors_bw, ylabel="APFD", yticks_range=range(0,101,20), fontsize=26, title=None):
    plt.figure(figsize=(16, 8))
    plt.gcf().subplots_adjust(bottom=0.27)
    #plt.style.use(u'ggplot')
    sns.set_style("ticks")
    #fontsize = 26
    plotobjList = [plotobj[t] for t in order]
    bp = plt.boxplot(plotobjList, labels=order, widths=0.75, patch_artist=True)
    medianValues = []
    for ind,box in enumerate(bp['boxes']):
        box.set(color='black')
        box.set(facecolor = colors_bw[ind])
    for ind,med in enumerate(bp['medians']):
        med.set(color='black', lw=4)
        medianValues.append(med.get_xydata()[1][1])
    for ind,wh in enumerate(bp['whiskers']):
        wh.set(color='black')
    for ind,wh in enumerate(bp['fliers']):
        wh.set(mew=2)
        
    plt.ylabel(ylabel, fontsize=fontsize)
    if len(plotobjList) > 3: #2:
        plt.xticks(fontsize=fontsize, rotation=30, ha='right')
    else:
        plt.xticks(fontsize=fontsize) # do not rotate x ticks
    if yticks_range is not None:
        plt.yticks(yticks_range, fontsize=fontsize)
    else:
        plt.yticks(fontsize=fontsize)
    plt.tight_layout()
    ybot, ytop = plt.gca().get_ylim()
    ypad = (ytop - ybot) / 50
    #ypad = 2
    plt.gca().set_ylim(ybot - ypad, ytop + ypad)
    #sns_plot.set_title('APFD - '+allkonly)
    plt.savefig(imagefile+".pdf", format='pdf')
    plt.close('all')
    return medianValues
#~ def plotBoxes()
