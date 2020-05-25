import matplotlib
matplotlib.use('Agg') # Plot when using ssh -X (avoid server ... error)

import matplotlib.pyplot as plt
import seaborn as sns

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
    plt.xticks(range(1, maxx+1, maxx/10) if (maxx % 10 == 0 or type(maxx) == int) else np.arange(1,maxx+1, maxx/10.0), fontsize=fontsize-5)
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
