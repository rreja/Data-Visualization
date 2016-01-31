import sys, os, math
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib import lines
import matplotlib.cm as cm
import matplotlib.colors as cl
from matplotlib import colorbar as clb
from mpl_toolkits.axes_grid1 import make_axes_locatable


mapping_stats = {}
pp_stats_s5e10 = {}
pp_stats_s20e40 = {}
pp_stats_s20e40_ppno = {}
pp_stats_s5e10_ppno = {}


def parse_mapping_stats(mapFile):
    input = open(mapFile,'rt')
    for line in input:
        data = line.rstrip().split("\t")
        ## Uniquely read map percent
        #mapping_stats[data[5]] = float(data[15][:-1])
        ## Uniquely mapped reads
        data[14] = data[14].replace('"','') # Strings are immutable hence this format.
        data[14] = data[14].replace(",","")
        mapping_stats[data[5]] = math.log(int(data[14]),10)
        
        
def parse_peakpair_stats(ppFile):
    input = open(ppFile,"rt")
    plot =  os.path.join(os.path.dirname(os.path.abspath(ppFile)), "uniq_reads_vs_frip.png")
    ylab = "Fraction of reads in peak-pairs (FRIP) x 10"
    for line in input:
        data = line.rstrip().split("\t")
        
        if os.path.splitext(data[0])[0].split("_")[-1] == "s5e10F1":
            key = os.path.splitext(data[0])[0][:-8]
            pp_stats_s5e10[key[2:]] = float(data[6])*10
            pp_stats_s5e10_ppno[key[2:]] = int(int(data[2])/2)
        elif os.path.splitext(data[0])[0].split("_")[-1] == "s20e40F1":
            key = os.path.splitext(data[0])[0][:-9]
            pp_stats_s20e40[key[2:]] = float(data[6])*10
            pp_stats_s20e40_ppno[key[2:]] = int(int(data[2])/2)
    return(plot,ylab)
                    
                    
def get_xy(ref,new1,new2):
    x = []
    y1 = []
    y2 = []
    count = 0   
    c1 = [] 
    c2 = []
    
    for k,v in ref.items():
        if new1.has_key(k) and new2.has_key(k):
            x.append(v)
            y1.append(new1[k])
            y2.append(new2[k])
            c1.append(pp_stats_s5e10_ppno[k])
            c2.append(pp_stats_s20e40_ppno[k])
            count +=1
    
    print str(count) + " Points plotted"       
    return(x,y1,y2,c1,c2)
            
def plot_graph(x,y1,y2,c1,c2,plot,ylab,LOC):
    
    fig = Figure(figsize=(6,6))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    ax.set_aspect('equal',adjustable='box',anchor='C')
    # Set the title.
    #ax.set_title("XY plot",fontsize=10)
    # Set the X Axis label.
    ax.set_xlabel("Uniquely mapped reads (log10)",fontsize=8)
    # Set the Y Axis label.
    
    ax.set_ylabel(ylab,fontsize=8)   
    # Display Grid.
    #ax.set_xlim(min(x),max(x))
    
    cmap = cm.get_cmap('gist_rainbow')
    cmap.set_over(color='black')
    #cmap.set_under(color='w')
    bounds = [0,1,10,50,200,1000,10000,100000]
    norm = cl.BoundaryNorm(bounds, ncolors=256, clip = False)
    
    
    #ax.grid(True,which='minor',linestyle='-',color='0.75')
    ax.scatter(x,y1,c=c1,cmap=cmap,norm=norm,marker='o',s=5,clip_path=None,label=u"s5e10",edgecolors='none')
    f = ax.scatter(x,y2,c=c2,cmap=cmap,norm=norm,marker='o',s=1,clip_path=None,label=u"s20e40",edgecolors='none')
    
    
    #ax.axhline(y=0.9244,linestyle='--',color='red',linewidth=0.5,label=u"Mac1_scrambled_s20e40")
    #ax.axhline(y=0.0828,linestyle='-',color='red',linewidth=0.5,label=u"Pdc2_scrambled_s20e40")
    #ax.axhline(y=3.895,linestyle='--',color='green',linewidth=0.5,label=u"Reb1-rep2_s20e40")
    
    lfp = FontProperties()
    lfp.set_size('xx-small')
    ax.legend(loc=LOC,prop=lfp)
    
    # Adjust the colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right",size="5%",pad=0.05)
    fig.colorbar(f,cax=cax)
    
    # Save the generated Scatter Plot to a PNG file.
    canvas.print_figure(plot,dpi=500, bbox_inches='tight')
    
    
        
def run():
    x = []
    y = []
    if not len(sys.argv) == 3:
        print "Please provide four input files"
        sys.exit(1)
        
    parse_mapping_stats(sys.argv[1])    
    plot,ylab = parse_peakpair_stats(sys.argv[2])
    x,y1,y2,c1,c2 = get_xy(mapping_stats,pp_stats_s5e10,pp_stats_s20e40)
    plot_graph(x,y1,y2,c1,c2,plot,ylab,2)
    
    
if __name__ == "__main__":
    run() 
