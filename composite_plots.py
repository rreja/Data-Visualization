import sys, os, re, difflib, csv
from operator import add
from itertools import izip, cycle, tee
from optparse import OptionParser , IndentedHelpFormatter
from pylab import *
import numpy as np
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as plticker
from scipy import stats


### color schema for Figure 2A
#colors = ["#FBDB0C","#FF00FF","#0000FF","#336600","#FF0000","#000000","#00FF00","#660066"] #"#663300", "#660066", "#009999"

## color schema for Figure 5A, the following order: Sua7, Taf10, Taf12, Taf1,Taf2,Taf4,Taf5,Taf8,Toa2
#colors = ["#00FF00","#CC00FF","#FFCC00","#0000FF","#006600","#CC0000","#FF9900","#66CCFF","#00FFFF", "#663300", "#660066", "#009999"]


# Setting 20 colors
# These are the "Tableau 20" colors as RGB.    
colors = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(colors)):    
    r, g, b = colors[i]    
    colors[i] = (r / 255., g / 255., b / 255.)  


# Color schema for Figure 3D, Sua7 shift
#colors = ["#000000","#CFD4CF","#CFD4CF","#CFD4CF","#00FF00","#A3FFC2","#A3FFC2"]
#colors = ["#FF3399","#FF0000","#CFD4CF","#CFD4CF","#CFD4CF","#A3FFC2","#A3FFC2","#000000","#00FF00"]
#colors = ["#66CCFF","#0000FF","#FF3399","#FF0000","#00CC00","#00FF00"]
#colors = ["#FF3399","#FF0000","#00CC00","#00FF00"]
#colors = ["#00CC00","#00FF00"]

#colors = ["#0000FF","#336600","#000000","#00FF00","#660066"]
### color schema for Figure 3E
#colors = ["#0000FF","#336600","#000000"]

## H2A color scheme. fig S2
#colors = ["#FFd699","#FF9900","r"]
# H2B
#colors = ["#EBAD99","#CC3300","r"]
#H3
#colors = ["#B2D1FF","#0066FF","r"]
#H4
#colors = ["#B2E0C2","#009933","r"]

### color schema for HS/MHS panel
#colors = ["#0000FF","#66CCFF","#009999"]
#colors = ["#000000","#FF0000"]

# Alwins figure color schme
#colors = ["#FF3399","#66CCFF","r","b"]
#colors = ["#00FF00","b"]
#colors = ["b","r"]


list1 = []
def  process_files(sense,anti,options,output_folder,ax,count,color_count):
    print "processing "+sense+" and "+anti
    #print count
    
    label = os.path.basename(sense).split("_")[0]
    X =[]
    # Process Sense file first
    in_sense = open(sense,"rt")
    for line in in_sense:
        if line.startswith("Uniqe"):# or line.startswith("ID"):
            tmp = line.rstrip().split("\t")[2:]
            X = [int(x) for x in tmp]
            xmin = min(X)
            xmax = max(X)
            Y1 = [0]*len(X)
            continue
        if line.startswith("EWEIGHT"):
            continue
        tmplist = line.rstrip().split("\t")[2:]
        newList = [float(x) for x in tmplist]
        Y1 = map(add,Y1,newList)
    
    # Process anti-sense file
    in_anti = open(anti,"rt")
    for line in in_anti:
        if line.startswith("Uniqe"):# or line.startswith("ID"):
            #X = line.rstrip().split("\t")[2:]
            Y2 = [0]*len(X)
            continue
        if line.startswith("EWEIGHT"):
            continue
        tmplist = line.rstrip().split("\t")[2:]
        newList = [float(x) for x in tmplist]
        Y2 = map(add,Y2,newList)
    
  

    #### To print out the X and Y coordiantes.
    out1 = open(os.path.join(output_folder,label+"_sense.txt"),"w")
    out2 = open(os.path.join(output_folder,label+"_anti.txt"),"w")
    
    for x,y in zip(X,Y1):
        if x >= -200 and x <= 200:
            out1.write(str(x)+"\t"+str(y)+"\n")
            
    for x,y in zip(X,Y2):
        if x >= -200 and x <= 200:
            out2.write(str(x)+"\t"+str(y)+"\n")
            
    
    
    plot_graph(X,Y1,Y2,xmin,xmax,options,ax,label,count,color_count)
    ##plot_graph(X,newY1,newY2,xmin,xmax,options,count)
    

def process_onestrand_files(infile,options,output_folder,ax,count,color_count):
    print "processing "+infile
    X = []
    label = os.path.basename(infile).split("_")[0]
    # Process the only CDT file 
    in_sense = open(infile,"rt")
    line_count = 0
    for line in in_sense:
        if line.startswith("Uniqe") or line.startswith("ID") or line.startswith("gene"):
            tmp = line.rstrip().split("\t")[2:]
            X = [int(x) for x in tmp]
            xmin = min(X)
            xmax = max(X)
            Y = [0]*len(X)
            continue
        if line.startswith("EWEIGHT") or line.startswith("chr1:1-100"):
            continue
        
        #tmplist = line.rstrip().split("\t")[602:1402]
        tmplist = line.rstrip().split("\t")[2:]
        newList = [float(x) for x in tmplist]
        #print len(Y)
        #print len(newList)
        Y = map(add,Y,newList)
        line_count = line_count + 1
    ##list1.append(smoothListGaussian(Y,50))
    #print line_count
    #Y = [x/line_count for x in Y]
    #Y = [float(x)/max(Y) for x in Y]
    
    ####### Create excel files for each factor #######
    #out3 = open(os.path.join(output_folder,label+".txt"),"w")
    #for x,y in zip(X,Y):
    #    if x >= -500 and x <= 500:
    #        out3.write(str(x)+"\t"+str(y)+"\n")
    #######
    plot_graph(X,Y,0,xmin,xmax,options,ax,label,count,color_count)
    

def plot_graph(X,Y1,Y2,xmin,xmax,options,ax,label,count,color_count):
    
    #print "The X-coordinate at the max Y value = "+str(X[Y1.index(max(Y1))])+" "+str(X[Y2.index(max(Y2))])
    # Converting the arrays to numpy array
    #out = open(os.path.join("/Users/rohitreja/Desktop",label+".txt"),"w")
    X = np.array(X)
    Y1 = np.array(Y1)
    if not Y2 == 0:
        
        Y2 = np.array(Y2)
        ## Smoothing using moving average, REQUIRED WHEN PLOTTING COMPOSITE ON +Y AND  -Y
        X = movingaverage(X,options.window)
        Y1 = movingaverage(Y1,options.window)
        Y2 = movingaverage(Y2,options.window)
        Y1 = [float(x)/max(Y1) for x in Y1]
        Y2 = [float(x)/max(Y2) for x in Y2]
        Y2 = [-1*x for x in Y2]
            
        ## Smoothing using moving average, GOOD FOR COMBINED COMPOSITE.
        ##X = smoothListGaussian(X,options.window)
        ##Y1 = smoothListGaussian(Y1,options.window)
        ##Y2 = smoothListGaussian(Y2,options.window)
        
        #Y3 = map(add, Y1, Y2)
        #Y3 = [float(x)/np.median(Y3) for x in Y3]
        #for i in zip(X,Y3):
        #    out.write(str(i[0])+"\t"+str(i[1])+"\n")
        #out.close()
        
        
        if label == "plus1":
        #    print "Plotting Hmo1 on right axis"
            ax.plot(X, Y1, color='#AAAAAA',label=label,lw=3.0)
            ax.plot(X, Y2, color='#AAAAAA',lw=3.0)
            ax.fill_between(X,Y1,0,facecolor='#AAAAAA',edgecolor='#AAAAAA')
            ax.fill_between(X,Y2,0,facecolor='#AAAAAA',edgecolor='#AAAAAA')
            #ax.set_xlim(-80,80)
        #    for tick in ax2.yaxis.get_major_ticks():
        #        tick.label.set_fontsize(12)
        #    
        else:
        #    
        #    ax.plot(X, Y3, color=colors[count],label=label,lw=2.0)
            ax.plot(X, Y1, color=colors[count],label=label,lw=3.0)
            ax.plot(X, Y2, color=colors[count],lw=3.0)
            #ax.set_xlim(-80,80)
        
        
        
        ##ax[count-1].fill_between(X,+Y1,0,facecolor='#9999ff',edgecolor='#9999ff')
        ##ax[count-1].fill_between(X,-Y2,0,facecolor='#ff9999',edgecolor='#ff9999')
    
        #ax.set_xlim(xmin,xmax)
        #ax.set_xlim(-80,80)
        # Uncomment here to make y-axis same for all graphs
        #ax[count-1].set_ylim(-500,300)
        
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)
            
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
    else:
        Y1 = movingaverage(Y1,options.window)
        Y1 = [float(x)/max(Y1) for x in Y1]
       
        X = movingaverage(X,options.window)
        
        
        ## To print the area under the curve
        #print sum(Y1)
        #Y1 = [float(x)/max(Y1) for x in Y1]
        
        
        #index = Y1.index(max(Y1)) - 500
        #print "Max value at "+str(index)
        #X = smoothListGaussian(X,options.window)
        #Y1 = smoothListGaussian(Y1,options.window)
        
        
        #if label == "Hmo1":
        #    
        #    X = smoothListGaussian(X,30)
        #    Y1 = smoothListGaussian(Y1,30)
        #    Y1 = [float(x)/max(Y1) for x in Y1]
        #    ax.plot(X, Y1,lw=2.0,color='#AAAAAA')
        #    ax.fill_between(X,Y1,0,facecolor='#AAAAAA',edgecolor='#AAAAAA')
        #    ax.set_xlim(-50,350)
        #    
        #if label == "Normal":
        #    
        #    ax.plot(X, Y1,lw=3.0,color='#AAAAAA')
        #    ax.fill_between(X,Y1,0,facecolor='#AAAAAA',edgecolor='#AAAAAA')
        #    ax.set_xlim(-500,1000)
        #    
        #    
        #
        #else:
        
        if label == "Nap1-IF-140-180-Rep2" or label == "plus1":
            ax.plot(X, Y1, color="#AAAAAA",label=label,lw=3.0)
            ax.fill_between(X,Y1,0,facecolor='#AAAAAA',edgecolor='#AAAAAA')
        
        #elif label == "H3-MHS":
        #    ax.plot(X, Y1, color="blue",label=label,lw=3.0)
        #    ax.fill_between(X,Y1,0,facecolor='blue',edgecolor='blue')
            #ax.set_xlim(-500,500)
            #ax.set_ylim(0,1)
        else:
            ax.plot(X, Y1, color=colors[count],label=label,lw=3.0)
            #ax.set_xlim(-250,250)
            #ax.set_ylim(0,1)
        ax.tick_params('x', length=10, width=1, which='major')
        #ax1.tick_params('both', length=10, width=1, which='minor')
        #ax.plot(X, Y1, color="#AAAAAA",label=label,lw=3.0)
        #ax.fill_between(X,Y1,0,facecolor='#AAAAAA',edgecolor='#AAAAAA')
        ax.set_xlim(-100,100)
        #ax.set_ylim(0,1)
            
        #ax.set_ylim(0,1)
        #ax.axes.get_xaxis().set_visible(False)
        #ax.yaxis.tick_right()
        #ax.yaxis.set_ticks_position('both')
        #loc = plticker.MultipleLocator(base=100.0)
        #ax.yaxis.set_major_locator(loc)
        #ax.set_ylim(0,max(Y1))
        #for tick in ax.yaxis.get_major_ticks():
        #    tick.label.set_fontsize(8)
        #    
        #for tick in ax.xaxis.get_major_ticks():
        #    tick.label.set_fontsize(8)
        #
    
def movingaverage(interval, window_size):
    # suggestion from: http://argandgahandapandpa.wordpress.com/2011/02/24/python-numpy-moving-average-for-data/
    #extended_interval = np.hstack([[interval[0]] * (window_size- 1), interval])
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'valid')
    #return np.convolve(extended_interval, window)[window_size-1:-(window_size-1)]


def smoothListGaussian(list,degree):  

     window=degree*2-1  
     weight=np.array([1.0]*window)  
     weightGauss=[]  
     for i in range(window):  
         i=i-degree+1  
         frac=i/float(window)  
         gauss=1/(np.exp((4*(frac))**2))  
         weightGauss.append(gauss)  
     weight=np.array(weightGauss)*weight  
     smoothed=[0.0]*(len(list)-window)  
     for i in range(len(smoothed)):  
         smoothed[i]=sum(np.array(list[i:i+window])*weight)/sum(weight)
     #print smoothed
     new_smooth = [int(i) for i in smoothed]
     return new_smooth  





usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python composite_for_many_factors_one_plot.py /usr/local/folder_containing_CDT_files
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-w', action='store', type='int', dest='window', default = 5,
                      help='Window size of moving average., Default=5')
    #parser.add_option('-d', action='store', type='int', dest='down',default=20,
    #                  help='Downstream distance to go from the peak-pair mid_point, Default=20')
    #parser.add_option('-p', action='store', type='string', dest='ppDir',
    #                  help='Directory containing peak-pairs.')
    #parser.add_option('-g', action='store', type='string', dest='genome_file',
    #                  help='File with chromosome lengths.')
    #parser.add_option('-i', action='store', type='int', dest='ilength', default = 40,
    #                  help='interval length., Default=40')

    (options, args) = parser.parse_args()
    
    if not args:
        parser.print_help()
        sys.exit(1)
        
    antisense_files = []
    sense_files = []
    
    output_folder = os.path.join(args[0],'_composite/') 
    if not os.path.exists(output_folder): os.makedirs(output_folder)
    outfile = os.path.join(output_folder,"composite_plot_all_factors.svg")
    
    
    if not os.path.exists(args[0]):
        parser.error('Path %s does not exist.' % args[0])
    if os.path.isdir(args[0]):
        for fname in os.listdir(args[0]):
            if fname.endswith("cdt"): # or fname.endswith("txt"):
                # intellignetly join paths without worrying about '/'
                fpath = os.path.join(args[0], fname)
                matchobj = re.match(r'(.*)anti(.*)',fname)
                if matchobj:
                    antisense_files.append(fpath)
                else:
                    sense_files.append(fpath)
                    
        # No of subplots to plot
        nof = len(sense_files)
         # Declaring plotting parameters

        f,ax = matplotlib.pyplot.subplots(1,1,sharex=True)
        count = -1
        color_count = 0
        
        
        if len(antisense_files) == 0:
            for f in sense_files:
                count = count + 1
                process_onestrand_files(f,options,output_folder,ax,count,color_count)
                color_count = color_count + 2
            #print stats.ks_2samp(list1[0],list1[1])
            lfp = FontProperties()
            lfp.set_size(12)
            ax.legend(prop=lfp)
            #ax.grid()
            savefig(outfile)
        else:   
            for f1 in antisense_files:
                count = count + 1
                Oratio = 0
                file1_anti = ""
                file2_sense = ""
                for f2 in sense_files:
                    # Calculating the distance ration between two file names and grouping together based upon the highest ratio
                    # Requirement: Antisense file should have anti word in it and both sense and antisense files should have similar names
                    Nratio = difflib.SequenceMatcher(None,f1,f2).ratio()
                    if(Nratio >= Oratio):
                        file1_anti = f1
                        file2_sense = f2
                        Oratio = Nratio
                    
                process_files(file2_sense,file1_anti,options,output_folder,ax,count,color_count)
                color_count = color_count + 2
            lfp = FontProperties()
            #lfp.set_size('small')
            lfp.set_size(12)
            ax.legend(prop=lfp)
            #ax.grid()
            savefig(outfile)
    
    
if __name__ == "__main__":
    run() 