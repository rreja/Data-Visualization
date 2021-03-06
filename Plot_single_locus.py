import sys, os, re, difflib, csv
from operator import add
from itertools import izip, cycle, tee
from optparse import OptionParser , IndentedHelpFormatter
from pylab import *
import numpy as np


# Color schema, Good for upto plotting 6 factors in a single plot, add more colors 
colors = ["#CC0000","#FF0000","#CC00CC","#FF00FF","#336600","#339900","#000000","#666666","#3399FF","#ADD6FF","#336600","#99B280"]


def process_twostrand_cdt_files(sense,anti,options,output_folder,count,ax,color_count):
    print "processing "+sense+" and "+anti
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
    
    # Process Anti-sense file
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
    
    plot_graph(X,Y1,Y2,xmin,xmax,options,count,ax,color_count)
    
    

def process_onestrand_cdt_files(infile,options,output_folder,count,ax,color_count):
    print "processing "+infile
    X = []
    # Process the CDT file 
    in_sense = open(infile,"rt")
    for line in in_sense:
        if line.startswith("Uniqe") or line.startswith("ID"):
            tmp = line.rstrip().split("\t")[2:]
            X = [int(x) for x in tmp]
            xmin = min(X)
            xmax = max(X)
            Y = [0]*len(X)
            
            continue
        if line.startswith("EWEIGHT"):
            continue
        tmplist = line.rstrip().split("\t")[2:]
        newList = [float(x) for x in tmplist]
        Y = map(add,Y,newList)
    #print X,Y
        
    plot_graph(X,Y,0,xmin,xmax,options,count,ax,color_count)
    

def plot_graph(X,Y1,Y2,xmin,xmax,options,count,ax,color_count):
    
    #print "The X-coordinate at the max Y value = "+str(X[Y1.index(max(Y1))])+" "+str(X[Y2.index(max(Y2))])
    # Converting the arrays to numpy array
    X = np.array(X)
    Y1 = np.array(Y1)
    if not Y2 == 0:
        
        Y2 = np.array(Y2)
        # Smoothing using moving average.
        #Y1 = movingaverage(Y1,options.window)
        #Y2 = movingaverage(Y2,options.window)
        
        #print Y1
        # Smoothing using gaussian
        X = smoothListGaussian(X,options.window)
        Y1 = smoothListGaussian(Y1,options.window)
        #print Y1
        Y2 = smoothListGaussian(Y2,options.window)
        
        Y2 = [x * -1 for x in Y2] 
        #ax[count-1].plot(X, Y1, color='#9999ff')
        #ax[count-1].plot(X, Y2, color='#ff9999')
        ax[count-1].plot(X, Y1, color=colors[color_count])
        ax[count-1].plot(X, Y2, color=colors[color_count+1])
        #ax[count-1].fill_between(X,Y1,0,facecolor='#9999ff',edgecolor='#9999ff')
        #ax[count-1].fill_between(X,Y2,0,facecolor='#ff9999',edgecolor='#ff9999')
        ax[count-1].fill_between(X,Y1,0,facecolor=colors[color_count],edgecolor=colors[color_count])
        ax[count-1].fill_between(X,Y2,0,facecolor=colors[color_count+1],edgecolor=colors[color_count+1])
    
        ax[count-1].set_xlim(xmin,xmax)
        # Uncomment here to make y-axis same for all graphs
        #ax[count-1].set_ylim(0,100)
        #ax[count-1].set_ylim(-100,100)
        
        for tick in ax[count-1].yaxis.get_major_ticks():
            tick.label.set_fontsize(8)
            
        for tick in ax[count-1].xaxis.get_major_ticks():
            tick.label.set_fontsize(8)
    else:
        
        X = movingaverage(X,options.window)
        Y1 = movingaverage(Y1,options.window)
        
        ax[count-1].plot(X, +Y1, color=colors[color_count])
        ax[count-1].fill_between(X,+Y1,0,facecolor=colors[color_count],edgecolor=colors[color_count])
        ax[count-1].set_xlim(xmin,xmax)
        #ax[count-1].set_ylim(0,100)
        
        
        for tick in ax[count-1].yaxis.get_major_ticks():
            tick.label.set_fontsize(8)
            
        for tick in ax[count-1].xaxis.get_major_ticks():
            tick.label.set_fontsize(8)
    
    
    


def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'valid')


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
python composite_plots_for_many_factor.py /usr/local/folder_containing_CDT_files
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
    
    output_folder = os.path.join(args[0],'Images/') 
    if not os.path.exists(output_folder): os.makedirs(output_folder)
    outfile = os.path.join(output_folder,"composite_combined.svg")
    
    
    if not os.path.exists(args[0]):
        parser.error('Path %s does not exist.' % args[0])
    if os.path.isdir(args[0]):
        for fname in os.listdir(args[0]):
            if fname.endswith("cdt"):
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
        f,ax = subplots(nof,1,sharex=True)
        f.subplots_adjust(hspace=0)
        count = 0
        color_count = 0
        
        
        if len(antisense_files) == 0:
            for f in sense_files:
                count = count + 1
                
                process_onestrand_cdt_files(f,options,output_folder,count,ax,color_count)
                color_count = color_count + 1
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
                    
                process_twostrand_cdt_files(file2_sense,file1_anti,options,output_folder,count,ax,color_count)
                color_count = color_count + 2
            savefig(outfile)
    
    
if __name__ == "__main__":
    run() 