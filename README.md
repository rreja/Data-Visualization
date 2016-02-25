# Data-Visualization
Collection of scripts for genomics data visualization

Introduction
-------------

The scripts below can be used to plot single gene locus screenshot plot and a average plot (composite plot) to visualize the ChIP-exo.

Requirements
------------

- The script requires Python packages: Matplotlib, Scipy, and Numpy.

Running the script
-----------------

To create a ChIP-exo tag distribution around single gene locus, use the following::

    $ python Plot_single_locus.py /usr/local/folder_containing_CDT_files -h
    $
    $ Options:
    $     -h, --help  show this help message and exit
    $     -w WINDOW   Window size of moving average., Default=5
    $
    $ Example output can be seen in the folder test_data/two_strand_cdt/Images


To create an average ChIP-exo tag distribution around reference point, use the following::

    $ python composite_plots.py /usr/local/folder_containing_CDT_files -h
    $
    $ Options:
    $     -h, --help  show this help message and exit
    $     -w WINDOW   Window size of moving average., Default=5
    $
    $ Example output can be seen in the folder test_data/one_strand_cdt/_composite

