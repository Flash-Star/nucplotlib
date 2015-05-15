# nucplotlib
Provides classes for storing and plotting the results of a nuclear reaction network calculation from, e.g. TORCH.

An example data input file is provided by cat_wn_0.dat. To subdivide the time interval (seconds) [0,0.02]
into 10 timesteps and plot those, run the following command:
>python nucplot.py -te 1.0e-2 -n 10 -png cat_wn_0.dat

To see the available command line options, run:
>python nucplot.py --help 
