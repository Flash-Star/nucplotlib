"""
nucplot is a script for plotting the results of a 
nuclear reaction network calculation from, e.g. Frank Timmes' TORCH code. 

execute 'python nucplot.py --help' for command line options.

TORCH commonly produces lots of files with a subset of the total number of 
isotopes in each file for the same calculation. Before using these files
as input, you should concatenate them into one file using the readcat and 
writecat functions defined in the Torch class. For an example, see the 
catTorch.py script provided. Use the concatenated file as the input for this 
script.

Copyright 2015 Donald E. Willcox

This file is part of nucplotlib.

    nucplotlib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    nucplotlib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with nucplotlib.  If not, see <http://www.gnu.org/licenses/>.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from Isotopes import *
from Torch import TorchZND
import argparse

# Parse command line parameters
parser = argparse.ArgumentParser()
parser.add_argument("input_file",type=str,
		help="Name of the input file containing time and nuclide tracks.")
parser.add_argument("-o","--outputprefix",type=str,help="Prefix of the output plots (default: 'nuclides'). The time step will be appended with the appropriate file extension.")
parser.add_argument("-ts","--starttime",type=float,help="Start plotting at given simulation time (default: plot all times in dataset).")
parser.add_argument("-te","--endtime",type=float,help="End plotting at given simulation time (default: plot all times in dataset).")
parser.add_argument("-n","--numplots",type=int,help="Generates the given number of plots evenly spaced in time. Cannot be used along with --stride. Note that if you specify numplots so high that the time domain must be subdivided more finely than the data points in the dataset, data points will not be duplicated in the output, so the number of plots saved will be less than numplots.")
parser.add_argument("-s","--stride",type=int,help="Generates plots every specified number of time points (not necessarily evenly spaced!!). Cannot be used along with --numplots. Default stride is 1.")
parser.add_argument("-png","--png",action="store_true",help="Save plots as pngs. (Default is eps).")
parser.add_argument("-dpi","--resolution",type=int,help="Use the specified dpi and output pngs. (Default is eps).")
parser.add_argument("-pdf","--pdf",action="store_true",help="Save plots as pdf. (Default is eps).")
parser.add_argument("-eps","--eps",action="store_true",help="Save plots as eps. (In addition to other output options).")
parser.add_argument("-title","--title",type=str,help="Use the specified plot title text with the time corresponding to the plot. (Default is 'TORCH Nuclides').")
args = parser.parse_args()

# Quickly sanity check the input arguments
if args.stride and args.numplots:
	print 'Error: please use only one of either the --numplots or --stride options!'
	exit()

# Plotting parameters
#n_range = [-0.5,100]
#z_range = [-0.5,100]
box_widths = 1
cmap = cm.get_cmap('Blues')

# Data handling classes
tznd = TorchZND()
nucs = Nuclides()

# Read the TORCH dataset and fill the Nuclides data structure
tznd.readcat([args.input_file])
data = tznd.get_numpy_data()
nucs.load_dataset(data)

# Get plot range from dataset
nz_range = nucs.get_range_nz()
n_range = nz_range['n']
n_range[1] = n_range[1] + 1
z_range = nz_range['z']
z_range[1] = z_range[1] + 1
if n_range[0]==0:
	n_range[0]=-0.5
if z_range[0]==0:
	z_range[0]=-0.5

# Set up the iteration over the dataset by interpreting arguments
def round_time_index(d,t,n):
	# Round the time index down if necessary
	if n == 0:
		return 0
	elif n == len(d['time']):
		return n-1
	elif (t-d['time'][n-1] < d['time'][n]-t):
		return n-1
	else:
		return n

if args.starttime:
	plt_itime_begin = np.searchsorted(data['time'],args.starttime)
	plt_itime_begin = round_time_index(data,args.starttime,plt_itime_begin)	
else:
	plt_itime_begin = 0

if args.endtime:
	plt_itime_end = np.searchsorted(data['time'],args.endtime)
	plt_itime_end = round_time_index(data,args.endtime,plt_itime_end)
else:
	plt_itime_end = len(data['time'])-1

if args.numplots:
	plt_time_values = np.linspace(data['time'][plt_itime_begin],data['time'][plt_itime_end],args.numplots)
	plt_time_indices = np.searchsorted(data['time'],plt_time_values)
	for pt in range(len(plt_time_values)):
		plt_time_indices[pt] = round_time_index(data,plt_time_values[pt],plt_time_indices[pt])

elif args.stride:
	plt_time_indices = range(plt_itime_begin,plt_itime_end,args.stride)
	plt_time_indices.append(plt_itime_end)
else:
	plt_time_indices = range(plt_itime_begin,plt_itime_end)
	plt_time_indices = append(plt_itime_end)

if args.outputprefix:
	plt_prefix = args.outputprefix
else:
	plt_prefix = 'nuclides'

if args.title:
	plt_title = args.title
else:
	plt_title = 'TORCH Nuclides'

# Determine how many decimal places we need to represent the time indices
plt_dec_places = int(np.floor(np.log10(plt_time_indices[-1])))+1

# Iterate over the time steps in the dataset and plot
for t_n in plt_time_indices:
	# Setup the plot
	fig = plt.gcf()
	fig.clf()

	axes = fig.add_axes([0.1,0.1,0.8,0.8])
	# Plot isotope grid
	for i in nucs.isodata:
		square = plt.Rectangle((i.n-0.5*box_widths,i.z-0.5*box_widths),
				box_widths,box_widths,facecolor=cmap(i.x[t_n]),
				edgecolor='black')
		axes.add_patch(square)
	
	# Set Plot Limits
	plt.ylim(z_range)
	plt.xlim(n_range)
	
	# Labels/Titles and Show/Save
	plt.title(plt_title + ', t = '+'{0:0.6e}'.format(data['time'][t_n]) 
			+ ' s')
	plt.xlabel('N')
	plt.ylabel('Z')
	
	# Setup colorbar
	axcb = mpl.colorbar.make_axes(axes,fraction=0.05)[0]
	cbar = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=mpl.colors.Normalize(vmin=0.0, vmax=1.0))
	
	#plt.show()
	t_n_spec = '{0:>0' + str(plt_dec_places) + '}'
	t_n_filename = t_n_spec.format(t_n)
	if args.resolution:
		plt.savefig(plt_prefix + '_'+t_n_filename+'.png',dpi=args.resolution)
	elif args.png:
		plt.savefig(plt_prefix + '_'+t_n_filename+'.png',dpi=300)
	if args.pdf:
		plt.savefig(plt_prefix + '_'+t_n_filename+'.pdf')
	if args.eps or not (args.resolution or args.png or args.pdf):
		plt.savefig(plt_prefix + '_'+t_n_filename+'.eps')

print "Plotting complete!"
print "If you would like to combine pngs into a video file (e.g. mp4), "
print "and have ffmpeg installed, consider using a command such as:"
print "./ffmpeg -framerate 30 -pattern_type glob -i 'nuclides_*.png' -c:v libx264 -r 30 -pix_fmt yuv420p nuclides.mp4"
