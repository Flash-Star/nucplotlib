"""
The Torch module provides a class for reading the data produced by Frank 
Timmes' TORCH code (http://cococubed.asu.edu/code_pages/net_torch.shtml).

At the moment, TorchZND is the only class present, written to handle data
produced by the ZND detonation mode in TORCH.

It may very well work with all other modes of TORCH, with little or no 
modifications.

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
from __future__ import print_function
from collections import OrderedDict
import numpy as np

class TorchZND:
	def __init__(self):
		self.data = OrderedDict([])
		self.npdata = None
		self.file_objs = []
		self.header = []

	def readcat(self,fnames):
		self.__init__()
		# fnames: list of Torch output files to concatenate.
		
		for fn in fnames:
			self.file_objs.append(open(str(fn),'r'))
		
		# Consume the headers
		for f in self.file_objs:
			f.readline() # Throwaway empty line 1
			self.header.append(f.readline().strip())
		
		if (self.header == self.header[0]*len(self.header)):
			self.header = [self.header[0]]

		# Create the data fields and read the data
		for f in self.file_objs:
			keys = f.readline().strip().lstrip('*').split()
			dkeys = self.data.keys()
			uniquekeys = []
			for k in keys:
				if k not in dkeys:
					self.data[k] = []
					uniquekeys.append(k)

			## WARNING: TORCH in some file(s) seems to write an incomplete data series at time 0 before
			## writing another series at time 0, this latter one complete. I throw away the first.
			## IF YOU ARE ADAPTING THIS SCRIPT FOR NON-ZND MODES OF TORCH, CHECK THIS BEHAVIOR!!!

			# Read file data into the self.data structure
			fbuffer = []
			for l in f:
				lss = l.strip().split()
				fbuffer.append(lss)
			
			startindx = 0
			if fbuffer[0][0] == fbuffer[1][0]:
				startindx = 1
			
			for values in fbuffer[startindx:]:
				fdict = dict(zip(keys,values[1:]))
				for uk in uniquekeys:
					self.data[uk].append(fdict[uk])

		# Close the input files
		for f in self.file_objs:
			f.close()

	def writecat(self,ofname):
		ofile = open(ofname,'w')
		ofile.write('*  \n')
		ofile.write(self.header[0] + '\n')
		ofile.write('*  ')
		for k in self.data.keys():
			ofile.write(k + '  ')
		ofile.write('\n')

		mlen = []
		for k,v in self.data.iteritems():
			mlen.append(len(v))
		
		if not mlen==[mlen[0]]*len(mlen):
			print('Unequal data array lengths for different keys!')
			exit()
		else:
			mlen = mlen[0]

		for i in xrange(mlen):
			ofile.write(str(i+1) + '  ')
			for k,v in self.data.iteritems():
				ofile.write(v[i] + '  ')
			ofile.write('\n')

		ofile.close()

	def get_numpy_data(self):
		if (self.npdata):
			return self.npdata
		else:
			self.npdata = OrderedDict([])
			for k,v in self.data.iteritems():
				self.npdata[k] = np.array([float(vi) for vi in v])
			return self.npdata
				
			
			
		

		
