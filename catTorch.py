"""
This script will concatenate all torch files piped into it into 'cat_torch.dat'

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

from Torch import TorchZND
from sys import argv

# Pipe the torch files to concatenate into this script.
torch_files = argv[1:]

tznd = TorchZND()
tznd.readcat(torch_files)
tznd.writecat("cat_torch.dat")

