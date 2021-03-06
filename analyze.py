#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
This file is part of P0008.6.

P0008.6 is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

P0008.6 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with P0008.6.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
from analysis import helpers, parse

dm = parse.getDataMatrix(cacheId='data')
dm = parse.filter(dm, cacheId='filter')
for arg in sys.argv:
	if hasattr(helpers, arg):
		print '*** Calling helpers.%s() ***' % arg
		retVal = getattr(helpers, arg)(dm)
		if retVal != None:
			dm = retVal
