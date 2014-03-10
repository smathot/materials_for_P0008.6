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

import numpy as np
from exparser.Cache import cachedDataMatrix
from exparser.CsvFolderReader import CsvFolderReader

columns = ['subject_nr', 'response_time', 'correct', 'targetRot', 'practice', \
	'postRotDelay', 'cue_dur', 'target_dur', 'soa', 'rot_onset', 'rot_offset']

@cachedDataMatrix
def getDataMatrix():

	"""
	Reads and merges the participant data files. Also selects a subset of
	columns.

	Returns:
	A DataMatrix.
	"""

	dm1 = CsvFolderReader(path='data/1', quote='"').dataMatrix()
	dm1['subject_nr'] += 1000
	dm2 = CsvFolderReader(path='data/2', quote='"').dataMatrix()
	dm2['subject_nr'] += 2000
	dm3 = CsvFolderReader(path='data/3', quote='"').dataMatrix()
	dm3['subject_nr'] += 3000
	dm = dm1 + dm2 + dm3
	dm = dm.selectColumns(columns)
	return dm

@cachedDataMatrix
def filter(dm):

	"""
	Filters the data.

	Arguments:
	dm		--	A DataMatrix.

	Returns:
	A filtered DataMatrix.
	"""

	# Pre-filter N
	print 'N (before) = %d' % dm.count('subject_nr')
	print 'Trials (before) = %d' % len(dm)
	# Remove participants with RT more than 2 SD above/ below the grand mean RT
	rottenApples = []
	cm = dm.collapse(['subject_nr'], 'response_time')
	gm = cm['mean'].mean()
	gstd = cm['mean'].std()
	cm = cm.addField('dev', dtype=float)
	cm['dev'] = np.abs(cm['mean']-gm)/gstd
	rottenApples += list(cm.select('dev > 2')['subject_nr'])
	print 'Rotten apples: %s' % rottenApples
	for ra in rottenApples:
		dm = dm.select('subject_nr != %d' % ra)
	# Remove practice trials
	dm = dm.select('practice == "no"')
	# Post-filter N
	print 'N (after) = %d' % dm.count('subject_nr')
	print 'Trials (after) = %d' % len(dm)
	# Recode conditions
	dm = dm.addField('cond', dtype=str)
	dm = dm.addField('valid', dtype=str)
	dm = dm.addField('logRt', dtype=float)
	dm = dm.addField('iRt', dtype=float)
	dm = dm.addField('rot_dur', dtype=float)
	dm['rot_dur'] = dm['rot_offset'] - dm['rot_onset']
	dm['logRt'] = np.log(dm['response_time'])
	dm['iRt'] = 1 / dm['response_time']
	dm['cond'] = 'object-based'
	dm['cond'][np.where(dm['targetRot'] == 0)[0]] = 'spatial'
	dm['cond'][np.where(dm['targetRot'] == 2)[0]] = 'spatial'
	dm['valid'] = '{0}valid'
	dm['valid'][np.where(dm['targetRot'] == 2)[0]] = '{1}invalid'
	dm['valid'][np.where(dm['targetRot'] == 3)[0]] = '{1}invalid'
	return dm
