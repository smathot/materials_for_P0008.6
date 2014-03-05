#!/usr/bin/env python

import sys
import numpy as np
from scipy.stats import ttest_rel
from exparser import Constants
from exparser.DataMatrix import DataMatrix
from exparser.PivotMatrix import PivotMatrix
from exparser.AnovaMatrix import AnovaMatrix
from exparser.CsvFolderReader import CsvFolderReader
from matplotlib import pyplot as plt
import helpers

columns = ['subject_nr', 'response_time', 'correct', 'targetRot', 'practice', \
	'postRotDelay', 'cue_dur', 'target_dur', 'soa', 'rot_onset', 'rot_offset']
stdF = 2.5
Constants.fontFamily = 'arial'
Constants.fontSize = 10
Constants.plotLineStyles = ['-', ':']
Constants.plotLineColors = ['black'] * 2
Constants.plotLineWidth = 1
Constants.palette = ['#ffffff', '#777777']

if '--parse' in sys.argv:
	dm1 = CsvFolderReader(path='data/1', quote='"').dataMatrix()
	dm1['subject_nr'] += 1000
	dm2 = CsvFolderReader(path='data/2', quote='"').dataMatrix()
	dm2['subject_nr'] += 2000
	dm3 = CsvFolderReader(path='data/3', quote='"').dataMatrix()
	dm3['subject_nr'] += 3000
	dm = dm1 + dm2 + dm3
	dm = dm.selectColumns(columns)
	dm.save('data.npy')
else:
	dm = DataMatrix('data.npy')

# The effective N
print 'N (before) = %d' % dm.count('subject_nr')

# Remove participants with RT more than 2 SD above/ below the grand mean RT
rottenApples = []
cm = dm.collapse(['subject_nr'], 'response_time')
gm = cm['mean'].mean()
gstd = cm['mean'].std()
cm = cm.addField('dev', dtype=float)
cm['dev'] = np.abs(cm['mean']-gm)/gstd
rottenApples += list(cm.select('dev > 2')['subject_nr'])
#cm = dm.collapse(['subject_nr'], 'correct')
#gm = cm['mean'].mean()
#gstd = cm['mean'].std()
#cm = cm.addField('dev', dtype=float)
#cm['dev'] = np.abs(cm['mean']-gm)/gstd
#rottenApples += list(cm.select('dev > 2')['subject_nr'])
print 'Rotten apples: %s' % rottenApples
for ra in rottenApples:
	dm = dm.select('subject_nr != %d' % ra)

# The effective N
print 'N (after) = %d' % dm.count('subject_nr')
print 'Trials (after) = %d' % len(dm)

# Summarize RT (incl. incorrect) and accuracy
print "RT\t", dm['response_time'].mean(), dm['response_time'].std()
print "Accuracy\t", dm['correct'].mean(), dm['correct'].std()

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

print dm.collapse(['valid', 'cond'], 'targetRot')

dm = dm.select('practice == "no"')
#dm = dm.selectByStdDev(['subject_nr'], 'response_time', 2.5)
dm_cor = dm.select('correct == 1')

dm_cor = dm_cor.select('subject_nr >= 3000')

# Some more info
print "Correct RT\t%.2f\t%.2f" % (dm_cor['response_time'].mean(), \
	dm_cor['response_time'].std())
print "Accuracy\t%.2f\t%.2f" % (dm['correct'].mean(), dm['correct'].std())
print "SOA short\t%.2f\t%.2f" % (dm_cor.select('postRotDelay == 0')['soa'] \
	.mean(), dm_cor.select('postRotDelay == 0', verbose=False)['soa'].std())
print "SOA long\t%.2f\t%.2f" % (dm_cor.select('postRotDelay == 1000')['soa'] \
	.mean(), dm_cor.select('postRotDelay == 1000', verbose=False)['soa'].std())
print "target_dur\t%.2f\t%.2f" % (dm_cor['target_dur'].mean(), \
	dm_cor['target_dur'].std())
print "cue_dur\t%.2f\t%.2f" % (dm_cor['cue_dur'].mean(), \
	dm_cor['cue_dur'].std())
print "rot_dur\t%.2f\t%.2f" % (dm_cor['rot_dur'].mean(), \
	dm_cor['rot_dur'].std())
print

# Call the helper functions for the actual analysis!
for arg in sys.argv:
	if hasattr(helpers, arg):
		print '*** Calling helpers.%s() ***' % arg
		getattr(helpers, arg)(dm)
