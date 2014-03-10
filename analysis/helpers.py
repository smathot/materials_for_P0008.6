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
import numpy as np
from scipy.stats import ttest_rel, linregress, ttest_1samp
import numpy.fft.fftpack as F
from exparser import Constants, Fitting
from exparser.PivotMatrix import PivotMatrix
from exparser.AnovaMatrix import AnovaMatrix
from exparser.DataMatrix import DataMatrix
from matplotlib import pyplot as plt
from exparser.TangoPalette import *
from exparser.CsvReader import CsvReader
from scipy.optimize import curve_fit
from scipy.special import erf

N = 42

plt.rc("font", family='Arial')
plt.rc("font", size=10)

def desc(dm):

	"""
	Save descriptive statistics.

	Arguments:
	dm		--	A DataMatrix.
	"""

	dm_cor = dm.select('correct == 1')

	pm = PivotMatrix(dm, ['postRotDelay'], ['subject_nr'], \
	'soa', colsWithin=False)
	pm.save('output/soa.csv')
	pm = PivotMatrix(dm, ['subject_nr'], ['subject_nr'], \
		'correct', colsWithin=False)
	pm.save('output/accuracy.csv')
	pm = PivotMatrix(dm_cor, ['subject_nr'], ['subject_nr'], \
		'response_time', colsWithin=False)
	pm.save('output/rt.csv')
	pm = PivotMatrix(dm, ['postRotDelay', 'targetRot'], ['subject_nr'], \
		'correct', colsWithin=False)
	pm.save('output/acc.targetrot.postRotDelay.csv')
	pm = PivotMatrix(dm_cor, ['postRotDelay', 'targetRot'], ['subject_nr'], \
		'response_time', colsWithin=False)
	pm.save('output/rt.targetrot.postRotDelay.csv')
	# Do the split analysis on RTs
	pm = PivotMatrix(dm_cor, ['cond', 'postRotDelay', 'valid'], \
		['subject_nr'], dv='response_time')
	pm.save('output/rt.cond.postRotDelay.valid.pm.csv')
	# Do the full split on RTs after applying 2.5 SD
	pm = PivotMatrix(dm_cor.selectByStdDev(['subject_nr'], 'response_time'), \
		['cond', 'postRotDelay', 'valid'], ['subject_nr'], dv='response_time')
	pm.save('output/rt.2.5sd.cond.postRotDelay.valid.pm.csv')
	# Do the full split on iRts
	pm = PivotMatrix(dm_cor, ['cond', 'postRotDelay', 'valid'], \
		['subject_nr'], dv='iRt')
	pm.save('output/iRt.cond.postRotDelay.valid.pm.csv')

def lmeCorrect(dm):

	"""
	Runs the overall accuracy LME.

	Arguments:
	dm		--	A DataMatrix.
	"""

	lme(dm, dv='correct')

def lme(dm, dv='iRt'):

	"""
	Runs the overall LME.

	Arguments:
	dm		--	A DataMatrix.

	Keyword arguments:
	dv		--	The dependent variable. (default='iRt')
	"""

	assert(dv in ['iRt', 'correct'])
	if dv == 'iRt':
		dm = dm.select('correct == 1')
	from exparser.RBridge import RBridge
	R = RBridge()
	R.load(dm)

	# Run all lmers
	lm = R.lmer('%s ~ postRotDelay*valid*cond + (1|subject_nr)' % dv, \
		lmerVar='lmerFull')
	lm._print(sign=5)
	lm.save('output/lme.%s.full.csv' % dv)

	lm = R.lmer('%s ~ postRotDelay*valid + (1|subject_nr)' % dv, \
		lmerVar='lmerNoCond')
	lm._print(sign=5)
	lm.save('output/lme.%s.noCond.csv' % dv)

	lm = R.lmer('%s ~ valid*cond + (1|subject_nr)' % dv, \
		lmerVar='lmerNoDelay')
	lm._print(sign=5)
	lm.save('output/lme.%s.noDelay.csv' % dv)

	lm = R.lmer('%s ~ postRotDelay*cond + (1|subject_nr)' % dv, \
		lmerVar='lmerNoValid')
	lm._print(sign=5)
	lm.save('output/lme.%s.noValid.csv' % dv)

	# Do the model comparisons
	am = R.anova('lmerNoCond', 'lmerFull')
	am._print(sign=5)
	am.save('output/anova.%s.noCond.csv' % dv)

	am = R.anova('lmerNoDelay', 'lmerFull')
	am._print(sign=5)
	am.save('output/anova.%s.noDelay.csv' % dv)

	am = R.anova('lmerNoValid', 'lmerFull')
	am._print(sign=5)
	am.save('output/anova.%s.noValid.csv' % dv)

	R.load(dm.select('postRotDelay == 0'))
	lm = R.lmer('%s ~ valid + (1|subject_nr)' % dv)
	lm._print(sign=5)
	lm.save('output/lme.%s.0.noCond.csv' % dv)

	R.load(dm.select('postRotDelay == 1000'))
	lm = R.lmer('%s ~ valid + (1|subject_nr)' % dv)
	lm._print(sign=5)
	lm.save('output/lme.%s.1000.noCond.csv' % dv)

def lmePlotCorrect(dm):

	"""
	Plots the accuracy graph for the overall LME.

	Arguments:
	dm		--	A DataMatrix.
	"""

	lmePlot(dm, dv='correct')

def lmePlot(dm, dv='iRt'):

	"""
	Plots the graph for the overall LME.

	Arguments:
	dm		--	A DataMatrix.

	Keyword arguments:
	dv		--	The dependent variable. (default='iRt')
	"""

	assert(dv in ['iRt', 'correct'])
	if dv == 'iRt':
		dm = dm.select('correct == 1')
	from exparser.RBridge import RBridge
	R = RBridge()
	# Now plot!
	if dv == 'iRt':
		fig = plt.figure(figsize=(5,3))
	else:
		fig = plt.figure(figsize=(5,1.5))
	plt.subplots_adjust(wspace=0, bottom=.15)
	i = 1
	for postRotDelay in (0, 1000):
		_dm = dm.select('postRotDelay == %s' % postRotDelay)
		_dmObj = _dm.select('cond == "object-based"')
		_dmSpa = _dm.select('cond == "spatial"')
		_dmVal = _dm.select('valid == "{0}valid"')
		_dmInv = _dm.select('valid == "{1}invalid"')
		_dmObjVal = _dmObj.select('valid == "{0}valid"')
		_dmObjInv = _dmObj.select('valid == "{1}invalid"')
		_dmSpaVal = _dmSpa.select('valid == "{0}valid"')
		_dmSpaInv = _dmSpa.select('valid == "{1}invalid"')

		R.load(_dm)
		__dm = R.lmer('%s ~ valid + (1|subject_nr)' % dv)
		__dm._print(sign=5)
		if dv == 'iRt':
			mVal = 1 / _dmVal[dv].mean()
			mInv = 1 / _dmInv[dv].mean()
			mObjVal = 1 / _dmObjVal[dv].mean()
			mObjInv = 1 / _dmObjInv[dv].mean()
			mSpaVal = 1 / _dmSpaVal[dv].mean()
			mSpaInv = 1 / _dmSpaInv[dv].mean()
			# Get the errorbars!
			lo = ((__dm['est'][0]+__dm['ci95lo'][1])**-1 - \
				(__dm['est'][0]+__dm['est'][1])**-1) / 2
			up = ((__dm['est'][0]+__dm['ci95up'][1])**-1 - \
				(__dm['est'][0]+__dm['est'][1])**-1) / 2
		elif dv == 'correct':
			mVal = 100. - 100.*_dmVal[dv].mean()
			mInv = 100. - 100.*_dmInv[dv].mean()
			mObjVal = 100. - 100.*_dmObjVal[dv].mean()
			mObjInv = 100. - 100.*_dmObjInv[dv].mean()
			mSpaVal = 100. - 100.*_dmSpaVal[dv].mean()
			mSpaInv = 100. - 100.*_dmSpaInv[dv].mean()
			# Get the errorbars!
			lo = ((__dm['est'][0]+__dm['ci95lo'][1]) - \
				(__dm['est'][0]+__dm['est'][1])) / 2
			up = ((__dm['est'][0]+__dm['ci95up'][1]) - \
				(__dm['est'][0]+__dm['est'][1])) / 2
			lo *= 100.
			up *= 100.
		eVal = [lo, up]
		eInv = [lo, up]
		plt.subplot(1,2,i)
		plt.errorbar([0,1], [mVal, mInv], yerr=[eVal, eInv], fmt='o-', \
			label='Preferred model', color='black')
		plt.plot([0,1], [mObjVal, mObjInv], '--', label='Object-centered', \
			color='black')
		plt.plot([0,1], [mSpaVal, mSpaInv], ':', label='Retinotopic', \
			color='black')
		plt.xlim(-.2, 1.2)
		plt.xticks([0,1], ['Valid', 'Invalid'])
		if dv == 'correct':
			plt.ylim(9, 15)
			plt.yticks([10,12,14])
		else:
			plt.ylim(595, 665)
		plt.xlabel('Cue Validity')
		if i == 2:
			plt.title('Long SOA')
			plt.gca().yaxis.set_ticklabels([])
			if dv != 'correct':
				plt.legend(frameon=False)
		else:
			plt.title('Short SOA')
			if dv == 'correct':
				plt.ylabel('Error rate (%)')
			else:
				plt.ylabel('Response time (ms)')
		i += 1
	plt.savefig('plot/lme.%s.png' % dv)
	plt.savefig('plot/lme.%s.svg' % dv)
	plt.show()

def corr(dm):

	"""
	Plots the between-subjects correlation between IOR and facilitation in
	object-centered and retinotopic conditions.

	Arguments:
	dm		--	A DataMatrix.
	"""

	#dv = 'response_time'
	dv = 'iRt'
	dm = dm.select('correct == 1')

	l = [ ['ObjCue0', 'ObjCue1000', 'SpaCue0', 'SpaCue1000'] ]
	for _dm in dm.group('subject_nr'):

		_dm0 = _dm.select('postRotDelay == 0', verbose=False)
		_dm1000 = _dm.select('postRotDelay == 1000', verbose=False)

		Obj0Inv = _dm0.select('cond == "object-based"', verbose=False) \
			.select('valid == "{1}invalid"', verbose=False)[dv] \
			.mean()
		Obj0Val = _dm0.select('cond == "object-based"', verbose=False) \
			.select('valid == "{0}valid"', verbose=False)[dv] \
			.mean()
		Obj1000Inv = _dm1000.select('cond == "object-based"', verbose=False) \
			.select('valid == "{1}invalid"', verbose=False)[dv] \
			.mean()
		Obj1000Val = _dm1000.select('cond == "object-based"', verbose=False) \
			.select('valid == "{0}valid"', verbose=False)[dv] \
			.mean()
		Spa0Inv = _dm0.select('cond == "spatial"', verbose=False) \
			.select('valid == "{1}invalid"', verbose=False)[dv] \
			.mean()
		Spa0Val = _dm0.select('cond == "spatial"', verbose=False) \
			.select('valid == "{0}valid"', verbose=False)[dv] \
			.mean()
		Spa1000Inv = _dm1000.select('cond == "spatial"', verbose=False) \
			.select('valid == "{1}invalid"', verbose=False)[dv] \
			.mean()
		Spa1000Val = _dm1000.select('cond == "spatial"', verbose=False) \
			.select('valid == "{0}valid"', verbose=False)[dv] \
			.mean()

		if dv == 'iRt':
			Obj0Inv = 1./Obj0Inv
			Obj0Val = 1./Obj0Val
			Obj1000Inv = 1./Obj1000Inv
			Obj1000Val = 1./Obj1000Val
			Spa0Inv = 1./Spa0Inv
			Spa0Val = 1./Spa0Val
			Spa1000Inv = 1./Spa1000Inv
			Spa1000Val = 1./Spa1000Val

		ObjCue0 = Obj0Inv - Obj0Val
		ObjCue1000 = Obj1000Inv - Obj1000Val
		SpaCue0 = Spa0Inv - Spa0Val
		SpaCue1000 = Spa1000Inv - Spa1000Val

		l.append( [ObjCue0, ObjCue1000, SpaCue0, SpaCue1000] )
		print '%s\t%s\t%s\t%s' % (ObjCue0, ObjCue1000, SpaCue0, SpaCue1000)

	_dm = DataMatrix(l)
	_dm.save('output/corr.%s.csv' % dv)

	fig = plt.figure(figsize=(8,8))
	plt.subplots_adjust(wspace=.3, hspace=.3)

	plt.subplot(221)
	regressplot(_dm['ObjCue0'], _dm['SpaCue0'])
	plt.xlim(-200, 200)
	plt.ylim(-200, 200)
	plt.xlabel('Obj. cuing effect / short SOA (ms)')
	plt.ylabel('Ret. cuing effect / short SOA (ms)')

	plt.subplot(222)
	regressplot(_dm['ObjCue1000'], _dm['SpaCue1000'])
	plt.xlim(-200, 200)
	plt.ylim(-200, 200)
	plt.xlabel('Obj. cuing effect / long SOA (ms)')
	plt.ylabel('Ret. cuing effect / long SOA (ms)')

	plt.subplot(223)
	regressplot(_dm['ObjCue0'], _dm['ObjCue1000'])
	plt.xlim(-200, 200)
	plt.ylim(-200, 200)
	plt.xlabel('Obj. cuing effect / short SOA (ms)')
	plt.ylabel('Obj. cuing effect / long SOA (ms)')

	plt.subplot(224)
	regressplot(_dm['SpaCue0'], _dm['SpaCue1000'])
	plt.xlabel('Ret. cuing effect / short SOA (ms)')
	plt.ylabel('Ret. cuing effect / long SOA (ms)')
	plt.xlim(-200, 200)
	plt.ylim(-200, 200)

	plt.savefig('plot/corr.%s.png' % dv)
	plt.savefig('plot/corr.%s.svg' % dv)
	plt.show()

def regressplot(x, y, title=None):

	"""
	Create a regression plot.

	Arguments:
	x		--	The X data.
	y		--	The Y data.
	title	--	The figure title.
	"""

	if isinstance(x, list):
		x = np.array(x)
	if isinstance(y, list):
		y = np.array(y)
	plt.axvline(linestyle=':', color='black')
	plt.axhline(linestyle=':', color='black')
	plt.plot(x, y, '.', color='black')
	xerr = x.std() / np.sqrt(len(x))
	yerr = y.std() / np.sqrt(len(y))
	s, i, r, p, se = linregress(x, y)
	xFit = np.array( [-200, 200] )
	yFit = s * xFit + i
	plt.plot(xFit, yFit, '-', color='black')
	if title == None:
		plt.title('r = %.4f, p = %.4f' % (r, p))
	else:
		plt.title('%s [r = %.2f, p = %.2f]' % (title, r, p))
	return r

def timing(dm):

	"""
	Validates the timing and presentation durations. Note that the first two
	sessions were affected by a clock issue, affecting the logged timestamps but
	not the actual timing, and we therefore analyze only the last session.

	Arguments:
	dm		--	A DataMatrix.
	"""

	dm = dm.select('subject_nr >= 3000')
	print "SOA short\t%.2f\t%.2f" % (dm.select('postRotDelay == 0')['soa'] \
		.mean(), dm.select('postRotDelay == 0', verbose=False)['soa'].std())
	print "SOA long\t%.2f\t%.2f" % (dm.select('postRotDelay == 1000')['soa'] \
		.mean(), dm.select('postRotDelay == 1000', verbose=False)['soa'].std())
	print "target_dur\t%.2f\t%.2f" % (dm['target_dur'].mean(), \
		dm['target_dur'].std())
	print "cue_dur\t%.2f\t%.2f" % (dm['cue_dur'].mean(), \
		dm['cue_dur'].std())
	print "rot_dur\t%.2f\t%.2f" % (dm['rot_dur'].mean(), \
		dm['rot_dur'].std())

