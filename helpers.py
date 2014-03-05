import sys
import numpy as np
from scipy.stats import ttest_rel, linregress
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

def plots(dm):

	"""
	Create plots and do separate anova's for each postRotDelay.

	Arguments:
	dm		--	A DataMatrix.
	"""

	dm_cor = dm.select('correct == 1')
	# Plots
	i = 1
	fig = plt.figure(figsize=(12,12))
	for postRotDelay in (0, 1000):

		_dm_cor = dm_cor.select('postRotDelay == %d' % postRotDelay)
		_dm = dm.select('postRotDelay == %d' % postRotDelay)

		# RTs
		plt.subplot(2,2,i)
		pm = PivotMatrix(_dm_cor, ['cond', 'valid'], ['subject_nr'], \
			'response_time', colsWithin=True)
		pm.save('output/rt.cond.valid.%d.csv' % postRotDelay)
		plt.title('postRotDelay = %d' % postRotDelay)
		Constants.capSize = 5
		pm.linePlot(fig=fig, yLabel='Response time (ms)', xLabels= \
			['Valid', 'Invalid'], lLabels= \
			['Object-centered', 'Retinotopic'], xLabel='Cue')

		# AnovaMatrices
		am = AnovaMatrix(_dm_cor, ['valid', 'cond'], \
			dv='response_time', subject='subject_nr')
		am._print(title='RT')
		am.save('output/aov.rt.%d.csv' % postRotDelay)

		# Accuracy
		plt.subplot(2,2,i+2)
		pm = PivotMatrix(_dm, ['cond', 'valid'], ['subject_nr'], \
			'correct', colsWithin=True)
		pm.save('output/correct.cond.valid.%d.csv' % postRotDelay)
		plt.title('postRotDelay = %d' % postRotDelay)
		Constants.capSize = 5
		pm.linePlot(fig=fig, yLabel='Accuracy (prop.)', xLabels= \
			['Valid', 'Invalid'], lLabels= \
			['Object-centered', 'Retinotopic'], xLabel='Cue')

		# AnovaMatrices
		am = AnovaMatrix(_dm, ['valid', 'cond'], \
			dv='correct', subject='subject_nr')
		am._print(title='RT')
		am.save('output/aov.correct.%d.csv' % postRotDelay)

		i += 1

	plt.savefig('plot/results.line.png')

def aov(dm):

	"""
	Run overall anova.

	Arguments:
	dm		--	A DataMatrix.
	"""

	# AnovaMatrices
	am = AnovaMatrix(dm_cor, ['valid', 'cond', 'postRotDelay'], \
		dv='response_time', subject='subject_nr')
	am._print(title='RT')
	am.save('output/aov.rt.csv')

	am = AnovaMatrix(dm, ['valid', 'cond', 'postRotDelay'], dv='correct', \
		subject='subject_nr')
	am._print(title='Accuracy')
	am.save('output/aov.correct.csv')

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
	R.write('postRotDelay <- factor(postRotDelay)')

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

	am = R.anova('lmerNoCond', 'lmerFull')
	am._print(sign=5)
	am.save('output/anova.%s.noCond.csv' % dv)

	am = R.anova('lmerNoDelay', 'lmerFull')
	am._print(sign=5)
	am.save('output/anova.%s.noDelay.csv' % dv)

	am = R.anova('lmerNoValid', 'lmerFull')
	am._print(sign=5)
	am.save('output/anova.%s.noValid.csv' % dv)

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
	fig = plt.figure(figsize=(5,3))
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
		eVal = [lo, up]
		eInv = [lo, up]
		plt.subplot(1,2,i)
		if dv != 'correct':
			plt.errorbar([0,1], [mVal, mInv], yerr=[eVal, eInv], fmt='o-', \
				label='Preferred model', color='black')
		plt.plot([0,1], [mObjVal, mObjInv], '--', label='Object-centered', \
			color='black')
		plt.plot([0,1], [mSpaVal, mSpaInv], ':', label='Retinotopic', \
			color='black')
		plt.xlim(-.2, 1.2)
		plt.xticks([0,1], ['Valid', 'Invalid'])
		if dv == 'correct':
			plt.ylim(8, 18)
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

	dv = 'response_time'
	dm = dm.select('correct == 1')

	l = [ ['ObjCue0', 'ObjCue1000', 'SpaCue0', 'SpaCue1000'] ]
	for _dm in dm.group('subject_nr'):

		_dm0 = _dm.select('postRotDelay == 0', verbose=False)
		_dm1000 = _dm.select('postRotDelay == 1000', verbose=False)

		ObjCue0 = _dm0.select('cond == "object-based"', verbose=False) \
			.select('valid == "{1}invalid"', verbose=False)[dv] \
			.mean() - \
			_dm0.select('cond == "object-based"', verbose=False) \
			.select('valid == "{0}valid"', verbose=False)[dv] \
			.mean()
		ObjCue1000 = _dm1000.select('cond == "object-based"', verbose=False) \
			.select('valid == "{1}invalid"', verbose=False)[dv] \
			.mean() - \
			_dm0.select('cond == "object-based"', verbose=False) \
			.select('valid == "{0}valid"', verbose=False)[dv] \
			.mean()
		SpaCue0 = _dm0.select('cond == "spatial"', verbose=False) \
			.select('valid == "{1}invalid"', verbose=False)[dv] \
			.mean() - \
			_dm0.select('cond == "object-based"', verbose=False) \
			.select('valid == "{0}valid"', verbose=False)[dv] \
			.mean()
		SpaCue1000 = _dm1000.select('cond == "spatial"', verbose=False) \
			.select('valid == "{1}invalid"', verbose=False)[dv] \
			.mean() - \
			_dm0.select('cond == "object-based"', verbose=False) \
			.select('valid == "{0}valid"', verbose=False)[dv] \
			.mean()

		if dv == 'iRt':
			ObjCue0 = 1./ObjCue0
			ObjCue1000 = 1./ObjCue1000
			SpaCue0 = 1./SpaCue0
			SpaCue1000 = 1./SpaCue1000

		l.append( [ObjCue0, ObjCue1000, SpaCue0, SpaCue1000] )
		print '%s\t%s\t%s\t%s' % (ObjCue0, ObjCue1000, SpaCue0, SpaCue1000)

	_dm = DataMatrix(l)
	_dm.save('output/corr.csv')

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

	plt.savefig('plot/corr.png')
	plt.savefig('plot/corr.svg')
	plt.show()

def regressplot(x, y, title=None):

	"""
	Create a regression plot.

	Arguments:
	x		--	The X data.
	y		--	The Y data.
	title	--	The figure title.
	"""

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

def timecourse(dm, nBin=10, subject=None):

	"""
	Plots the timecourse for object-based versus retinotopic cuing and IOR.

	Arguments:
	dm		--	A DataMatrix.

	Keyword arguments:
	nBin	--	The number of bins. (default=10)
	"""

	if subject != None:
		dm = dm.select('subject_nr == %d' % subject)
	dm = dm.addField('perc', dtype=float)
	print 'Calculating percentile scores ...'
	dm = dm.calcPerc('response_time', 'perc', ['subject_nr', 'cond', 'valid', \
		'postRotDelay'], nBin=nBin)
	print 'Done!'
	i = 0
	print 'Collapsing ...'
	dm = dm.collapse(['subject_nr', 'cond', 'postRotDelay', 'valid', 'perc'], \
		'response_time')
	print 'Done'
	dm.save('output/perc.csv')
	plt.clf()
	fig = plt.figure(figsize=(12,6))
	for postRotDelay in (0, 1000):
		i += 1
		plt.subplot(1,2,i)
		if subject == None:
			plt.ylim(-50, 70)
		colors = [blue[1], orange[1]]
		for cond in ('object-based', 'spatial'):
			col = colors.pop()
			plt.title('%s' % postRotDelay)
			_dm = dm.select('cond == "%s"' % cond).select('postRotDelay == %d' \
				% postRotDelay)
			_dmVal = _dm.select('valid == "{0}valid"')
			_dmInv = _dm.select('valid == "{1}invalid"')
			y = []
			x = []
			yerr = []
			xerr = []
			for b in _dmVal.unique('perc'):
				rtVal = _dmVal.select('perc == %f' % b)['mean']
				rtInv = _dmInv.select('perc == %f' % b)['mean']
				x.append(.5 * rtVal.mean() + .5 * rtInv.mean())
				y.append(rtInv.mean()-rtVal.mean())
				xerr.append( np.std(.5 * rtVal + .5 * rtInv) / np.sqrt(N) )
				yerr.append( np.std(rtInv-rtVal) / np.sqrt(N) )
			x = np.array(x)
			y = np.array(y)
			xerr = np.array(xerr)
			yerr = np.array(yerr)
			plt.fill_between(x, y-yerr, y+yerr, color=col, alpha=.2)
			plt.plot(x, y, 'o-', label=cond, color=col)
		plt.axhline(linestyle=':', color='black')
		plt.legend(loc='lower left', frameon=False)
	plt.savefig('plot/timecourse-%s.png' % subject)

def timecourseSubject(dm, nBinAll=2, nBinSubject=2):

	"""
	Plots the timecourse for object-based versus retinotopic cuing and IOR,
	separately for each subject.

	Arguments:
	dm		--	A DataMatrix.

	Keyword arguments:
	nBin	--	The number of bins. (default=10)
	"""

	for subject in [None] + list(dm.unique('subject_nr')):
		if subject == None:
			nBin = nBinAll
		else:
			nBin = nBinSubject
		timecourse(dm, subject=subject, nBin=nBin)

def deltaPlot(dm):

	"""
	Creates a deltaplot for validly and invalidly cued conditions, separately
	for each conditon and subject.

	Arguments:
	dm		--	A DataMatrix.
	"""

	p0 = [-1./600, 1000]

	if '--parseDelta' in sys.argv:
		l = [ ['subject_nr', 'cond', 'postRotDelay', 'peak'] ]
		blacklist = []
		for subject_nr in [None] + list(dm.unique('subject_nr')):
			if subject_nr == None:
				_dm = dm
			else:
				_dm = dm.select('subject_nr == %d' % subject_nr)
			for cond in ('object-based', 'spatial'):
				__dm = _dm.select('cond == "%s"' % cond, verbose=False)
				for postRotDelay in (0, 1000):
					___dm = __dm.select('postRotDelay == %d' % postRotDelay, \
						verbose= False)
					dmVal = ___dm.select('valid == "{0}valid"', verbose=False)
					dmInv = ___dm.select('valid == "{1}invalid"', verbose=False)
					dmVal.sort('response_time')
					dmInv.sort('response_time')
					plt.subplot(211)
					plt.xlim(-0.004, 0)
					yData = np.linspace(0, 1, len(dmVal))
					xData = -1./dmVal['response_time']
					plt.plot(xData, yData, '.', color=green[1], label= \
						'Valid')
					err, pVal = Fitting.fit(xData, yData, func=Fitting.sigmoid, \
						p0=p0, color='green', plot=True)
					if err == None:
						blacklist.append(subject_nr)
						continue
					yData = np.linspace(0, 1, len(dmInv))
					xData = -1./dmInv['response_time']
					plt.plot(xData, yData, '.', color=red[1], label='Invalid')
					err, pInv = Fitting.fit(xData, yData, func=Fitting.sigmoid, \
						p0=p0, color='red', plot=True)
					if err == None:
						blacklist.append(subject_nr)
						continue

					print pInv

					N = 100
					xData = np.linspace(-1./100, -1./3000, N)
					fVal = Fitting.sigmoid(xData, *pVal)
					fInv = Fitting.sigmoid(xData, *pInv)

					# This will take the horizontal difference
					#lDiff = []
					#for y in fVal:
						## Get the first
						#iVal = np.where(fVal >= y)[0][0]
						#iInv = np.where(fInv >= y)[0]
						#if len(iInv) == 0:
							#break
						#iInv = iInv[0]
						##print iVal, iInv
						#iRtVal = xData[iVal]
						#iRtInv = xData[iInv]
						##print iRtVal, iRtInv
						#lDiff.append(iRtInv - iRtVal)
					#fDiff = np.array(lDiff)
					#fMean = xData[:len(fDiff)]

					# This will take the vertical difference
					fDiff = fVal - fInv
					fMean = xData #(fVal+fInv)/2


					plt.subplot(212)
					plt.xlim(-0.004, 0)
					plt.plot(fMean, fDiff)
					if postRotDelay == 0:
						i = np.argmax(fDiff)
					else:
						i = np.argmin(fDiff)
					plt.axvline(fMean[i])

					rt = 1./fMean[i]
					print subject_nr, cond, postRotDelay, rt

					l.append( [subject_nr, cond, postRotDelay, rt] )

					plt.savefig('plot/delta/%s.%s.%s.png' % (cond, \
						postRotDelay, subject_nr))
					#plt.show()
					plt.clf()

		dm = DataMatrix(l)
		for subject_nr in blacklist:
			dm = dm.select('subject_nr != %s' % subject_nr)
		dm.save('output/peaks.csv')
	else:
		dm = CsvReader('output/peaks.csv').dataMatrix()
	print dm.ttest(['cond', 'postRotDelay'], 'peak')

def showLogRt(dm):

	plt.subplot(311)
	plt.hist(dm['logRt'], bins=100)
	plt.subplot(312)
	plt.hist(1/dm['response_time'], bins=100)
	plt.subplot(313)
	plt.hist(dm['response_time'], bins=100)
	plt.show()
