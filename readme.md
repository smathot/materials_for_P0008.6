# P0008.6

Experimental materials that accompany the following manuscript:

[Theeuwes, J.](http://ems.psy.vu.nl/userpages/theeuwes/), [Mathôt, S](http://www.cogsci.nl/smathot)., & [Grainger, J](http://gsite.univ-provence.fr/gsite/document.php?pagendx=2044&project=lpc). (in preparation). *Object-centered Orienting and IOR*.

## Experiment

The experimental script can be found at

- `resources/0035F1.opensesame.tar.gz`

Dependencies:

- [OpenSesame](http://osdoc.cogsci.nl/)

## Analysis

The analysis scripts can be found at

- `analysis/*.py`
- `analyze.py`

To run the analysis, start `analyze.py` and specify one or more analyses to conduct. To execute the full analysis chain, run:

	python analyze lme # Perform RT analysis
	python analyze lmeCorrect # Perform error-rate analysis
	python analyze lmePlot # Create plot for RT analysis
	python analyze lmePlotCorrect # Create plot for error-rate analysis
	python analyze corr # Correlation analysis
	python analyze timing # Verify stimulus timing

A more detailed description of what each of these analyses does can be found in the docstrings of `analysis/helpers.py`.

Output will be saved in the folders `output` and `plot`, which need to be created beforehand.

Dependencies:

- [exparser](https://github.com/smathot/exparser)
- For more dependencies, see the `exparser` documentation.

## Data

Participant data can be found in plain-text `.csv` format in the folders `data/1`, `data/2`, and `data/3`.

The file `data/3/subject-8.csv` was corrupted. The last few fields from the final trial were missing, and this was corrected by copy-pasting the missing fields from the second-to-last trial. The corrupted fields did not contain relevant (for the purpose of the analysis) data.

## License

Participant data is available under the [Creative Commons Attribution-ShareAlike 4.0 International][CC-by-SA] license. This license applies to the following files:

	data/*

All other files are available under the [GNU General Public License 3][gpl].

[CC-by-SA]: http://creativecommons.org/licenses/by-sa/4.0/
[gpl]: https://www.gnu.org/copyleft/gpl.html
