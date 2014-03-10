# P0008.6

Experimental materials that accompany:

Theeuwes, J., Mathôt, S., & Grainger, J. (in preparation). *Object-centered Orienting and IOR*.

# Experiment

The experimental script can be found at

- `resources/0035F1.opensesame.tar.gz`

Dependencies:

- [OpenSesame](http://osdoc.cogsci.nl/)

# Analysis

The analysis scripts can be found at

- `analysis/*.py`
- `analyze.py`

To run the analysis, start `analyze.py` and specify one or more analyses to conduct. To execute the full analysis chain, run:

	python analyze lme lmeCorrect lmePlot lmePlotCorrect corr timing

A description of what each of these analyses do can be found in the docstrings of `analysis/helpers.py`.

Dependencies:

- [`exparser`](https://github.com/smathot/exparser)
- For more dependencies, see the `exparser` documentation.

# Data

## Corruption

`data/3/subject-8.csv` from the third run was corrupted. The last few fields from the final trial were missing, and this was corrected by copy-pasting the missing fields from the second-to-last trial.

## License

Participant data is available under the [Creative Commons Attribution-ShareAlike 4.0 International][CC-by-SA] license. This license applies to the following files:

	data/*

All other files are available under the [GNU General Public License 3][gpl].

[CC-by-SA]: http://creativecommons.org/licenses/by-sa/4.0/
[gpl]: https://www.gnu.org/copyleft/gpl.html
