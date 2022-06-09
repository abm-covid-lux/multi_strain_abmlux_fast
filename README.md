# Fast Multi-Strain ABMLUX

This is a version of the [multi-strain ABMlux](https://github.com/abm-covid-lux/multi_strain_abmlux) stochastic agent-based model of COVID-19, which is in turn a version of the original [ABMlux](https://github.com/abm-covid-lux/abmlux) model.

![ABMLUX Logo](abmlux_logo.jpg)

## Overview
Acting on a config, the multi-strain ABMlux model can be used to generate state files, representing a world of agents and locations. In the case of the provided config, the world represents Luxembourg, exactly as in the original ABMlux model and as used in the publication Thompson, J. and Wattam, S. "Estimating the impact of interventions against COVID-19: from lockdown to vaccination", 2021, PLOS ONE, https://doi.org/10.1371/journal.pone.0261330.

The file [abm_sim_extractor.py](https://github.com/abm-covid-lux/multi_strain_abmlux_fast/blob/main/abm_sim_extractor.py) extracts data from these state files and converts them into a new format. Taking this formatted data as input, the C++ code contained in the file [multi_strain_abm.cpp](https://github.com/abm-covid-lux/multi_strain_abmlux_fast/blob/main/multi_strain_abm.cpp) is then able to run a basic version of the multi strain model featuring the baseline and lockdown scenarios. The performance improvements reduce the runtime of each simulation to only a few minutes, allowing for high resolution parameter sweeps, for example of the strain parameters.

## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

Read the full text for details, but basically this means:
 * No commercial exploitation ([contact us](https://www.ms_abmlux.org) for another license in this case);
 * You must re-publish the source if you modify the application.

We would like this work to be useful to non-profit and academic users without significant effort.  If the license is an impediment to you using the work, please get in touch with us to discuss other licensing options.
