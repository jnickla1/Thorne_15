Code for evaluating various methods for determining a long-term temperature trend for the paper "How will we know when 1.5°C of global warming has been reached?" (preprint DOI to be added soon).



## installation

### requirements
- `anaconda` for `python3`
- `python>=3.11` (or higher depending on Xarray requirements)
- `R>=4.4.0` and `cmake>=3.2`
- `Matlab >= R2016a` to access over the command line, just for Bayes_Sequential_Change_Point (intermediate output already saved in this repo)

First, create the `conda` environment for `python`. You will require an installation of `anaconda`. Other folks recommend [miniconda](https://docs.anaconda.com/miniconda/) or [miniforge](https://github.com/conda-forge/miniforge). After installing `conda`, run the following commands
 
```
conda env create -f cleanpy_environment.yml
conda activate cleanpy
```

As of around June 2024, some `R` environments may no longer install properly from `conda`, so we'll do this manually.

```
R
> install.packages("mgcv")
```

To make this code work on your machine, you need to also install a climate_data directory hosted on [Harvard Dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/YN5GHP). I've already compressed down lots of netcdf files into only the averages that we need (mostly using [NCO](https://nco.sourceforge.net/)  - I've left the code in this external repo as well as lots of smaller files. But I didn't include copies of the original big files from the [ESGF Grid](https://esgf.github.io/nodes.html)  - you can download them yourself if you want.

Then we need to overwrite a bunch of file paths saved throughout the code, relative to the home directory. Lots of these files may not matter, but I'm replacing all the paths regardless. Remove the extra '' below if you're not on a Mac, (-exec sed -i 's|data/), and if you are on a Mac, you might need to add this to the beginnning of each command: LC_CTYPE=C find ...)

```
find . -type f -name '*.py'  -exec sed -i '' 's|data/jnickla1/climate_data|climate_data|g' '{}' \;
find . -type f -name '*.py'  -exec sed -i '' 's|data/jnickla1/Thorne_15_codefigurestats|Downloads/Thorne_15_codefigurestats|g' '{}' \;
find . -type f -name '*.m'  -exec sed -i '' 's|data/jnickla1/climate_data|climate_data|g' '{}' \;
find . -type f -name '*.m'  -exec sed -i '' 's|data/jnickla1/Thorne_15_codefigurestats|Downloads/Thorne_15_codefigurestats|g' '{}' \;
find . -type f -name '*.pyc'  -exec sed -i '' 's|data/jnickla1/climate_data|climate_data|g' '{}' \;
find . -type f -name '*.pyc'  -exec sed -i '' 's|data/jnickla1/Thorne_15_codefigurestats|Downloads/Thorne_15_codefigurestats|g' '{}' \;
find . -type f -name '*.ipynb'  -exec sed -i '' 's|data/jnickla1/climate_data|climate_data|g' '{}' \;
find . -type f -name '*.ipynb'  -exec sed -i '' 's|data/jnickla1/Thorne_15_codefigurestats|Downloads/Thorne_15_codefigurestats|g' '{}' \;
 ```

Once you have that installed, the main files are run from the code root directory like this. Negative indicates to only run that ensemble member once, positive means to run a batch of 2 (or 10 if you change the code) starting with that number. futplotcomb generates only one combination figure. For futplotcomb ESM1-2-LR the output depends on the inputted ensemble member: SSP370 and SSP126 plotted with that member highlighted, but futplotcomb for NorESM is not affected by anything - it resets both the secenario and ensemble member parameters.

```
python3 hist_evaluation_script_elim.py
python3 fut_evaluation_script_elim.py fut_ESM1-2-LR_SSP370_constVolc -9
python3 fut_evaluation_script.py fut_ESM1-2-LR_SSP245_constVolc 9
python3 fut_evaluation_script_elim.py futplotcomb_ESM1-2-LR_SSP126_constVolc -15
python3 fut_evaluation_script_elim.py futplotcomb_NorESM_RCP45_Volc 20
python3 fut_evaluation_script_elim.py fut_NorESM_RCP45_VolcConst -28
 ```
Other important figure-generating scripts are plot_fut_results.py and are found in /Illustrative_Figures. The batch script for running things in parallel on a large computing cluster is run_parBATCH.sh (and other versions of it like run_parBATCH2.sh).

The current results for the [historical hindcase](https://docs.google.com/spreadsheets/d/10izz9VruI9L1pNT3pwKLlNPVhzrvGRdYk3VxvdQ1es8/edit?usp=sharing), [MPI-ESM-1-2-LR](https://docs.google.com/spreadsheets/d/1eWAeL1HHHSqyL1YF2IYaQwgMQuh8Y1RJ/edit?usp=sharing&ouid=101500668294780806861&rtpof=true&sd=true), and [NorESM1](https://docs.google.com/spreadsheets/d/1gHNtpZ4MVIw_NYp62kjtuCHZ2kZcWYwo/edit?usp=sharing&ouid=101500668294780806861&rtpof=true&sd=true) comparisons are linked

Progress / note sheet with more [detailed citation information](https://docs.google.com/spreadsheets/d/1iShljXO2rmPHpPjPkBGbPSwZc_7XVtCU32O3sn6sjTc/edit?gid=0#gid=0)
