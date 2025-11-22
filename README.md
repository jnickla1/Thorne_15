Code for evaluating various methods for determining a long-term temperature trend for the paper "How will we know when 1.5°C of global warming has been reached?" (preprint DOI to be added soon).



## installation

### submodules
Before you do anything else, you must install a bunch of submodules git repositories that were originally developed by other authors, several of which I forked. The ones that are not forked are not used, but I investigated trying to turn them into potential methods (and possibly you will succeed in adding them to the list).
```
git submodule update --init --recursive
```

### requirements
- `anaconda` for `python3` (note do not use mamba - I know conda is slower but mamba produces a weird error.)
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
> install.packages("this.path","mgcv","KCC","abind")
```
### external climate data (compressed to just necessary files)
To make this code work on your machine, you need to also install a climate_data directory [also hosted on github](https://github.com/jnickla1/climate_data). I've already compressed down lots of netcdf files into only the averages that we need (mostly using [NCO](https://nco.sourceforge.net/)  - I've left the code in this external repo. But I didn't include copies of the original big files from the [ESGF Grid](https://esgf.github.io/nodes.html)  - you can download them yourself if you want.


### specifying the file paths
Then we need to specify file paths saved throughout the code. To do this, you MUST create a new `config.py` file that you will save at the top level of the codebase directory (that's this Thorne_15 github directory or whatever you want to rename it to). This completes the installation.

Let's say that you install the climate_data here: `~/climate_data`, and the codebase in your Downloads folder `~/Downloads/Throne_15`. Python will automatically expand out the home directory part as follows with 'os.path'. This is all that needs to be in config.py:
```
# config.py
#all are relative to home already so ~/Downloads is just Downloads
import os
CODEBASE_PATH = os.path.expanduser('~/')+"Downloads/Thorne_15"
CLIMATE_DATA_PATH = os.path.expanduser('~/')+"climate_data"
 ```
 
 Remember, this config.py file is saved here:
 ```
 Thorne_15/
│
├── config.py  # Contains the paths
├── .gitignore  # Ignores the config.py file so pulls dont break other installs
│...
├── hist_evaluation_script.py
│...
├── Methods/
├── Results/
 ```
 ##Running

Once you have that installed, the main files are run from the code root directory like this. The hist_evaluation_script.py has no elim because we are evaluating it on all methods and then eliminating none. A negative number at the very end indicates to only run that ensemble member once, positive means to run a batch of 2 (or 10 if you change the code) starting with that number. futplotcomb generates only one combination figure. For futplotcomb ESM1-2-LR the output depends on the inputted ensemble member: SSP370 and SSP126 plotted with that member highlighted, but futplotcomb for NorESM is not affected by anything - it resets both the secenario and ensemble member parameters.

```
### Run Order (so everything is created correctly), assume you have a slurm parallel computing cluster.
```
python3 hist_fitprob.py -historical
python3 hist_evaluation_script.py
sbatch run_parBATCH_future.sh
python3 plot_fut_results.py

python3 fut_evaluation_script_elim.py futplotcomb_ESM1-2-LR_SSP126_constVolc -11
python3 fut_evaluation_script_elim.py futplotcomb_NorESM_RCP45_Volc -11

sbatch run_histens.sh
python3 plot_heads_tails_probs.py
python3 hist_evaluation_script_rates.py

sbatch run_parBATCH_hist_fitprob_ff.sh
sbatch run_LOO2.sh
python3 create_combination_method_LOOffcombine.py -2
sbatch run_LOO16.sh
python3 create_combination_method_LOOffcombine.py -16
sbatch run_LOO8.sh
python3 create_combination_method_LOOffcombine.py -8

 ```
Other important figure-generating scripts are plot_fut_results.py and are found in /Illustrative_Figures.

Make sure hist_regen=True and all of the Methods subfolders are selected in the hist_evaluation_script.py (apologies for not keeping this perfectly tidy).

The current results for the [historical hindcase](https://docs.google.com/spreadsheets/d/10izz9VruI9L1pNT3pwKLlNPVhzrvGRdYk3VxvdQ1es8/edit?usp=sharing), [MPI-ESM-1-2-LR](https://docs.google.com/spreadsheets/d/1eWAeL1HHHSqyL1YF2IYaQwgMQuh8Y1RJ/edit?usp=sharing&ouid=101500668294780806861&rtpof=true&sd=true), and [NorESM1](https://docs.google.com/spreadsheets/d/1gHNtpZ4MVIw_NYp62kjtuCHZ2kZcWYwo/edit?usp=sharing&ouid=101500668294780806861&rtpof=true&sd=true) comparisons are linked

Progress / note sheet with more [detailed citation information](https://docs.google.com/spreadsheets/d/1iShljXO2rmPHpPjPkBGbPSwZc_7XVtCU32O3sn6sjTc/edit?gid=0#gid=0)
