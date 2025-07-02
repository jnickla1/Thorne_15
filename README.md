Code for evaluating various methods for determining a long-term temperature trend for the paper "How will we know when 1.5°C of global warming has been reached?" (preprint DOI to be added soon).

To make this code work on your machine, you need to also install a climate_data directory "somewhere" (link will be added soon). Then we need to overwrite a bunch of file paths saved throughout the code, relative to the home directory:

find . -type f -name '*.py'  -exec sed -i '' 's|data/jnickla1/climate_data|climate_data|g' '{}' \;
find . -type f -name '*.py'  -exec sed -i '' 's|data/jnickla1/Thorne_15_codefigurestats|Downloads/Thorne_15_codefigurestats|g' '{}' \;
find . -type f -name '*.m'  -exec sed -i '' 's|data/jnickla1/climate_data|climate_data|g' '{}' \;
find . -type f -name '*.m'  -exec sed -i '' 's|data/jnickla1/Thorne_15_codefigurestats|Downloads/Thorne_15_codefigurestats|g' '{}' \;
 

The current results for the historical hindcase comparisons are linked here: https://docs.google.com/spreadsheets/d/16_MBRPqwyCllO3Q3OCW7JB84l3aNQV9klEszcq0WGpU/edit?gid=271378293#gid=271378293
Progress / note sheet with more detailed citation information is here: https://docs.google.com/spreadsheets/d/1iShljXO2rmPHpPjPkBGbPSwZc_7XVtCU32O3sn6sjTc/edit?gid=0#gid=0
