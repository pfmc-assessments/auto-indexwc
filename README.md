The objective of this repository is to automate index standardization for the West Coast Bottom Trawl Survey (WCGBTS) using Github Actions. 

All model fitting is performed using the `indexwc` package. The script can be accessed in the `/code` folder and is [here](https://github.com/pfmc-assessments/auto-indexwc/blob/main/code/01_calculate_indices.R)

There are two workflows that run models and process data. Model fitting is done using the [Run SDMs workflow](https://github.com/pfmc-assessments/auto-indexwc/blob/main/.github/workflows/R-run-models.yml) and output is processed with the [Commit results workflow](https://github.com/pfmc-assessments/auto-indexwc/blob/main/.github/workflows/R-commit-results.yml).

Model fitting: [![R-run-models1](https://github.com/pfmc-assessments/auto-indexwc/actions/workflows/R-run-models.yml/badge.svg)](https://github.com/pfmc-assessments/auto-indexwc/actions/workflows/R-run-models.yml)


# Accessing output

There are three ways to access output:

1. Output and diagnostics from converged models can be visualized here
[https://eric-ward.shinyapps.io/index-viewer/](https://eric-ward.shinyapps.io/index-viewer/)

2. The majority of output files (csvs, figures) can be found here
[https://github.com/pfmc-assessments/auto-indexwc/tree/autogen-results/output](https://github.com/pfmc-assessments/auto-indexwc/tree/autogen-results/output)

3. The only files not included in the `output` folder are the fitted model object and data files. Because of file size, these remain in the artifcats and can be accessed individually by examining the results of the Github actions. For example,
[https://github.com/pfmc-assessments/auto-indexwc/actions/runs/27374606207/artifacts/7577427344](https://github.com/pfmc-assessments/auto-indexwc/actions/runs/27374606207/artifacts/7577427344)

# Additional notes on workflow setup:
Github actions currently timeout after 6 hours, so by dividing the total model runs across runners, the run time can be greatly reduced and all jobs will finish. The setup used is based on https://docs.github.com/en/actions/writing-workflows/choosing-what-your-workflow-does/running-variations-of-jobs-in-a-workflow. For this application, there's a little over 40 models that need running. More runners can be added -- if so, they also need to be increased here: https://github.com/pfmc-assessments/autogen-indices/blob/e90905dab568a829406d877297bbc6ca2383ff09/code/01_calculate_indices.R#L12.

The second issue with the actions is the processing of results. Because of previous conflicts with the workflows yml, the approach here is to save results to artifacts: https://docs.github.com/en/actions/writing-workflows/choosing-what-your-workflow-does/storing-and-sharing-data-from-a-workflow. The second workflow, [https://github.com/pfmc-assessments/autogen-indices/blob/main/.github/workflows/R-commit-results.yml](https://github.com/pfmc-assessments/auto-indexwc/blob/main/.github/workflows/R-commit-results.yml), processes these artifacts, extracts .csv files, and adds them to a new branch. 
