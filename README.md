shark-bay-cazymes
==============================

Analysis of CAZyme distribution 84 metagenome-assembled genomes from a Shark Bay pustular mat

Project Organization
------------

    ├── LICENSE
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering)
    |                          and a short `-` delimited description, e.g.
    │                         `01-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    └── src                <- Source code for use in this project.
        ├── __init__.py    <- Makes src a Python module
        │
        ├── config.py      <- Stores paths and variables used throughout the project
        │
        ├── data           <- Scripts to download or generate data
        │   └── make_dataset.py
        │
        ├── features       <- Scripts to turn raw data into features for modeling
        │   └── build_features.py
        │
        └── visualization  <- Scripts to create exploratory and results oriented visualizations
            └── visualize.py
--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>

## About the data
The raw data for this project is the output of two annotation programs: dbCAN and SignalP v.5. 

One step of this project was manual annotation of data tables created from the raw data. In this workflow, data in `raw` is used to create sets of data tables stored in `interim`. These interim data tables were then manually annotated before a final processing step resulting in the finished dataset stored in `processed`. 

To create the interim data used for manual annotation, `cd` to the root directory and run `make interim_data`. To create the final dataset, run `make processed_data`. 
