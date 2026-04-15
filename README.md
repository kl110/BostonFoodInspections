# Chain and Independent Food Service Business Performance According to Income and Neighborhood in Boston
## Author: Kamila Lim
### Institution: Wellesely College
### Project for Data Science Capstone Class Spring 2026

### Project
Boston Food Inspections dataset used to investigate whether there is a 
relationship between food inspection violations and income and whether or not the business is a chain. Project has 
potential implications for food safety in national chains, where practices and level of service should be standardized.
ACS income data from 2024 and census tract boundaries from 2023 are used with Food Inspection data from 2026. Boston Food
Inspections data is too large (over 300MB) and is instead linked [here](https://data.boston.gov/dataset/food-establishment-inspections).

### Repository Structure
- Census shape files: tl_2023_25_tract. Shape file (tl_2023_25_tract.shp) used in python file.
- ACS income data: ACSST5Y2024. Income csv (ACSST5Y2024.S1901-Data.csv) used in python file.
- Poster used for research findings: LimK.PosterV1.pdf
- Modelling and analyses: LimK.Project.R
- Data cleaning and visualization: LimK.Project.ipynb

### Reproduction Steps
- Update file paths in python files for census shape files, income data, and food inspection data
- Packages needed in LimK.Project.ipynb: pandas, seaborn, re, matplotlib.pyplot, requests, Beautiful Soup,
  geopandas, Point, folium, MarkerCluster, json, matplotlib.fontmanager
- Run LimK.Project.ipynb to clean and visualize data
- Packages needed in LimK.Project.R: lme4, lmerTest, MASS, tidyverse, performance, glmmTMB, ggplot2, dplyr, broom.mixed
- Run LimK.Project.R for models, model summaries, and model visualizations.
