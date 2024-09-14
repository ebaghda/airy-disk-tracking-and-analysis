# airy-disk-tracking-and-analysis
MATLAB code for fluorescent particle localization from video data and subsequent mean-squared-displacement analysis

Code requires several packages:
- the CU Boulder BioFrontiers Imaging Toolbox by Jian Tay at https://github.com/Biofrontiers-ALMC/bioformats-matlab
- the LocalGradient localization algorithm published by Kashchuk et al. at DOI: 10.1021/acsnano.2c09787. The package requires one updated localization script (detect_trj2D_eb)
- the msdanalyzer package by tinevez at https://github.com/tinevez/msdanalyzer
