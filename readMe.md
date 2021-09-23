# pacMan
(ParticleAnalysisCodeMadebyanAmatureNoob)

This code is for the analysis of XCT data.
The main functionalities include:
- Reading and cropping XCT data
- Filtering the XCT data
- Segmenting the XCT data (binarization and watershed)
- Computing particle size from segmented data
- Computing particle size distribution
- Computing relative breakage parameters
- Computing interparticle contact orientations
- Computing fabric tensors
- Generating pretty plots

There are 5 modules that do these things:
- Reader
- Filter
- Segment
- Measure
- Plot

## Reader
This module reads XCT file data and extracts cubical subregions for analysis.

## Filter
This module iteratively filters the XCT data using non-local means filter.

## Segment
This module segments the XCT data using a watershed algo and corrects over segmentation.

## Measure
This module measures the particle size parameters and interparticle contacts.

## Plot
This module plots the data generated from the analysis of the XCT data

### Requirements
pacMan uses the following libraries:
- numpy
- scipy
- statistics
- scikit-image
- tifffile
- uncertainties
- matplotlib
- spam
- glob

# Speeding up the code
- Numba for numpy manipulations (future)


