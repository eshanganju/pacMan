# pacMan

(ParticleAnalysisCodeMadebyanAmatureNoob)

This code is for the analysis of tomography data.

The code can read, segment, and analyze tomography data from a multitude of
sources.

No oops just simple functional programming - might update in future.

The main functionalities include:
- Reading and cropping XCT data
- Filtering the XCT data
- Segmenting the XCT data (binarization and watershed)
- Computing feature size from segmented data
- Computing feature size distribution
- Computing feature contact orientations
- Computing fabric tensors

There are 5 modules that do these things:
- Read
- Filter
- Segment
- Measure
- Plot

## Read
This module reads tomography data and extracts cubical subregions for analysis.

## Filter
This module iteratively filters the tomography data using non-local means filter.

## Segment
This module segments the data using a watershed algo and corrects over segmentation.

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
- numba
- tomopy

# Speeding up the code
- Numba for numpy manipulations (future)


# Future updates
