# pacMan
(ParticleAnalysisCodeMostlyforAmaturesandNoobs)

This code is for static analysis of XCT. The objectives are:
- Read the XCT files
- Filter the XCT data
- Extract subregions from the XCT data
- Segment the XCT data (binarization and watershed)
- compute particle size for segmented data
- obtain particle size distribution and its parameters
- obtain relative breakage parameters
- obtain interparticle contact orientation
- compute fabric tensors
- generate pretty plots

There are 5 modules that do these things:
- Reader
- Filter
- Segment
- Measure
- Plot


## Reader
This module reads XCT file sequences and extracts subregions from the data

## Filter
This module filters the XCT data

## Segment
This module segments the data using WS and also corrects over segmentation

## Measure
This module measures the particle size parameters and interparticle contact

## Plot
This module plots the data generated from the analysis of the XCT data


