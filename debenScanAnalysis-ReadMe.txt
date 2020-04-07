This contains all the data for the calibrations for the deben tests and code notes at the end of the file
_________________________________________________
OTC
  Starting Details:
    D50 = 0.72 mm
    emin = 0.50
    emax = 0.75
    Gs = 2.65
    max D50 = 7.80 D50s (We carry out all analysis till 7D50 for REV)
    max cube = 5.65 mm
    Images are vertically flipped
    vertical useable lentha at Z = 513 = 7.5 mm
    Center voxel Y = 457 px
    Center voxel X = 507 px

  Load = 0 N
    Width: 1,000 px (11,930.80 µm)
    Height: 1,010 px (12,050.11 µm)
    Depth: 1,024 px (12,217.14 µm)
    Time steps: 1
    Voxels: 1,034,240,000
    Size: 1973 Mb  Data type: USHORT
    Volume: 1,756,426,686,733.08 µm³
    e_measured = 0.541 (2 significant figures + 1 extra)

  Load = 500 N
    Width: 1,000 px (11,930.80 µm)
    Height: 1,010 px (12,050.11 µm)
    Depth: 1,024 px (12,217.14 µm)
    Time steps: 1
    Voxels: 1,034,240,000
    Size: 1973 Mb  Data type: USHORT
    Volume: 1,756,426,686,733.08 µm³

  Load = 1500 N
    Width: 1,000 px (11,930.80 µm)
    Height: 1,010 px (12,050.11 µm)
    Depth: 1,024 px (12,217.14 µm)
    Time steps: 1
    Voxels: 1,034,240,000
    Size: 1973 Mb  Data type: USHORT
    Volume: 1,756,426,686,733.08 µm³

  Load = 4500 N
    Width: 1,004 px (11,978.52 µm)
    Height: 1,010 px (12,050.11 µm)
    Depth: 1,024 px (12,217.14 µm)
    Time steps: 1
    Voxels: 1,038,376,960
    Size: 1981 Mb  Data type: USHORT
    Volume: 1,763,452,393,480.01 µm³
_________________________________________________
OGF



_________________________________________________
2QR

_________________________________________________
Misc notes
  1. We choose a smaller cube to prevent effects of edge
  2. Cube should be big enough to accomodate the REV
  3. REV should be based on the stabilization of the GSDs

Code notes:
  --------
  2020-03-31
    The code works well till segmentation - the edges are removed and the dilation error that creeps up in the contactList function in
    spam can be corrected by applying a padding around the sample.

    There is undersegmentation of the particle in the edges
    I see some particles within particles - this may be a sequencing issue

  --------
  2020-04-06

    #Fixes:
      - Dilation error fixed by adding and removing a 2px padding around the sample
      - Updated relative breakage code - gives out values but not agreeying with manual method

    #Priority:
      Gradation --> Relative Breakage --> Fabric

      - Gradation
          Otsu overestimates particle size needs density based thresholding
          Particle within particle - rectify by controlling the sequence in which particle are merged
          Implement no contact particle correction in code - prevent getting stuck in loop for particle touching nothing
          Add original gradation as a check (line or spline curve?) - standard CSV file I think will work.

      - Relative breakage
          Br is spitting out smaller area values - area calculation is a little lower for the code than for the manual method - need to check
          Which is incorrect? - Manual or Auto - Manual check uses lesser numerb of points gradation.

      - Fabric
          Contact detemination --> Equal area projection / Rose diagram --> Fabric tensor

    #Secondary issues:
      - Parallelize the particle size estimation - NOT URGENT
      - Create log file in output folder

    #Add:
      - Void ratio values at different load levels for all sands
      - Calibration values for all sands

  --------
  2020-04-07

    #Fixes:
      Otsu overestimates particle size needs density based thresholding in deben analysis - Done

    #Priority:
      Gradation --> Relative Breakage --> Fabric

      - Gradation
          *Implement no contact particle correction in code - prevent getting stuck in loop for particle touching nothing
          *Particle within particle - rectify by controlling the sequence in which particle are merged
          *Add original gradation as a check (line or spline curve?) - standard CSV file I think will work.

      - Relative breakage
          *Br is spitting out smaller area values - area calculation is a little lower for the code than for the manual method - need to check
          *Which is incorrect? - Manual or Auto - Manual check uses lesser numerb of points gradation.

      - Fabric
          *Contact detemination --> Equal area projection / Rose diagram --> Fabric tensor

    #Secondary issues:
      - Parallelize the particle size estimation - NOT URGENT
      - Create log file in output folder

    #Add:
      - Void ratio values at different load levels for all sands
      - Calibration values for all sands
_________________________________________________
