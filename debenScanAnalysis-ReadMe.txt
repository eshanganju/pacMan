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
      Updated Segment now asks for method to follow for binarization. void ratio is passed from start, if not user is asked;
      Updated Segment to check for no particle contact when fixing oversegmentation - if not contact - moves to next label;
      Particle within particle - rectify by controlling the sequence in which particle are merged;

    #Priority:
      Gradation --> Relative Breakage --> Fabric

      - Gradation
          *Add original gradation as a check (line or spline curve?) - standard CSV file I think will work.
          Label 1 to 2 and label 2 to 1 contact areas are diferent - slightly.

      - Relative breakage
          *Br is spitting out smaller area values - area calculation is a little lower for the code than for the manual method - need to check
          *Which is incorrect? - Manual or Auto - Manual check uses lesser numerb of points gradation.

      - Fabric
          *Contact detemination --> Equal area projection / Rose diagram --> Fabric tensor

    #Secondary issues:
      - Parallelize the particle size estimation - NOT URGENT
      - Create log file in output folder
      - Get rid og user inputs in final version

    #Add:
      - Void ratio values at different load levels for all sands
      - Calibration values for all sands

  --------
  2020-04-10

    #Fixes:

    #Priority:
    Gradation --> Relative Breakage --> Fabric

      - Gradation
          *Add original gradation as a check (line or spline curve?) - standard CSV file I think will work.
          Label 1 to 2 and label 2 to 1 contact areas are diferent - slightly. This maybe an outcome of the dilation and particle surface

      - Relative breakage
          *Br is spitting out smaller area values - area calculation is a little lower for the code than for the manual method - need to check
          *Which is incorrect? - Manual or Auto - Manual check uses lesser number of points gradation.
          Check the CA particle size terms - are they giving right numbers? Check with a particle whose parameters are known
            write down the method.

      - Fabric
          *Contact detemination --> Equal area projection / Rose diagram --> Fabric tensor
          Is the correction for small contact implemented correctly in the measureContactNormalsSpamRW()
          Which orientation do we need to correct for in the projections?
          Rose diagrams instead of EA projections

    #Secondary issues:
      - Parallelize the particle size estimation - NOT URGENT
      - Create log file in output folder
      - Get rid og user inputs in final version

    #Add:
      - Void ratio values at different load levels for all sands
      - Calibration values for all sands

  --------
  2020-04-15

    #Fixes:
    Fixed relative breakage calculation to calculate correct ertive breakage - which agrees with origin code

    #Priority:
    Gradation --> Relative Breakage --> Fabric

      - Gradation
          *Add original gradation as a check (line or spline curve?) - standard CSV file I think will work.
          Label 1 to 2 and label 2 to 1 contact areas are diferent - slightly. This maybe an outcome of the dilation and particle surface

      - Relative breakage
          Check the CA particle size terms - are they giving right numbers? Check with a particle whose parameters are known
            write down the method.

      - Fabric
          *Contact detemination --> Equal area projection / Rose diagram --> Fabric tensor
          Is the correction for small contact implemented correctly in the measureContactNormalsSpamRW()

          Which orientation do we need to correct for in the projections?
            The x orientation is positive in the code. Z and Y orientations can be negative. This is true when looking at contact (weibicke's code)
            When plotting the orts, the Z direction is flipped.
            I am not sure why this is done, but it doesnt affect the results.
            All I need to make sure that we are plotting the same things as the ground reality.

          Vectorize the checking for small contact areas
          Rose diagrams instead of EA projections

    #Secondary issues:
      - Parallelize the particle size estimation - NOT URGENT
      - Create log file in output folder
      - Get rid og user inputs in final version
      - Split getParticleSize into sub functions

    #Add:
      - Void ratio values at different load levels for all sands
      - Calibration values for all sands

  --------
  2020-04-16

    #Fixes:
    Fixed relative breakage calculation to calculate correct ertive breakage - which agrees with origin code

    #Priority:
    Gradation --> Relative Breakage --> Fabric

      - Gradation

      - Relative breakage

      - Fabric
        Contact detemination --> Equal area projection / Rose diagram --> Fabric tensor
          Add code for plotting
            EAP plots
            2D rose diagrams

          Vectorize the checking for small contact areas
          Rose diagrams instead of EA projections

    #Secondary issues:
      - Parallelize the particle size estimation - NOT URGENT
      - Create log file in output folder
      - Get rid og user inputs in final version
      - Split getParticleSize into sub functions
      - Check CA calculation steps + write down the steps

    #Add:
      * Void ratio values at different load levels for all sands
      * Calibration values for all sands
_________________________________________________
