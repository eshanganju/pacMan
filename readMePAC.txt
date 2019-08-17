Particle analysis code
  Simeple but [accurate] analysis of particles in a static CT scan
  Extracion of particle size and contact normal distribution

Approach
  Baseline data is taken from DEM example of Kalisphera
  The data has the centres of the circles and the radii of each circle. 
  This gives us the gradation and the contact norm distribution.
  
  The input is first run through Kalisphera code to get edge blur (PVE)
  Then a controlled gaussian blur is added to the data (maybe overkill?)
  Then a controlled random noise is added to the data
  These are done to simulate real CT data features
  
  The main things we need to do are: 
    1. Filtering of the CT data and check fidelity
    2. Segmentation of the CT data and check fidelity
    3. Compute measures of particle size and morphology from CT data
    4. Identify and compute contact normals from CT data

Structure
  Main code: PAC
  PAC supporting classes
    1. Particle class       (Particle.py)
    2. Aggregate class      (Aggregate.py)
    3. Filter class         (Filter.py)
    4. Segment class        (Segment.py)
    5. Measure class        (Measure.py)
    6. Visualization class  (LemmeC.py)
 
