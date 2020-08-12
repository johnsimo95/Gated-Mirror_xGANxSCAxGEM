# QEM_GANxGEMxSCA
## Code for the work on the gated electron mirror, including code to simulate our circuitry, the frequency response of different loads, and ray tracing deflections of different loads 

Code segments include:
- COMSOL files for RF simulations of various loads exporting
  1. S-Parameters for interacting loads with the LTSpice simulations
  2. Transfer functions for taking the port voltages and extracting relevant field quantities at points of interest
- MATLAB Files for:
  1. Generating Spice netlists
  2. Convolving LTSpice temporal outputs with the COMSOL transfer functions to get the full system response
- PSpice files for:
  1. Interacting GaNFETs with various loads
  2. Sweeping Parameters
  
- Also included in this are auxilary code segments for
  1. Estimaing parasitic and resonant quantities for various systems
  2. CAD files for different load types
  3. Generating figures and plots.
  4. Reading and plotting oscilloscope data
  
This code heavily references a generic MATLAB library we have been developing:
https://github.com/johnsimo95/QEM-MATLAB.git

  TESTING CHANGES
