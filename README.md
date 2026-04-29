# mc-photon-transport-nanoparticle-suspension
Monte Carlo simulation of photon transport in nanoparticle suspensions, including absorption, scattering, and boundary reflections. Computes volumetric heating rates in a cuvette geometry with Gaussian or circular beam injection.

## Overview
This MATLAB code implements a Monte Carlo model for photon transport in a homogeneous nanoparticle suspension. Photon packets are propagated through a finite cuvette volume and undergo absorption, scattering, and Fresnel reflection/transmission at boundaries. The simulation outputs spatially resolved absorption events, which are converted into volumetric heating rates.

## Physical model
- Photon free paths are sampled from an exponential distribution with mean free path 1/(extinction coefficient)
- Absorption and scattering are treated as competing stochastic processes
- Scattering angles follow the Henyey–Greenstein phase function
- Boundary interactions are handled using Fresnel coefficients
- Total internal reflection is included at the liquid–air interface

## Geometry
- Rectangular cuvette defined by (0,x_length), (0,y_length), (0,z_length)
- Uniform Cartesian cubic voxel grid
- Output stored as heating rate q_{abs}(x,y,z) in mW/cm^3

## Beam definition
Two beam profiles are supported:
- Gaussian beam (defined by 1/e² radius)
- Circular (top-hat) beam
The beam propagates along the +x direction.

## Output graph
- 2D slices of heating rate q_{abs} (mW/cm³)
- Histogram of scattering events before absorption

## Usage
Run the script in MATLAB. If input variables are not predefined, a dialog will prompt for:
- Optical parameters (cross-sections, concentration, anisotropy)
- Geometry and grid resolution
- Beam parameters
- Input power

## Reference
The implementation is inspired by the following publications:
- https://doi.org/10.1016/0169-2607(95)01640-f
- https://doi.org/10.1021/nl5016975
- https://doi.org/10.1021/acsnano.8b03929
