# Bipolar_fat_water_separation

This repository contains code to perform fat-water separation for data acquired with bipolar readout gradient pulses.

Refer to the publication for more details: A flexible approach for fat-water separation with bipolar readouts and correction of gradient-induced phase and amplitude effects.

## Before starting
This code requires the Optimization Toolbox, which is available on MATLAB 2026a onwards.

This code uses the crameri perceptually uniform scientific colormaps to display some images. Add this toolbox to your matlab directory if needed https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

This code uses tools from the ISMRM fat-water toolbox. For this code to work, make sure to include the toolbox in matlab's directory.

(Optional) This code uses a phase unwrapping (used only to display some results without phase wraps). Performing phase unwrapping is optional, but if used, include in the code directory the function 'qualityGuidedUnwrapping' from Fortier and Levesque, DOI 10.1002/mrm.26989, 2017, https://gitlab.com/veronique_fortier/Quality_guided_unwrapping

## Code usage

This code presents an example for performing fat-water separation with unipolar and bipolar readout gradient pulses using a graph-cut fat-water separation approach based on Hernando et al. technique doi: 10.1002/mrm.22177. For datasets acquired with bipolar gradients, the graph-cut approach was modified to eliminate phase and amplitude effects induced by the gradient

Run 'Example_fat_water_separation_mat_files_phantom.m' to perform fat-water separation in phantom data

## Data availability 

MRI data (.mat files) for phantom experiments is available at:  hHps://osf.io/bavk7/ 

## License
Copyright Jorge Campos 2025 - MIT License

