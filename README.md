ratiometric-fret
================
This github.org repository contains all code necessary for the ratiometric analysis of FRET efficiency based on the measurement of fluorescence intensities from fluorescent protein based FRET standards. The example files are included and can be exchanged for a complete set of new ones.

Voltage-2_PC-DNA.fcs :  Control cells transfected with empty plasmid for background subtraction.
Voltage-2_gfp.fcs    :  Control cells transfected with FRET donor only (GFP) for bleed-through correction.
Voltage-2_rfp.fcs    :  Control cells transfected with FRET acceptor only (RFP) for bleed-through correction.
Voltage-2_7aa.fcs    :  Cells transfected with FRET standard featuring 7 aminoacid linker.
Voltage-2_19aa.fcs   :  Cells transfected with FRET standard featuring 7 aminoacid linker.
Voltage-2_32aa.fcs   :  Cells transfected with FRET standard featuring 7 aminoacid linker.

A minimum of two different FRET standard samples are required.

Source code:

FRETeff.m       :  Matlab code for FRET efficiency analysis
dscatter.m      :  Support Matlab code for plotting 2D scatter plots (Robert Henson)
fca_readfcs.m   :  Support Matlab code for importing FCS files (Laszlo Balkay)
num2eng.m       :  Support Matlab code for converting number into engineering style
FRETresults.tex :  TEX source for creating report file


Acknowledgements:
This repository contains code from 
* FCS data reader, by Laszlo Balkay
    (https://www.mathworks.com/matlabcentral/fileexchange/9608-fcs-data-reader)
* Flow Cytometry Data Reader and Visualization, by Robert Henson
    (https://www.mathworks.com/matlabcentral/fileexchange/8430-flow-cytometry-data-reader-and-visualization)

Distributed under the BSD license
Jakub Nedbal, 2014
King's College London


