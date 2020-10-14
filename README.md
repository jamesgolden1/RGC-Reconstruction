# RGC-Reconstruction

A simulation of the responses of RGCs to stimulation by an electrode array and reconstruction of the perceived phosphenes.
For use with the 'phosphenePaper' [branch of isetbio](https://github.com/isetbio/isetbio/tree/phosphenePaper).

The [RGC-Reconstruction repo](https://github.com/jamesgolden1/RGC-Reconstruction) is also required.

[Simulation of visual perception and learning with a retinal prosthesis](https://www.biorxiv.org/content/10.1101/206409v4), Journal of Neural Engineering

James R. Golden, Cordelia Erickson-Davis, Nicolas P. Cottaris, Nikhil Parthasarathy, Fred Rieke, David H. Brainard, Brian A. Wandell, E.J. Chichilnisky

-------

1) In dat/ unzip WNMovie.zip  
2) Create output directory output/ in main repo  
3) Go into scripts/ and run "runReconstruct.m" to generate filters and reconstructions of the white noise stimuli.  
It also generates STAs for testing.  
4) Run testRecon.m and filterAnalyze.m to get some tests with a bar and grating stimulus,  
and tests looking at the actual filters/comparing to STAs. (This is still experimental.. will be edited to be  
more complete and accurate).  

Notes: if filters are already trained you can run "reconsFromFilt.m" to generate test reconstructions, using  
input filters.  
