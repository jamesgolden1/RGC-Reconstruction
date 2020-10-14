# RGC-Reconstruction
For use with the 'phosphenePaper' [branch of isetbio](https://github.com/isetbio/isetbio/tree/phosphenePaper).

Collaborative repo for working on reconstruction from RGC spike responses

1) In dat/ unzip WNMovie.zip  
2) Create output directory output/ in main repo  
3) Go into scripts/ and run "runReconstruct.m" to generate filters and reconstructions of the white noise stimuli.  
It also generates STAs for testing.  
4) Run testRecon.m and filterAnalyze.m to get some tests with a bar and grating stimulus,  
and tests looking at the actual filters/comparing to STAs. (This is still experimental.. will be edited to be  
more complete and accurate).  

Notes: if filters are already trained you can run "reconsFromFilt.m" to generate test reconstructions, using  
input filters.  
