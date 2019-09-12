# DTM Generation (dsm2dtm)

Here you find a naive and thus simple Matlab implementation of the Digital Terrain Model (DTM) generation algorithm proposed in our papers [1,2]. From a given input Digital Surface Model (DSM) a DTM is extracted according to Algorithm 1 [2]. An optimized version can be found within the commercial software Remote Sensing Software Graz (RSG) [3]. There the version using two passes (cf. Figure 16, [2]) with only one label image is implemented.

If you use the code or any results in your work then refer to our papers [1,2].

# References

[1] Roland Perko, Hannes Raggam, Karlheinz Gutjahr, and Mathias Schardt. Advanced DTM generation from very high resolution satellite stereo images. In ISPRS Annals of the Photogrammetry, Remote Sensing and Spatial Information Sciences, volume II-3/W4, pages 165-172, Munich, Germany, March 2015.

@InProceedings{perko2015advanced,  
  title                    = {Advanced {DTM} generation from very high resolution satellite stereo images},  
  author                   = {Roland Perko and Hannes Raggam and Karlheinz Gutjahr and Mathias Schardt},  
  booktitle                = {ISPRS Annals of the Photogrammetry, Remote Sensing and Spatial Information Sciences},  
  year                     = {2015},  
  address                  = {Munich, Germany},  
  month                    = {March},  
  pages                    = {165-172},  
  volume                   = {II-3/W4},  
  doi                      = {10.5194/isprsannals-II-3-W4-165-2015},  
  timestamp                = {2015-03-26}  
}  

[2] Roland Perko, Hannes Raggam, and Peter M. Roth. Mapping with Pleiades - End-to-end workflow. Remote Sensing, Special Issue on 3D Reconstruction Based on Aerial and Satellite Imagery, 11(17):2052, September 2019.

@Article{perko2019mapping,  
  author    = {Roland Perko and Hannes Raggam and Peter M. Roth},  
  title     = {Mapping with {P}l{\'e}iades -- {E}nd-to-End Workflow},  
  journal   = {Remote Sensing},  
  year      = {2019},  
  volume    = {11},  
  number    = {17},  
  pages     = {2052},  
  issn      = {2072-4292},  
  doi       = {10.3390/rs11172052},  
  url       = {https://www.mdpi.com/2072-4292/11/17/2052},   
  month     = {September},  
  timestamp = {2019-09-01},  
}  

[3] https://www.remotesensing.at/remote-sensing-software/
