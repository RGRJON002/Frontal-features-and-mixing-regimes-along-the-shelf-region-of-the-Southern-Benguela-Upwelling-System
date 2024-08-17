Author: Jonathan James Rogerson
Date:   May 2024

A series of shell, MATLAB and python scripts used to generate the various visuals and figures in the manuscript: 
"Frontal features and mixing regimes along the shelf region of the Southern Benguela Upwelling System".

The MATLAB scripts call several external packages:

(1) Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 2016. True colors of oceanography: 
    Guidelines for effective and accurate colormap selection. Oceanography 29(3):9–13. <http://dx.doi.org/10.5670/oceanog.2016.66>

(2) Kristjan Onu, Florian Huhn and George Haller, "LCS Tool : A Computational Platform for Lagrangian Coherent Structures",
    <https://github.com/jeixav/LCS-Tool-Article/>.

(3) Frederic Moisy (2024). EzyFit 2.44 (<https://www.mathworks.com/matlabcentral/fileexchange/10176-ezyfit-2-44>),
    MATLAB Central File Exchange. Retrieved August 17, 2024. 

(4) CROCO_TOOLS: refer to Penven et al. (2008) or download from  <https://www.croco-ocean.org/download/>

List of Scripts:
_____________________________________________________________________________________________________________________

**Finite Time Lyapunov Exponents**:    
-    FTLE_compute.m 
-    getvarANIM.m
-    myFTLE.m
-    myLCS.m
-    plotFTLE.m

**Lagrangian dyanmics and residence times**: 
-    cname.sh
-    Ftime.m
-    FU2day.m
-    Fto3d.m
-    TRACK_FLOAT.m
-    diag_float.m
-    probmap.m

**Canny edge detection**:
-    FRONTS_daily.m
-    myCanny.py
-    pickle_to_mat_converter.py    





