# Hurst-estimation-BBC-TV
Code for paper J. D. B. Nelson; C. Nafornita; A. Isar, "Semi-local scaling exponent estimation with box-penalty constraints and total-variation regularisation," in IEEE Transactions on Image Processing, doi: 10.1109/TIP.2016.2551365


If you use this code please cite:

J. D. B. Nelson; C. Nafornita; A. Isar, "Semi-local scaling exponent estimation with box-penalty constraints and total-variation regularisation," in IEEE Transactions on Image Processing, doi: 10.1109/TIP.2016.2551365,
URL: http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=7448417 

-----
Free for research purposes only.
-----

DESCRIPTION:
[1] generate fractal Brownian surface with piecewise constant varying
Hurst parameter

[2] estimate Hurst using DTCWT method with configurable combination of:
energy downsampling; semi-local weights; bound constraints on beta;
total variation smoothing

Installation:
- prerequisites assumes presence of INRIA's fraclab (in directory
  'fraclab') and Nick Kingsbury's DTCWT toolbox in directory (in
'dtcwt'), cf.: addpath('dtcwt', genpath('fraclab'))

- run test.m

LIST:
test.m
synthmbf.m
bbctv.m
constr_ols.m
penalty2.m
bbc.m

CONTACT INFORMATION:
Corina Nafornita,    
Politehnica University of Timisoara,
Romania
email corina.nafornita at upt.ro


