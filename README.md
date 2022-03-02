# Fast-CRS
Fast estimation of CRS parameters using structure tensors. 


A MATLAB implementation of the method described in:

Fast and robust common-reflection-surface parameter estimation.<br>
Waldeland, Anders Ueland; Zhao, Hao; Faccipieri, J. H.; Solberg, Anne H Schistad; Gelius, Leiv-J.<br>
Geophysics. 2017, 83 (1), O1-O13, DOI: http://dx.doi.org/10.1190/geo2017-0113.1

Fast estimation of prestack common reflection surface parameters.<br>
Waldeland, Anders Ueland;  Tiago A. Coimbra;  Jorge H. Faccipieri;  Anne H. Schistad Solberg;  Leiv‐J. Gelius<br>
Geophysical Prospecting, 2019, online pre-paper, https://onlinelibrary.wiley.com/doi/abs/10.1111/1365-2478.12740

### Setup:
* Clone this repository
* Enter the utils/lininterpf/ directory and run build.m to compile the lininterp1f function
* Download SegyMAT and add it to the MATLAB-path (http://segymat.sourceforge.net)
* Download example data extract content into the Fast-CRS folder (https://www.dropbox.com/sh/0mi2j0ewlbmkxti/AACdLT3xKZkSDZ6nE8lGnORXa?dl=0)

### Main files:
* example_ZO_CRS.m - Example of how to run ZO CRS full search or the fast structure tensor method
* example_FO_CRS.m - Example of how to run FO CRS full search or the fast structure tensor method
* The folder nD_structure_tensor contains an n-dimensional implementation of the gradient structure tensor (GST)  and the quadratic gradient structure tensor (QST)
* The folder 2D_structure_tensor contains a 2-dimensional implementation of the gradient structure tensor which should be more easy to understand than the n-dimensional.
* structure_tensors.pdf - Details about the GST and QST.
* dip_estimation_3d.py - an (incomplete) example of how dip estimation can be done in python using the structure-tensor package (https://pypi.org/project/structure-tensor/)

# Credit:
Fast linear interpolation method by Umberto Picchini is included in this repositorty (https://se.mathworks.com/matlabcentral/fileexchange/8627-fast-linear-interpolation)

# Disclaimer:
This code was not used to generate the results in the referenced paper and has not been thoroughly tested. It is therefore not nessicarily representative for running times and results presented in the paper. The purpose of the code is rather to provide an educational presentation of the proposed method. 

# Contact:
Email: anders.u.waldeland@gmail.com
