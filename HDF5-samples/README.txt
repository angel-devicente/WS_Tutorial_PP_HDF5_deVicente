Files from in-house applications:

 + From Mancha3D code (MHD)
  + orszag_0051.h5
  + orszag_0051_bx.vtk - Bx dataset from orszag_0051.h5 in VTK format
    "h5tovtk orszag_0051.h5:bx -o orszag_0051_bx.vtk"
  + orszag_0051.B.vtk - B vector from  orszag_0051.h5 in VTK format
    "h5tovtk orszag_0051.h5:bx orszag_0051.h5:by orszag_0051.h5:bz -o orszag_0051.B.vtk"
   

 + From PORTA code (RT)
  + twolevel.h5 
  + twolevel_gz8.h5 - Same as above, but compressed with
    "h5repack -f GZIP=8 twolevel.h5 twolevel_gz8.h5"


 + From the sample Heat2D code for this tutorial
  + iac.h5 - Sample IAC logo in HDF5
  + heat_0000.h5 - First output file of Heat2D code (should be identical to
    iac.h5)
  + heat_0001.h5 - Seconde output file of Heat2D
   
 
Other formats that rely on HDF5:

 + Chombo
  Chombo/paramesh_chombo_0003.hdf5 (AMR output, generated with
   Mancha3D+PARAMESH)
  
 + CGNS
  CGNS/multi_hdf5.cgns (from http://cgns.sourceforge.net/CGNSFiles.html)
  
 + netCDF-4
  netCDF-4/OMI-Aura_L2-example.nc (from
   https://www.unidata.ucar.edu/software/netcdf/examples/files.html)
  
