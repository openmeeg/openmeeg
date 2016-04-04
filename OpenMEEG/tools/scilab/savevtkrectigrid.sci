function [output] = savevtkrectigrid(filename,x,y,z,field,fieldname)
// A routine to save data under a rectilinear grid format
// (c) Maureen Clerc, April 2008

// Ouput variables initialisation (not found in input variables)
output=[];

// Number of arguments in function call
[%nargout,%nargin] = argn(0)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

x = x(:)';
y = y(:)';
z = z(:)';

fid = mtlb_fopen(filename,"w");
fprintf(fid,"# vtk DataFile Version 2.0\n");
fprintf(fid,"File "+filename+"\n");
fprintf(fid,"ASCII\n");
fprintf(fid,"DATASET RECTILINEAR_GRID\n");
nx = length(x); ny = length(y); nz = length(z);
fprintf(fid,"DIMENSIONS %g %g %g\n",nx,ny,nz);
fprintf(fid,"X_COORDINATES %g float\n",nx);
mputl(string(x),fid);
fprintf(fid,"Y_COORDINATES %g float\n",ny);
mputl(string(y),fid);
fprintf(fid,"Z_COORDINATES %g float\n",nz);
mputl(string(z),fid);
if argn(2)>4
  sz = size(field);
  if min(sz)==1
  // scalar values
  if max(sz)== nx*ny*nz
    fprintf(fid,"POINT_DATA %g\n",nx*ny*nz);
    fprintf(fid,"SCALARS "+fieldname+" float\n");
    fprintf(fid,"LOOKUP_TABLE default\n");
    mputl(string(field),fid);
  else disp('Only POINT_DATA supported');
  end
  elseif min(sz)==3
  // vector values
  if max(sz)== nx*ny*nz
    fprintf(fid,"POINT_DATA %g\n",nx*ny*nz);
    fprintf(fid,"VECTORS "+fieldname+" float\n");
    fprintf(fid,"%g %g %g\n ",field);
  else disp('Only POINT_DATA supported');
  end
  end
end




output = mclose(fid);
endfunction
