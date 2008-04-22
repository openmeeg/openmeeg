function [output] = savevtk(filename,points,triangles,field,fieldname)
// function output = savevtk(filename,points,triangles,field,fieldname)
// A routine to save points and triangles in vtk format
// (c) Maureen Clerc, December 2004


// Ouput variables initialisation (not found in input variables)
output=[];

// Number of arguments in function call
[%nargout,%nargin] = argn(0)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


fid = mtlb_fopen(filename,"w");
fprintf(fid,"# vtk DataFile Version 2.0\n");
fprintf(fid,"File "+filename+"\n");
fprintf(fid,"ASCII\n");
fprintf(fid,"DATASET POLYDATA\n");
sz = size(points);
fprintf(fid,"POINTS %g float\n",sz(1));
fprintf(fid,"%g %g %g\n ",points);


sz = size(triangles);
triangles(:,2:4) = triangles-1;
triangles(:,1) = 3*ones(sz(1),1);
fprintf(fid,"POLYGONS %g %g\n",[sz(1),sz(1)*4]);
fprintf(fid,"%g %g %g %g\n ",triangles);

if %nargin>3 then
  field = field(:);
  if max(size(mtlb_double(field)))==sz(1) then
    fprintf(fid,"CELL_DATA %g\n",sz(1));
    fprintf(fid,"SCALARS "+fieldname+" float 1\n");
    fprintf(fid,"LOOKUP_TABLE default\n");
    fprintf(fid,"%g\n ",field);
  else
    sz = size(mtlb_double(points));
    fprintf(fid,"POINT_DATA %g\n",sz(1));
    fprintf(fid,"SCALARS "+fieldname+" float 1\n");
    fprintf(fid,"LOOKUP_TABLE default\n");
    fprintf(fid,"%g\n ",field);
  end;
end;

output = mclose(fid);
endfunction
