function [points,triangles] = loadvtk(filename)
// A routine to load points and triangles in vtk format
// (c) Maureen Clerc, December 2004
// Usage: [points,triangles] = loadvtk(filename)

// Ouput variables initialisation (not found in input variables)
points=[];
triangles=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);



fid = mopen(filename);
s = " ";
while  max(size(mtlb_findstr("POINTS",s)))==0 
 s = mgetl(fid,1); 
end;
npoints = msscanf(s,"POINTS %d");
points = mtlb_fscanf(fid,"%g %g %g\n",[3,npoints]);
points = points';
s = " ";
while max(size(mtlb_findstr("POLYGONS",s)))==0 
 s = mgetl(fid,1); 
end;
ntriangles = msscanf(s,"POLYGONS %d %d")
dim2 = ntriangles(2) ./ntriangles(1);
ntriangles = mtlb_e(ntriangles,1);
triangles = mtlb_fscanf(fid,"%g %g %g %g\n",[dim2,ntriangles]);
triangles = triangles';
if dim2==4 then
      triangles = triangles(:,2:4)+1;
else
  triangles = triangles+1;
end;
mclose(fid);
endfunction
