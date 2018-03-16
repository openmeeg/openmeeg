function [points,triangles] = loadtri(filename,option)
// A routine to load points and triangles in tri format
// (c) Maureen Clerc, May 2008
// Usage: [points,triangles] = loadtri(filename,option)
//option: 'full' means extract the 6 values from the point field
//          if not, only extract first 3 values
// Ouput variables initialisation (not found in input variables)
points=[];
triangles=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

// Number of arguments in function call
[%nargout,%nargin] = argn(0);

if %nargin<2
    option = 'light'; // an option to diregard the normal 
end

fid = mopen(filename);
s = " ";
while  max(size(mtlb_findstr("-",s)))==0 
 s = mgetl(fid,1);
end;
npoints = msscanf(s,"- %d");
for i=1:npoints
  points(i,:) = mfscanf(fid,"%g %g %g %g %g %g\n");
end
if option=='light'
    points = points(:,1:3);
end
s = " ";
while max(size(mtlb_findstr("-",s)))==0 
 s = mgetl(fid,1);
end;
ntriangles = msscanf(s,"- %d");
dim2 = 3;
for i=1:ntriangles
  triangles(i,:) = mfscanf(fid,"%g %g %g\n");
end
triangles = triangles+1;
mclose(fid);
endfunction
