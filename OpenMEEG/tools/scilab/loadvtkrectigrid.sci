function [x,y,z] = loadvtkrectigrid(filename)
// A routine to load data under a rectilinear grid format
// (c) Maureen Clerc, April 2008

// Ouput variables initialisation (not found in input variables)
x=[];
y=[];
z=[];
// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);
fid = mopen(filename);
s = " ";
while  max(size(mtlb_findstr("DIMENSIONS",s)))==0 
 s = mgetl(fid,1); 
end
[nx,ny,nz] = sscanf(s,"DIMENSIONS %d %d %d\n");
 s = " ";
 while max(size(mtlb_findstr("X_COORDINATES",s)))==0 
 s = mgetl(fid,1); 
 end
  nnx = sscanf(s,"X_COORDINATES %d float\n");
 if ~(nnx==nx) 
  disp('Error in x dimension!")
 else
  for i=1:nx
  x(i) = mfscanf(fid,"%g");
  end
 end
  s = [];
 while length(s)==0
 s = mgetl(fid,1); 
 end
 while max(size(mtlb_findstr("Y_COORDINATES",s)))==0 
 s = mgetl(fid,1); 
 end
 nny = sscanf(s,"Y_COORDINATES %d float\n");
 if ~(nny==ny) 
  disp('Error in y dimension!")
 else
  for i=1:ny
  y(i) = mfscanf(fid,"%g");
  end
 end
  s = [];
 while length(s)==0
 s = mgetl(fid,1); 
 end
 while max(size(mtlb_findstr("Z_COORDINATES",s)))==0 
 s = mgetl(fid,1); 
disp(s)
 end
 nnnz = sscanf(s,"Z_COORDINATES %d float\n");
 if ~(nnnz==nz) 
  disp('Error in z dimension!")
 else
  for i=1:nz
  z(i) = mfscanf(fid,"%g");
  end
 end
output = mclose(fid);
endfunction
