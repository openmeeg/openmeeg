function tritovtk(trifilename,vtkfilename)
// A routine to convert tri to vtk format
// (c) Maureen Clerc, May 2008
// Usage: tritovtk(trifilename[,vtkfilename])
//option: if vtkfilename is not specified, simply replaces the tri expansion
// by a vtk extension

// Number of arguments in function call
[%nargout,%nargin] = argn(0);

if %nargin==1
  vtkfilename = strsubst(trifilename,'.tri','.vtk');
end
[pt,tr] = loadtri(trifilename);
savevtk(vtkfilename,pt,tr)
