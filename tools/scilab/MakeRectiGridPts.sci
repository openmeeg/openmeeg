function [x,y,z]=MakeRectiGridPts(xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz,posdir,ptsname)
// (c) Maureen Clerc, April 2008
ptsfile=posdir+ptsname;
vtkfile=ptsfile+'.vtk';
ptsdxfile=ptsfile+'dx';
ptsdyfile=ptsfile+'dy';
ptsdzfile=ptsfile+'dz';

x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;
nx=length(x);
ny=length(y);
nz=length(z);
if nz==1
  xx = x'*ones(1,ny);
  yy = ones(nx,1)*y;
  zz = zmin*ones(nx,ny);
  eps = min(dx,dy)/10; // value of eps for finite difference current computation
elseif nx ==1
  yy = y'*ones(1,nz);
  zz = ones(ny,1)*z;
  xx = xmin*ones(ny,nz);
  eps = min(dy,dz)/10; // value of eps for finite difference current computation
elseif ny ==1
  zz = z'*ones(1,nx);
  xx = ones(nz,1)*x;
  yy = ymin*ones(nz,nx);
  eps=min(dx,dz)/10; // value of eps for finite difference current computation
else disp('Error: one of the three dimensions must be equal to 1.')
end
savevtkrectigrid(vtkfile,x,y,z);
// dimensions of the grid np1 * np2
[np1,np2] =  size(zz);
np = np1*np2
pt = zeros(np,3);
pt(:,1) = Matrix(xx,[np 1]);
pt(:,2) = Matrix(yy,[np 1]);
pt(:,3) = Matrix(zz,[np 1]);
fprintfMat(ptsfile,pt);
newpt = pt;
newpt(:,1) = newpt(:,1)+eps
fprintfMat(ptsdxfile,newpt);
newpt = pt;
newpt(:,2) = newpt(:,2)+eps;
fprintfMat(ptsdyfile,newpt);
newpt = pt;
newpt(:,3) = newpt(:,3)+eps;
fprintfMat(ptsdzfile,newpt);

