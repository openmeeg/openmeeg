function [pot,cur] = plotVolumePotential(model,cutplane,inj,extr)
// (c) Maureen Clerc, April 2008
// example: 
//   model = 'nerve0'
// inj = 1; extr = 2;
// cutplane = 'fichier_coupy';
getf('savevtk.sci')
getf('loadvtk.sci')
compdir ='../../data/Computations/';
injcode = ['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L'];
extrcode = ['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l'];
vtkpointfile = strcat([compdir model '/positions/' cutplane '.vtk']);
potfile = strcat([compdir model '/results/' model '.' cutplane '.pot' injcode(inj) extrcode(extr) '.vtk']);
curfile = strcat([compdir model '/results/' model '.' cutplane '.cur' injcode(inj) extrcode(extr) '.vtk']);
gain = fscanfMat(strcat([compdir model '/results/' model '.' cutplane '.gain']));
gaindx = fscanfMat(strcat([compdir model '/results/' model '.' cutplane 'dx.gain']));
gaindy = fscanfMat(strcat([compdir model '/results/' model '.' cutplane 'dy.gain']));
gaindz = fscanfMat(strcat([compdir model '/results/' model '.' cutplane 'dz.gain']));

injectedcurrent = zeros(12,1);
for i=1:length(inj)
  injectedcurrent(inj(i)) = 1/length(inj);
end
for i=1:length(extr)
  injectedcurrent(extr(i)) = -1/length(extr);
end
pot = gain*injectedcurrent;
//[pt,quad]=loadvtkquad(vtkpointfile);
//savevtkquad(resultsfile,pt,quad,pot,'pot');
[x,y,z]=loadvtkrectigrid(vtkpointfile);
savevtkrectigrid(potfile,x,y,z,pot,'pot');
if argn(1)>1
  // compute the current
  potdx = gaindx*injectedcurrent;
  potdy = gaindy*injectedcurrent;
  potdz = gaindz*injectedcurrent;
  cur(:,1) = potdx-pot;
  cur(:,2) = potdy-pot;
  cur(:,3) = potdz-pot;
  savevtkrectigrid(curfile,x,y,z,cur,'cur');
end

endfunction
