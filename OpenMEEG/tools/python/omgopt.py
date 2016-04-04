#!/usr/bin/python
import numpy
import re
import scipy.io
import scipy.optimize

alpha=.1   
#For partial differences. Let x1 and x2 be consecutive x values of grid points.
#Then x2-x1=dx, then we use a partial difference of alpha*dx.  Similarly for dy,dz.
#Thus, alpha is used to compute currents using finite differences on potentials.

conductivity=.0006
#This constant should be changed to match your model!

def GenerateCubicGrid(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz):
    if(xmin>xmax or nx<=0 or ymin>ymax or ny<=0 or zmin>zmax or nz<=0):
        print "Bad Arguments to MakeCubicGrid"
        return
    return numpy.resize(GenerateMGrid(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz),(nx*ny*nz,3))
def GenerateMGrid(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz):
    if(xmin>xmax or nx<=0 or ymin>ymax or ny<=0 or zmin>zmax or nz<=0):
        print "Bad Arguments to MakeCubicGrid"
        return
    return numpy.mgrid[xmin:xmax:1j*nx,ymin:ymax:1j*ny,zmin:zmax:1j*nz].swapaxes(0,1).swapaxes(1,2).swapaxes(2,3)
def CubicGridToMGrid(grid,nx,ny,nz):
    return numpy.resize(grid,(nx,ny,nz,3))
def SaveCubicGrid(grid,filename):
    file=open(filename,'w')
    scipy.io.write_array(file,grid)
    file.close()
def LoadCubicGrid(filename):
    file=open(filename,'r')
    return scipy.io.read_array(file)
def SaveGridFiles(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz,dir,name):
    y=re.compile('/\Z')
    if(y.match(dir) is None):
        dir=dir+"/"
    if(xmin>xmax or nx<=0 or ymin>ymax or ny<=0 or zmin>zmax or nz<=0):
        print "Bad Arguments to MakeCubicGrid"
        return
    grid=GenerateCubicGrid(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)
    SaveCubicGrid(grid,dir+name)
    SaveCubicGrid(CreateDGrid(grid,0),dir+name+"dx")
    SaveCubicGrid(CreateDGrid(grid,1),dir+name+"dy")
    SaveCubicGrid(CreateDGrid(grid,2),dir+name+"dz")
    SaveCubicGrid(CreateNegDGrid(grid,2),dir+name+"-dz")
    SaveGridVTK(grid,nx,ny,nz,dir+name+".vtk")
    
def CreateDGrid(grid,index):
    if(index not in set([0,1,2])):
        print "Index must be 0,1, or 2\n"
        sys.exit()
    delta=GetGridSpacing(grid[:,index])
    tileunit=[0,0,0]
    tileunit[index]=1
    vec=numpy.tile(numpy.array(tileunit),(len(grid),1))
    return grid+alpha*delta*vec
def CreateNegDGrid(grid,index):
    if(index not in set([0,1,2])):
        print "Index must be 0,1, or 2\n"
        sys.exit()
    delta=GetGridSpacing(grid[:,index])
    tileunit=[0,0,0]
    tileunit[index]=1
    vec=numpy.tile(numpy.array(tileunit),(len(grid),1))
    return grid-alpha*delta*vec
def SaveInjVTK(inj,filename):
    #Saves a VTK file with 12 glyphs representing injected current.  Locations are hardwired.
    file=open(filename,'w')
    N=12
    locations=numpy.array([
       [  1.11    ,   0.00    ,  -6.00    ],
       [  0.00    ,   1.11    ,  -6.00    ],
       [ -1.11    ,   0.00    ,  -6.00    ],
       [  0.00    ,  -1.11    ,  -6.00    ],
       [  1.11    ,   0.00    ,   0.00    ],
       [  0.00    ,   1.11    ,   0.00    ],
       [ -1.11    ,   0.00    ,   0.00    ],
       [  0.00    ,  -1.11    ,   0.00    ],
       [  1.11    ,   0.00    ,   6.00    ],
       [  0.00    ,   1.11    ,   6.00    ],
       [ -1.11    ,   0.00    ,   6.00    ],
       [  0.00    ,  -1.11    ,   6.00    ]])
    file.write("# vtk DataFile Version 2.0\n")
    file.write(filename+"\n")
    file.write("ASCII\n")
    file.write("DATASET POLYDATA\n")
    file.write("POINTS "+str(N)+" float\n")    
    scipy.io.write_array(file,locations,keep_open=True)
    file.write("POINT_DATA "+str(N)+"\n")
    file.write("SCALARS Injected_Current float 1\n")
    file.write("LOOKUP_TABLE default\n")
    scipy.io.write_array(file,inj)
    file.close()
def SaveGridVTK(grid,nx,ny,nz,filename):
    file=open(filename,'w')
    file.write("# vtk DataFile Version 2.0\n")
    file.write(filename+"\n")
    file.write("ASCII\n")
    file.write("DATASET STRUCTURED_GRID\n")
    file.write("DIMENSIONS "+str(nx)+" "+str(ny)+" "+str(nz)+"\n")
    file.write("POINTS "+str(nx*ny*nz)+" float\n")    
    scipy.io.write_array(file,grid.reshape((nx,ny,nz,3)).reshape((nx*ny*nz,3),order="F"))
    file.close()
def SavePolyVTK(grid,N,filename):
    file=open(filename,'w')
    file.write("# vtk DataFile Version 2.0\n")
    file.write(filename+"\n")
    file.write("ASCII\n")
    file.write("DATASET POLYDATA\n")
    file.write("POINTS "+str(N)+" float\n")    
    scipy.io.write_array(file,grid.reshape((N,3),order="F"))
    file.close()
def SaveTrimmedFieldVTK(gridxyz,field,filename,FieldName,epsilon):
    #Careful.  This function requires the field as a rank 2 array (that is, a matrix)
    #This way a scalar field is written as a vector field of 1 component vectors
    #In particular, it allows the same framework to apply to both potential and current
    m=len(field.transpose())
    indices=numpy.where(numpy.apply_along_axis(numpy.linalg.norm,1,field)>epsilon)
    N=numpy.size(indices)
    tempgrid=numpy.array(gridxyz[indices])
    tempfield=numpy.array(field[indices])
    
    SavePolyVTK(tempgrid,N,filename)
    file=open(filename,'a')
    file.write("POINT_DATA "+str(N)+"\n")
    if(m==1):
        file.write("SCALARS "+FieldName+" float 1\n")
        file.write("LOOKUP_TABLE default\n")
    else:
        file.write("VECTORS "+FieldName+" float\n")
    
    scipy.io.write_array(file,tempfield.reshape((N,m),order="F"))
def SaveFieldsVTK(inj,geom,grid,gains,fileprefix):
    #Saves a bunch of VTK files for visualization.
    #Ie current, current magnitude, potential
    pot,cur=GetPotentialAndCurrent(inj,geom,grid,gains)
    curmagn=GetCurrentMagnitude(cur)
    epsilon=1e-7
    SaveInjVTK(inj,fileprefix+"_inj.vtk")
    
    SaveTrimmedFieldVTK(grid,TrimFieldNerve(geom,grid,curmagn),fileprefix+"_cmag_nerve.vtk","Current_Magnitude",epsilon)
    SaveTrimmedFieldVTK(grid, TrimFieldFocus(geom,grid,curmagn),fileprefix+"_cmag_focus.vtk","Current_Magnitude",epsilon)
    SaveTrimmedFieldVTK(grid, TrimFieldNerve(geom,grid,cur),fileprefix+"_cur_nerve.vtk","Current",epsilon)
    SaveTrimmedFieldVTK(grid, TrimFieldNerve(geom,grid,pot),fileprefix+"_pot_nerve.vtk","Potential",epsilon)
    SaveTrimmedFieldVTK(grid, TrimFieldFocus(geom,grid,cur),fileprefix+"_cur_focus.vtk","Current_Focus",epsilon)
    
    #SaveTrimmedFieldVTK(grid, TrimFieldCore(geom,grid,curmagn),fileprefix+"_cmag_core.vtk")
    #SaveTrimmedFieldVTK(grid, TrimFieldCore(geom,grid,cur),fileprefix+"_cur_core.vtk")    
def GetDimensions(grid):
    #Get the number of grid points in each dimension of the grid
    nx=Get1Dimension(grid[:,0])
    ny=Get1Dimension(grid[:,1])
    nz=Get1Dimension(grid[:,2])
    return nx,ny,nz
def Get1Dimension(x):
    M=max(x)
    m=min(x)
    return int(round(1+(M-m)/GetGridSpacing(x)))
def LoadGain(geom,grid,filename):
    #Loads Gain matrix, then uses grid information to zero out all gain elements corresponding to 
    #grid locations outside the nerve.  We do this because solver gives undefined results outside nerve.
    return TrimFieldNerve(geom,grid,scipy.io.read_array(filename))
    #return scipy.io.read_array(filename)
def LoadGains(geom,grid,fileprefix):
    #fileprefix has form like "/somewhere/nerve1.mycut"
    gain=LoadGain(geom, grid,fileprefix+".gain")
    gaindx=LoadGain(geom, grid,fileprefix+"dx.gain")
    gaindy=LoadGain(geom, grid,fileprefix+"dy.gain")
    gaindz=LoadGain(geom, grid,fileprefix+"dz.gain")
    gainminusdz=LoadGain(geom, grid,fileprefix+"-dz.gain")
    return gain,gaindx,gaindy,gaindz,gainminusdz
def GetPotentialAndCurrent(inj,geom,grid,gains):
    [dx,dy,dz]=GetFiniteDifferenceDxDyDz(grid)
    pot=numpy.dot(gains[0],inj)
    curx=-conductivity*(numpy.dot(gains[1],inj)-pot)/dx
    cury=-conductivity*(numpy.dot(gains[2],inj)-pot)/dy
    curz=-conductivity*(numpy.dot(gains[3],inj)-pot)/dz
    return numpy.transpose(numpy.array([pot])),numpy.transpose(numpy.array([curx,cury,curz]))
def GetActivationFunction(inj,geom,grid,gains):
    [dx,dy,dz]=GetFiniteDifferenceDxDyDz(grid)
    activation=numpy.dot(gains[3]+gains[4]-2*gains[0],inj)/(dz*dz)
    return numpy.transpose(numpy.array([activation]))
def GetGridDxDyDz(grid):
    return [GetGridSpacing(grid[:,i]) for i in range (3)]
def GetGridSpacing(x):
    v=x-x[0]
    v=numpy.unique(abs(v))
    if(v[0]==0):
        return v[1]
    else:
        return v[0]
def GetFiniteDifferenceDxDyDz(grid):
    #These are the deltas used to calculate currents from potentials via finite differences and ohms law.
    return [alpha*GetGridSpacing(grid[:,i]) for i in range (3)]
def TrimFieldNerve(geom,grid,field):
    #If grid[i] is outside of the nerve region, we set field[i]=0 (or 0,0,0 for current)
    newfield=numpy.array([field[i]*IsInsideNerve(grid[i],geom) for i in range(len(field))])
    return newfield
def TrimFieldCore(geom,grid,field):
    #If grid[i] is outside of the core region, we set field[i]=0 (or 0,0,0 for current)
    newfield=numpy.array([field[i]*IsInsideCore(grid[i],geom) for i in range(len(field))])
    return newfield
def TrimFieldFocus(geom,grid,field):
    #If grid[i] is outside of the focus region, we set field[i]=0 (or 0,0,0 for current)
    newfield=numpy.array([field[i]*IsNearFocus(grid[i],geom) for i in range(len(field))])
    return newfield
def IsNearFocus(x,geom):
    x0=geom[2,0:3]
    r=geom[2,6]
    if numpy.linalg.norm(x-x0)<r:
        return float(True)
    else:
        return float(False)
def IsInsideNerve(x,geom):
    #rnerve=geom[1,4]
    #lnerve=geom[1,3]
    #if(numpy.linalg.norm([x[0],x[1]])<rnerve and (abs(x[2])<lnerve)):
    x0=geom[1,0]
    y0=geom[1,1]
    z0=geom[1,2]
    l=geom[1,3]
    r=geom[1,4]
    if(numpy.linalg.norm([x[0]-x0,x[1]-y0])<=r and (abs(x[2]-z0) <=l )):
        return float(True)
    else:
        return float(False)
def IsInsideCore(x,geom):
    x0=geom[0,0]
    y0=geom[0,1]
    z0=geom[0,2]
    l=geom[0,3]
    r=geom[0,4]
    if(numpy.linalg.norm([x[0]-x0,x[1]-y0])<=r and (abs(x[2]-z0) <=l )):
        return float(True)
    else:
        return float(False)
    
class workspace:
    def __init__(self,gridfilename,gainfilename):
        self.geom=numpy.array([[0.,0.,0.,5.,.3,0.0,0.0],[0.,0.,0.,12.,.95,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,.3]])
        
        #self.geom[2,3:6] = J0.  This MUST HAVE LENGTH 1 !
        #The geometry matrix has form (* denotes unused)
        #geom[0]=CORE= [x,y,z,l,r,*,*]
        #geom[1]=NERVE = [x,y,z,l,r,*,*]
        #geom[2]=Focus/Chi/Omega = [x,y,z,J0_x,J0_y,J0_z,sigma]
        self.grid=LoadCubicGrid(gridfilename)

        g0,g1,g2,g3,g4=LoadGains(self.geom,self.grid,gainfilename)   
        self.gains=numpy.array([g0,g1,g2,g3,g4])
        self.NumberOfElectrodes=len(self.gains[0,0])
        self.ConstrainedNumberOfElectrodes=self.NumberOfElectrodes-1
        self.SetRandomInj()
        self.GTol=.0005 #Tolerance (for norm of gradient) for when to stop optimization iterations.  
    def SetRandomInj(self): #randomize the injection current
        self.cinj=numpy.random.sample(self.ConstrainedNumberOfElectrodes)-.5 #Constrained injected current: only the first N-1 positions.
        self.inj=numpy.concatenate((self.cinj,[-sum(self.cinj)]))
        alpha=(1/numpy.linalg.norm(self.inj))
        self.cinj=self.cinj*alpha
        self.inj=self.inj*alpha
    def f_Phi(self,inj):
        return f_Phi(inj,self.geom,self.grid,self.gains)
    def Constrained_f_Phi(self,cinj):
        y=numpy.concatenate((cinj,[-sum(cinj)]))
        return self.f_Phi(y)
    def f_Omega(self,inj):
        return f_Omega(inj,self.geom,self.grid,self.gains)
    def Constrained_f_Omega(self,cinj):
        y=numpy.concatenate((cinj,[-sum(cinj)]))
        return self.f_Omega(y)
    def f_Chi(self,inj):
        return f_Chi(inj,self.geom,self.grid,self.gains)
    def Constrained_f_Chi(self,cinj):
        y=numpy.concatenate((cinj,[-sum(cinj)]))
        return self.f_Chi(y)
    def f_Ksi(self,inj):
        return f_Ksi(inj,self.geom,self.grid,self.gains)
    def Constrained_f_Ksi(self,cinj):
        y=numpy.concatenate((cinj,[-sum(cinj)]))
        return self.f_Ksi(y)
    def OptimizePhi(self):
        self.SetRandomInj()
        self.CurrentFunc=self.Constrained_f_Phi
        temp=scipy.optimize.fmin_bfgs(self.Constrained_f_Phi,self.cinj, callback=self.MyCallback,gtol=self.GTol)
        self.SetInj(temp)
        return temp
    def OptimizeOmega(self):
        self.geom[2,3:6]=(1/numpy.linalg.norm(self.geom[2,3:6]))*self.geom[2,3:6]
        self.SetRandomInj()
        self.CurrentFunc=self.Constrained_f_Omega
        temp=scipy.optimize.fmin_bfgs(self.Constrained_f_Omega,self.cinj,callback=self.MyCallback,gtol=self.GTol)
        self.SetInj(temp)
        return temp
    def OptimizeChi(self):
        self.SetRandomInj()
        self.CurrentFunc=self.Constrained_f_Chi
        temp=scipy.optimize.fmin_bfgs(self.Constrained_f_Chi,self.cinj,callback=self.MyCallback,gtol=self.GTol)
        self.SetInj(temp)
        return temp
    def OptimizeKsi(self):
        self.SetRandomInj()
        self.CurrentFunc=self.Constrained_f_Ksi
        temp=scipy.optimize.fmin_bfgs(self.Constrained_f_Ksi,self.cinj,retall=1,callback=self.MyCallback,gtol=self.GTol)
        return temp
    def OptimizeOmegaGeom(self):
        self.SetRandomInj()
        self.CurrentFunc=self.f_OmegaGeom
        x=numpy.concatenate((self.cinj,self.geom[2,0:3]))
        temp=scipy.optimize.fmin_bfgs(self.f_OmegaGeom,x,callback=self.MyCallback,gtol=self.GTol)
        self.SetInjGeom(temp)
        return temp
    def SetInj(self,cinj):
        self.inj[0:self.ConstrainedNumberOfElectrodes]=cinj
        self.inj[self.ConstrainedNumberOfElectrodes]=-sum(cinj)
    def SetInjGeom(self,x):
        self.cinj=x[0:self.ConstrainedNumberOfElectrodes]
        self.inj[0:self.ConstrainedNumberOfElectrodes]=self.cinj
        self.inj[self.ConstrainedNumberOfElectrodes]=-sum(self.cinj)
        self.geom[2,0:3]=x[self.ConstrainedNumberOfElectrodes:self.ConstrainedNumberOfElectrodes+3]
    def f_OmegaGeom(self,x):
        
        cinj=x[0:self.ConstrainedNumberOfElectrodes]
        inj=numpy.concatenate((cinj,[-sum(cinj)]))
        g2=self.geom[2]
        g2[0:3]=x[self.ConstrainedNumberOfElectrodes:self.ConstrainedNumberOfElectrodes+3]
        geom=numpy.array([self.geom[0],self.geom[1],g2])
        return f_Omega(inj,geom,self.grid,self.gains)
    def SetRandomOmegaGeom(self):
        self.geom[2,0:3]=(numpy.random.sample(3)-.5)*.5
        self.geom[2,2]=(numpy.random.sample(1)[0]-.5)*12.0
        self.geom[2,6]=.3
    def OptimizeChiGeom(self):
        self.SetRandomInj()
        self.SetRandomOmegaGeom()
        self.CurrentFunc=self.f_ChiGeom
        x=numpy.concatenate((self.cinj,self.geom[2,0:3]))
        temp=scipy.optimize.fmin_bfgs(self.f_ChiGeom,x,callback=self.MyCallback,gtol=self.GTol)
        self.SetInjGeom(temp)
        return temp    
    def f_ChiGeom(self,x):
        cinj=x[0:self.ConstrainedNumberOfElectrodes]
        inj=numpy.concatenate((cinj,[-sum(cinj)]))
        g2=self.geom[2]
        g2[0:3]=x[self.ConstrainedNumberOfElectrodes:self.ConstrainedNumberOfElectrodes+3]
        geom=numpy.array([self.geom[0],self.geom[1],g2])
        return f_Chi(inj,geom,self.grid,self.gains)
    def OptimizeKsiGeom(self):
        self.SetRandomInj()
        self.SetRandomOmegaGeom()
        self.CurrentFunc=self.f_KsiGeom
        x=numpy.concatenate((self.cinj,self.geom[2,0:3]))
        temp=scipy.optimize.fmin_bfgs(self.f_KsiGeom,x,callback=self.MyCallback,gtol=self.GTol)
        self.SetInjGeom(temp)
        return temp
    def f_KsiGeom(self,x):
        cinj=x[0:self.ConstrainedNumberOfElectrodes]
        inj=numpy.concatenate((cinj,[-sum(cinj)]))
        g2=self.geom[2]
        g2[0:3]=x[self.ConstrainedNumberOfElectrodes:self.ConstrainedNumberOfElectrodes+3]
        geom=numpy.array([self.geom[0],self.geom[1],g2])
        return f_Ksi(inj,geom,self.grid,self.gains)
    def MyCallback(self,x):
        print "Callback."
        print "params = ", x
        print self.CurrentFunc(x)
        
def f_Phi(inj,geom,grid,gains):
    a=PhiN(inj,geom,grid,gains)
    b=PhiC(inj,geom,grid,gains)
    c=VolumeNerve(geom,grid)
    d=VolumeCore(geom,grid)
    return a*d/(float(b)*float(c))

def PhiN(inj,geom,grid,gains):
    dv=numpy.product(GetGridDxDyDz(grid))
    cur=GetPotentialAndCurrent(inj,geom,grid,gains)[1]
    return dv*numpy.sum([numpy.dot(cur[x],cur[x]) for x in range(len(grid))])
def PhiC(inj,geom,grid,gains):
    dv=numpy.product(GetGridDxDyDz(grid))
    cur=GetPotentialAndCurrent(inj,geom,grid,gains)[1]
    return dv*numpy.sum([IsInsideCore(grid[x],geom)*numpy.dot(cur[x],cur[x]) for x in range(len(grid))])

def GetCurrentMagnitude(cur):
    func=lambda x: numpy.linalg.norm(x)
    return numpy.reshape(numpy.apply_along_axis(func,1,cur),(-1,1)) # -1 means unspecified value - inferred from the data
def GetCurSq(cur):
    func = lambda x: numpy.linalg.norm(x)**2.
    return numpy.reshape(numpy.apply_along_axis(func,1,cur),(-1,1))

def W(x,x0,sigma):
    return sigma**(-1)*(2*scipy.pi)**(-.5)*scipy.exp(-.5*((numpy.linalg.norm(x-x0)/sigma)**2))

def Chi(inj,geom,grid,gains):
    x0=geom[2,0:3]
    sigma=geom[2,6]
    dv=numpy.product(GetGridDxDyDz(grid))

    cur=GetPotentialAndCurrent(inj,geom,grid,gains)[1]
    Omega_i=lambda i: W(grid[i],x0,sigma)*numpy.dot(cur[i],cur[i])
    return dv*numpy.sum([Omega_i(i) for i in range(len(grid))])
def Ksi(inj,geom,grid,gains):
    x0=geom[2,0:3]
    sigma=geom[2,6]
    dv=numpy.product(GetGridDxDyDz(grid))
    activ=GetActivationFunction(inj,geom,grid,gains)
    Omega_i=lambda i: W(grid[i],x0,sigma)*(activ[i]**2)
    return dv*numpy.sum([Omega_i(i) for i in range(len(grid))])
def f_Chi(inj,geom,grid,gains):
    a=PhiN(inj,geom,grid,gains)
    b=Chi(inj,geom,grid,gains)
    return a/float(b)
def f_Omega(inj,geom,grid,gains):
    a=PhiN(inj,geom,grid,gains)
    b=Omega(inj,geom,grid,gains)
    return a/float(b)
def f_Ksi(inj,geom,grid,gains):
    a=PhiN(inj,geom,grid,gains)
    b=Ksi(inj,geom,grid,gains)
    return 1000000*a/float(b)
def VolumeNerve(geom,grid):
    return scipy.pi*geom[1,3]*geom[1,4]**2
def VolumeCore(geom,grid):
    return scipy.pi*geom[0,3]*geom[0,4]**2
def Omega(inj,geom,grid,gains):
    x0=geom[2,0:3]
    sigma=geom[2,6]
    J0=geom[2,3:6]
    dv=numpy.product(GetGridDxDyDz(grid))
    cur=GetPotentialAndCurrent(inj,geom,grid,gains)[1]
    omega_i=lambda i: W(grid[i],x0,sigma)*(numpy.dot(cur[i],J0)**2)
    return dv*numpy.sum([omega_i(i) for i in range(len(grid))])
def Normalize(inj,geom,grid,gains):
    #A way to scale injected currents.  This should produce comparable current densities throughout the nerve.
    return numpy.array((1/PhiN(inj,geom,grid,gains))**.5 *inj,float)
def InjFromCinj(cinj,NumberOfElectrodes):
    N=NumberOfElectrodes
    CN=N-1
    temp=numpy.zeros(N)
    temp[0:CN]=cinj
    temp[CN]=-sum(cinj)
    return numpy.array(temp,float)
def SymmetricalMatch(inj1,inj2):
    # try to find a transformation T (combining rotations and symmetries)
    # that minimizes the L2 norm between inj1 and T(inj2)
    bestinj2=inj2;
    bestnorm=numpy.norm(inj1-bestinj2,2)
    newinj2=inj2;
    for mirrortype in range(0,4):
        newinj2=Mirror(newinj2,mirrortype)
        for i in range(0,4):
            newinj2=Rotation(newinj2)
            newnorm=numpy.norm(inj1-newinj2,2)
            if (newnorm<bestnorm):
                bestnorm=newnorm
                bestinj2=newinj2
    return bestinj2

def Rotation(v):
    mat=numpy.diag(numpy.ones(11),1)
    mat[11,0]=1
    return mat*v
def Mirror(v,type):
    if type==0:
        #type 0: identity
        mat=numpy.diag(numpy.ones(12))
    elif type==1:
        # type 1: mirror through a horizontal plane containing z
        mat=numpy.zeros((12,12));
        for i in range(0,4):
            mat[i,4-i]=1
        for i in range(4,8):
            mat[i,8-i]=1
        for i in range(8,12):
            mat[i,12-i]=1
    elif type==2:
        #type 2: mirror through a vertical plane containing z
        mat=numpy.zeros((12,12));
        for i in range(0,6):
            mat[2*i,2*i]=1
        for i in range(1,6):
            mat[2*i-1,2*i+1]=1
            mat[2*i+1,2*i-1]=1
    elif type==3:
        # type 3: mirror through the central electrode (only if intercuff distances
        # are equal)
        mat=numpy.zeros((12,12))
        for i in range(0,4):
            mat[i,i+8]=1
        for i in range(4,8):
            mat[i,i]=1
        for i in range(8,12):
            mat[i,i-8]=1
    return mat*v
