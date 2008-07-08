#!/usr/bin/python
import numpy
import scipy.io
import scipy.optimize

alpha=.1   
#For partial differences. Let x1 and x2 be consecutive x values of grid points.
#Then x2-x1=dx, then we use a partial difference of alpha*dx.  Similarly for dy,dz.
#This will be used to compute currents using finite differences.
conductivity=.0006
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
    if(xmin>xmax or nx<=0 or ymin>ymax or ny<=0 or zmin>zmax or nz<=0):
        print "Bad Arguments to MakeCubicGrid"
        return
    grid=GenerateCubicGrid(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)
    SaveCubicGrid(grid,dir+name)
    SaveCubicGrid(CreateDGrid(grid,0),dir+name+"dx")
    SaveCubicGrid(CreateDGrid(grid,1),dir+name+"dy")
    SaveCubicGrid(CreateDGrid(grid,2),dir+name+"dz")
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
def SaveInjVTK(inj,filename):
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
    #Careful.  This function requires the field as a rank 2 array (that is, a Matrix)
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
    pot,cur=GetPotentialAndCurrent(inj,geom,grid,gains)
    curmagn=GetCurrentMagnitude(cur)
    epsilon=1e-7
#    cursq=GetCurSq(cur)
    SaveInjVTK(inj,fileprefix+"_inj.vtk")
#   dotsqwithJ0=GetDotSQWithJ0(inj,geom,grid,gains)
#    signdotJ0=GetSignDotWithJ0(inj,geom,grid,gains)
    normalcurrent=NormalCurrent(cur,geom,grid)
    
    SaveTrimmedFieldVTK(grid,TrimFieldNerve(geom,grid,curmagn),fileprefix+"_cmag_nerve.vtk","Current_Magnitude",epsilon)
    SaveTrimmedFieldVTK(grid, TrimFieldFocus(geom,grid,curmagn),fileprefix+"_cmag_focus.vtk","Current_Magnitude",epsilon)
    SaveTrimmedFieldVTK(grid, TrimFieldNerve(geom,grid,cur),fileprefix+"_cur_nerve.vtk","Current",epsilon)
    SaveTrimmedFieldVTK(grid, TrimFieldNerve(geom,grid,pot),fileprefix+"_pot_nerve.vtk","Potential",epsilon)
    SaveTrimmedFieldVTK(grid, TrimFieldNerve(geom,grid,normalcurrent),fileprefix+"_cur_normal.vtk","NormalCur",epsilon)
    #SaveTrimmedFieldVTK(grid, TrimFieldCore(geom,grid,curmagn),fileprefix+"_cmag_core.vtk")
    #SaveTrimmedFieldVTK(grid, TrimFieldCore(geom,grid,cur),fileprefix+"_cur_core.vtk")    
    #SaveTrimmedFieldVTK(grid, TrimFieldNerve(geom,grid,dotsqwithJ0),fileprefix+"_dotsq_nerve.vtk")
    #SaveTrimmedFieldVTK(grid, TrimFieldCore(geom,grid,dotsqwithJ0),fileprefix+"_dotsq_core.vtk")    
    #SaveTrimmedFieldVTK(grid, TrimFieldNerve(geom,grid,signdotJ0),fileprefix+"_dotsign_nerve.vtk")
    #SaveTrimmedFieldVTK(grid, TrimFieldCore(geom,grid,signdotJ0),fileprefix+"_dotsign_core.vtk")

def GetDimensions(grid):
    nx=Get1Dimension(grid[:,0])
    ny=Get1Dimension(grid[:,1])
    nz=Get1Dimension(grid[:,2])
    return nx,ny,nz
def Get1Dimension(x):
    M=max(x)
    m=min(x)
    return int(round(1+(M-m)/GetGridSpacing(x)))
def LoadGain(geom,grid,filename):
    #Loads Gain Matrix, then uses grid information to zero out all gain elements corresponding to 
    #grid locations outside the nerve.  We do this because solver gives undefined results outside nerve.
    return TrimFieldNerve(geom,grid,scipy.io.read_array(filename))
    #return scipy.io.read_array(filename)
def LoadGains(geom,grid,fileprefix):
    #fileprefix has form like "/somewhere/nerve1.mycut"
    gain=LoadGain(geom, grid,fileprefix+".gain")
    gaindx=LoadGain(geom, grid,fileprefix+"dx.gain")
    gaindy=LoadGain(geom, grid,fileprefix+"dy.gain")
    gaindz=LoadGain(geom, grid,fileprefix+"dz.gain")
    return gain,gaindx,gaindy,gaindz
def GetPotentialAndCurrent(inj,geom,grid,gains):
    [dx,dy,dz]=GetFiniteDifferenceDxDyDz(grid)
    pot=numpy.dot(gains[0],inj)
    curx=-conductivity*(numpy.dot(gains[1],inj)-pot)/dx
    cury=-conductivity*(numpy.dot(gains[2],inj)-pot)/dy
    curz=-conductivity*(numpy.dot(gains[3],inj)-pot)/dz
    return numpy.transpose(numpy.array([pot])),numpy.transpose(numpy.array([curx,cury,curz]))

def GetGridDxDyDz(grid):
    return [GetGridSpacing(grid[:,i]) for i in range (3)]
def GetFiniteDifferenceDxDyDz(grid):
    return [alpha*GetGridSpacing(grid[:,i]) for i in range (3)]
def GetGridSpacing(x):
    v=x-x[0]
    v=numpy.unique(abs(v))
    if(v[0]==0):
        return v[1]
    else:
        return v[0]

def TrimFieldNerve(geom,grid,field):
    #If grid[i] is outside of the nerve, we set field[i]=0 (or 0,0,0 for current)
    newfield=numpy.array([field[i]*IsInsideNerve(grid[i],geom) for i in range(len(field))])
    return newfield
def TrimFieldCore(geom,grid,field):
    #If grid[i] is outside of the core, we set field[i]=0 (or 0,0,0 for current)
    newfield=numpy.array([field[i]*IsInsideCore(grid[i],geom) for i in range(len(field))])
    return newfield
def TrimFieldFocus(geom,grid,field):
    #If grid[i] is outside of the core, we set field[i]=0 (or 0,0,0 for current)
    newfield=numpy.array([field[i]*IsNearFocus(grid[i],geom) for i in range(len(field))])
    return newfield
def IsNearFocus(x,geom):
    x0=geom[2,0:3]
    r=.0*geom[2,6]
    if numpy.linalg.norm(x-x0)<r:
        return float(True)
    else:
        return float(False)
def IsInsideNerve(x,geom):
    rnerve=geom[1,4]
    lnerve=geom[1,3]
    if(numpy.linalg.norm([x[0],x[1]])<rnerve and (abs(x[2])<lnerve)):
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
        #self.geom[2,3:6] = J_0.  This MUST HAVE LENGTH 1 !!!!!!
        self.grid=LoadCubicGrid(gridfilename)
        self.RandInj()
        g0,g1,g2,g3=LoadGains(self.geom,self.grid,gainfilename)
        self.gains=numpy.array([g0,g1,g2,g3])
    def RandInj(self):
        self.cinj=numpy.random.sample(11)-.5
        self.inj=numpy.concatenate((self.cinj,[-sum(self.cinj)]))
    def ObjPhi(self,inj):
        return ObjPhi(inj,self.geom,self.grid,self.gains)
    def ConstrainedObjPhi(self,cinj):
        y=numpy.concatenate((cinj,[-sum(cinj)]))
        return self.ObjPhi(y)
    def ObjOmega(self,inj):
        return ObjOmega(inj,self.geom,self.grid,self.gains)
    def ConstrainedObjOmega(self,cinj):
        y=numpy.concatenate((cinj,[-sum(cinj)]))
        return self.ObjOmega(y)
    def ObjChi(self,inj):
        return ObjChi(inj,self.geom,self.grid,self.gains)
    def ConstrainedObjChi(self,cinj):
        y=numpy.concatenate((cinj,[-sum(cinj)]))
        return self.ObjChi(y)
    def OptimizePhi(self):
        self.RandInj()
        self.CurrentFunc=self.ObjPhi
        temp=scipy.optimize.fmin_bfgs(self.ConstrainedObjPhi,self.cinj, callback=self.MyCallback)
        self.SetInj(temp)
        return temp
    def OptimizeOmega(self):
        self.geom[2,3:6]=(1/numpy.linalg.norm(self.geom[2,3:6]))*self.geom[2,3:6]
        self.RandInj()
        self.CurrentFunc=self.ConstrainedObjOmega
        temp=scipy.optimize.fmin_bfgs(self.ConstrainedObjOmega,self.cinj,callback=self.MyCallback)
        self.SetInj(temp)
        return temp
    def OptimizeChi(self):
        self.RandInj()
        self.CurrentFunc=self.ConstrainedObjChi
        temp=scipy.optimize.fmin_bfgs(self.ConstrainedObjChi,self.cinj,callback=self.MyCallback)
        self.SetInj(temp)
        return temp
    def OptimizeOmegaGeom(self):
        self.RandInj()
        self.CurrentFunc=self.ObjOmegaGeom
        x=numpy.concatenate((self.cinj,self.geom[2,0:3]))
        temp=scipy.optimize.fmin_bfgs(self.ObjOmegaGeom,x,callback=self.MyCallback)
        self.SetInjGeom(temp)
        return temp
    def SetInj(self,cinj):
        self.inj[0:11]=cinj
        self.inj[11]=-sum(cinj)
    def SetInjGeom(self,x):
        self.cinj=x[0:11]
        self.inj[0:11]=self.cinj
        self.inj[11]=-sum(self.cinj)
        self.geom[2,0:3]=x[11:14]
    def OptimizeOmegaGeomOnly(self):
        self.RandInj()
        self.CurrentFunc=self.ObjOmegaGeomOnly
        x=self.geom[2,0:3]
        print x
        temp=scipy.optimize.fmin_bfgs(self.ObjOmegaGeomOnly,x,callback=self.MyCallback)
        cinj=temp[0:11]
        inj=numpy.concatenate((cinj,[-sum(cinj)]))
        return temp
    
    def ObjOmegaGeom(self,x):
        cinj=x[0:11]
        inj=numpy.concatenate((cinj,[-sum(cinj)]))
        g2=self.geom[2]
        g2[0:3]=x[11:14]
        geom=numpy.array([self.geom[0],self.geom[1],g2])
        return ObjOmega(inj,geom,self.grid,self.gains)
    def ObjOmegaGeomOnly(self,x):
        g2=self.geom[2]
        g2[0:3]=x[0:3]
        geom=numpy.array([self.geom[0],self.geom[1],g2])
        return ObjOmega(self.inj,geom,self.grid,self.gains)
    def RandOmegaGeom(self):
        self.geom[2,0:3]=(numpy.random.sample(3)-.5)*.5
        self.geom[2,2]=(numpy.random.sample(1)[0]-.5)*12.0
        self.geom[2,6]=.3
    def OptimizeChiGeom(self):
        self.RandInj()
        self.RandOmegaGeom()
        self.CurrentFunc=self.ObjChiGeom
        x=numpy.concatenate((self.cinj,self.geom[2,0:3]))
        temp=scipy.optimize.fmin_bfgs(self.ObjChiGeom,x,callback=self.MyCallback)
        self.SetInjGeom(temp)
        return temp
    def OptimizeChiGeomZOnly(self):
        self.RandInj()
        self.RandOmegaGeom()
        self.geom[2,0]=0.
        self.geom[2,1]=0.
        self.CurrentFunc=self.ObjChiGeomZOnly
        x=numpy.concatenate((self.cinj,self.geom[2,2]))
        temp=scipy.optimize.fmin_bfgs(self.ObjChiGeomZOnly,x,callback=self.MyCallback)
        return temp
    def OptimizeChiGeomOnly(self):
        self.RandInj()
        self.CurrentFunc=self.ObjChiGeomOnly
        x=self.geom[2,0:3]
        temp=scipy.optimize.fmin_bfgs(self.ObjChiGeomOnly,x,callback=self.MyCallback)
        return temp    
    def ObjChiGeom(self,x):
        cinj=x[0:11]
        inj=numpy.concatenate((cinj,[-sum(cinj)]))
        g2=self.geom[2]
        g2[0:3]=x[11:14]
        geom=numpy.array([self.geom[0],self.geom[1],g2])
        return ObjChi(inj,geom,self.grid,self.gains)
    def ObjChiGeomZOnly(self,x):
        cinj=x[0:11]
        inj=numpy.concatenate((cinj,[-sum(cinj)]))
        g2=self.geom[2]
        g2[2]=x[12]
        geom=numpy.array([self.geom[0],self.geom[1],g2])
        return ObjChi(inj,geom,self.grid,self.gains)
    def ObjChiGeomOnly(self,x):
        g2=self.geom[2]
        g2[0:3]=x[0:3]
        geom=numpy.array([self.geom[0],self.geom[1],g2])
        return ObjChi(self.inj,geom,self.grid,self.gains)
    
    def MyCallback(self,x):
        print "Callback."
        print "params = ", x
        print self.CurrentFunc(x)

def ObjPhi(inj,geom,grid,gains):
    a=PhiN(inj,geom,grid,gains)
    b=PhiC(inj,geom,grid,gains)
    c=VolN(geom,grid)
    d=VolC(geom,grid)
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
    return numpy.reshape(numpy.apply_along_axis(func,1,cur),(-1,1))
def GetCurSq(cur):
    func = lambda x: numpy.linalg.norm(x)**2.
    return numpy.reshape(numpy.apply_along_axis(func,1,cur),(-1,1))
def GetCoreVar(inj,geom,grid,gains):
    cur=GetPotentialAndCurrent(inj,geom,grid,gains)[1]
    temp=GetCurSq(TrimFieldCore(geom,grid,cur)).flatten()
    temp2=[x for x in temp if x > 10e-7]
    return numpy.var(temp2)    

def W(x,x0,sigma):
    return sigma**(-1)*(2*scipy.pi)**(-.5)*scipy.exp(-.5*((numpy.linalg.norm(x-x0)/sigma)**2))

def Chi(inj,geom,grid,gains):
    x0=geom[2,0:3]
    sigma=geom[2,6]
    dv=numpy.product(GetGridDxDyDz(grid))

    cur=GetPotentialAndCurrent(inj,geom,grid,gains)[1]
    Omega_i=lambda i: W(grid[i],x0,sigma)*numpy.dot(cur[i],cur[i])
    return dv*numpy.sum([Omega_i(i) for i in range(len(grid))])
def ObjChi(inj,geom,grid,gains):
    a=PhiN(inj,geom,grid,gains)
    b=Chi(inj,geom,grid,gains)
    return a/float(b)
def ObjOmega(inj,geom,grid,gains):
    a=PhiN(inj,geom,grid,gains)
    b=Omega(inj,geom,grid,gains)
    return a/float(b)
def GetDotWithJ0(inj,geom,grid,gains):
    cur=GetPotentialAndCurrent(inj,geom,grid,gains)[1]
    Dot_i=lambda i: numpy.dot(cur[i],J0)
    return numpy.array([Dot_i(i) for i in range(len(grid))]).reshape((-1,1))
def GetDotSQWithJ0(inj,geom,grid,gains):
    cur=GetPotentialAndCurrent(inj,geom,grid,gains)[1]
    Dot_i=lambda i: numpy.dot(cur[i],J0)
    return numpy.array([Dot_i(i)**(2.0) for i in range(len(grid))]).reshape((-1,1))
def GetSignDotWithJ0(inj,geom,grid,gains):
    cur=GetPotentialAndCurrent(inj,geom,grid,gains)[1]
    Dot_i=lambda i: numpy.dot(cur[i],J0)
    return numpy.array([numpy.sign(Dot_i(i)) for i in range(len(grid))]).reshape((-1,1))
def VolN(geom,grid):
    return scipy.pi*geom[1,3]*geom[1,4]**2
def VolC(geom,grid):
    return scipy.pi*geom[0,3]*geom[0,4]**2

def Omega(inj,geom,grid,gains):
    x0=geom[2,0:3]
    sigma=geom[2,6]
    J0=geom[2,3:6]
    dv=numpy.product(GetGridDxDyDz(grid))

    cur=GetPotentialAndCurrent(inj,geom,grid,gains)[1]
    omega_i=lambda i: W(grid[i],x0,sigma)*(numpy.linalg.norm(cur[i])**2.0-numpy.linalg.norm(cur[i])*numpy.dot(cur[i],J0))
    return dv*numpy.sum([omega_i(i) for i in range(len(grid))])
def Normalize(inj,geom,grid,gains):
    return numpy.array((1e10/PhiN(inj,geom,grid,gains))**.5 *inj,float)
def InjFromCinj(cinj):
    temp=numpy.zeros(12)
    temp[0:11]=cinj
    temp[11]=-sum(cinj)
    return numpy.array(temp,float)
def NormalCurrent(cur,geom,grid):
    Proj_i=lambda i: numpy.dot(cur[i,0:2],grid[i,0:2])/numpy.linalg.norm(grid[i,0:2])**2.0
    return numpy.array([Proj_i(i) for i in range(len(grid))]).reshape((-1,1))
    