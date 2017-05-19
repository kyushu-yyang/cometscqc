# filename: XPostLoading.pyx

import math
import numpy as np
cimport numpy as np
cimport cython

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

## unit:
m  = 1.
cm = 1e-2
mm = 1e-3

sec = 1.
msec = 1e-3
usec = 1e-6

## enumerate name
cpdef enum GeoType:
    kConductor = 0
    kStrip = 1
    kShell = 2

cpdef enum PosType:
    kZ
    kPhi
    kR

cpdef enum MatType:
    kRRR
    kConductivity
    kCapacity
    kField
    kTemperature
    kDose


## class to contain the geometry info
cdef class XGeoInfo:

    ## private
    cdef int _node
    cdef GeoType _geo
    cdef np.ndarray _id
    cdef np.ndarray _pos
    cdef np.ndarray _len

    ## constructor
    def __init__(self):
        self._node = 0
        self._geo  = kConductor
        self._id   = np.array([0, 0, 0], dtype=int)
        self._pos  = np.array([0., 0., 0.])
        self._len  = np.array([0., 0., 0.])

    ## setup node id
    def SetNode(self, int node):
        self._node = node

    ## return the node id
    def GetNode(self):
        return self._node

    ## setup cell's 3d id
    def SetId(self, int z, int phi, int r):
        self._id[0] = z
        self._id[1] = phi
        self._id[2] = r

    ## return cell's 3d id
    def GetId(self, PosType dim):
        if dim==kZ:
            return self._id[0]
        elif dim==kPhi:
            return self._id[1]
        elif dim==kR:
            return self._id[2]
        else:
            raise
        
    ## setup position
    def SetPosition(self, double z, double phi, double r):
        self._pos[0] = z
        self._pos[1] = phi
        self._pos[2] = r

    ## return position
    def GetPosition(self, PosType dim):
        if dim==kZ:
            return self._pos[0]
        elif dim==kPhi:
            return self._pos[1]
        elif dim==kR:
            return self._pos[2]
        else:
            raise

    ## setup cell width
    def SetWidth(self, double z, double phi, double r):
        self._len[0] = z
        self._len[1] = phi
        self._len[2] = r

    ## return cell width
    def GetWidth(self, PosType dim):
        if dim==kZ:
            return self._len[0]
        elif dim==kPhi:
            return self._len[1]
        elif dim==kR:
            return self._len[2]
        else:
            raise

    ## setup geometry enumerate
    def SetGeometry(self, int geo):
        if geo==0:
            self._geo = kConductor
        elif geo==1:
            self._geo = kStrip
        elif geo==2:
            self._geo = kShell
        else:
            raise

    ## return the geometry enumerate
    def GetGeometry(self):
        return self._geo

    ## return geometry name
    def GetGeoName(self):
        if self._geo==kConductor:
            return "Conductor"
        elif self._geo==kStrip:
            return "Strip"
        elif self._geo==kShell:
            return "Shell"
        else:
            raise


## material property collection
cdef class XMatInfo:

    ## private
    cdef int    _id
    cdef double _T
    cdef double _RRR
    cdef double _C
    cdef double _B
    cdef double _dose
    cdef np.ndarray _k

    ## constructor
    def __init__(self):
        self._id   = 0
        self._T    = 0.
        self._RRR  = 0.
        self._C    = 0.
        self._B    = 0.
        self._dose = 0.
        self._k    = np.array([0., 0., 0.])

    ## setup node id
    def SetNode(self, int node):
        self._id = node

    ## return the node id
    def GetNode(self):
        return self._id

    ## setup temperature
    def SetTemperature(self, double T):
        self._T = T

    ## return the coil temperature
    def GetTemperature(self):
        return self._T

    ## setup coil magnetic field
    def SetField(self, double B):
        self._B = B

    ## return coil magnetic field
    def GetField(self):
        return self._B

    ## setup coil RRR value
    def SetRRR(self, double RRR):
        self._RRR = RRR

    ## return coil RRR
    def GetRRR(self):
        return self._RRR

    ## setup heat capacity
    def SetCapacity(self, double C):
        self._C = C

    ## return heat capacity
    def GetCapacity(self):
        return self._C

    ## setup thermal conductivity
    def SetConductivity(self, double kz, double kp, double kr):
        self._k[0] = kz
        self._k[1] = kp
        self._k[2] = kr

    ## return thermal conductivity
    def GetConductivity(self, PosType dim):
        if dim==kZ:
            return self._k[0]
        elif dim==kPhi:
            return self._k[1]
        elif dim==kR:
            return self._k[2]
        else:
            raise

    ## setup energy deposition
    def SetDose(self, double dose):
        self._dose = dose

    ## return energy deposition
    def GetDose(self):
        return self._dose



## class to load the geometry file
cdef class XGeoFileLoad:

    ## private
    cdef np.ndarray _gc
    cdef str _name
    cdef np.ndarray _msh

    ## constructor
    def __init__(self, char* filename):
        thisfile = open(filename, "r")
        dataline = thisfile.readlines()
        thisfile.close()

        self._name = ""
        self._gc = np.array([])
        self._msh = np.array([0, 0, 0], dtype=int)

        self.fillcollect(dataline)
        self._msh[0] = self.checkmax(kZ)
        self._msh[1] = self.checkmax(kPhi)
        self._msh[2] = self.checkmax(kR)
        

    ## fill the data into collection
    def fillcollect(self, data):
        cdef XGeoInfo geo
        for eachline in data:
            eachline.strip()
            item = eachline.split()
            geo = XGeoInfo()
            geo.SetNode( int(item[0]) )
            geo.SetId( float(item[1]), float(item[2]), float(item[3]) )
            geo.SetPosition( float(item[4])/mm, float(item[5])/mm, float(item[6])/mm )
            geo.SetWidth( float(item[7])/mm, float(item[8])/mm, float(item[9])/mm )
            geo.SetGeometry( int(item[10]) )
            self._gc = np.append( self._gc, geo )


    ## check the maximum id
    cpdef int checkmax(self, PosType dim):
        cdef np.ndarray collect = np.array([], dtype=int)
        for eachgeo in self._gc:
            collect = np.append( collect, eachgeo.GetId(dim) )
        return max(collect)


    ## return the mesh size of 3 dimensions
    def GetMesh(self, PosType dim):
        if dim==kZ:
            return self._msh[0]
        elif dim==kPhi:
            return self._msh[1]
        elif dim==kR:
            return self._msh[2]
        else:
            raise


    ## setup name of magnet
    def SetName(self, str name):
        self._name = name


    ## return the magnet name
    def GetName(self):
        return self._name


    ## return the geometry collection
    def GetCollection(self):
        return self._gc


    ## check the width limits for given dimension
    def CheckLimits(self, PosType dim):
        cdef double lower
        cdef double upper
        cdef int msh = self.GetMesh(dim)
        for eachgeo in self._gc:
            if eachgeo.GetId(dim)==msh:
                upper = eachgeo.GetPosition(dim) + 1.2*eachgeo.GetWidth(dim)
                break
        for eachgeo in self._gc:
            if eachgeo.GetId(dim)==1:
                lower = eachgeo.GetPosition(dim) - 0.2*eachgeo.GetWidth(dim)
                break
        return lower, upper


    ## return the collect class for given node id
    def GetCollect(self, int node):
        for geo in self._gc:
            if geo.GetNode()==node:
                return geo



## loading the material property file
cdef class XMatFileLoad:

    ## private
    cdef np.ndarray _mc

    ## constructor
    def __init__(self, str filename):
        thisfile = open(filename, "r")
        dataline = thisfile.readlines()
        thisfile.close()
        self._mc = np.array([])

        self.fillcollect(dataline)


    ## fill the data into collection
    def fillcollect(self, data):
        cdef XMatInfo mat
        for eachline in data:
            eachline.strip()
            item = eachline.split()
            mat = XMatInfo()
            mat.SetNode( int(item[0]) )
            mat.SetTemperature( float(item[1]) )
            mat.SetRRR( float(item[2]) )
            mat.SetField( float(item[3]) )
            mat.SetCapacity( float(item[4]) )
            mat.SetConductivity( float(item[5]), float(item[6]), float(item[7]) )
            mat.SetDose( float(item[8]) )
            self._mc = np.append(self._mc, mat)


    ## return the material collection
    def GetCollection(self):
        return self._mc

    
    ## return the collect class for given node id
    def GetCollect(self, int Id):
        for mat in self._mc:
            if mat.GetNode()==Id:
                return mat



## class to plot the zr 2d plot
cdef class XPost2dPlot:

    ## private
    cdef np.ndarray _patch
    cdef np.ndarray _dset
    cdef XGeoFileLoad _geof
    cdef XMatFileLoad _matf
    cdef int _phi
    cdef PosType _direct
    cdef MatType _info
    cdef char* _name
    cdef char* _save
    cdef np.ndarray _xlim
    cdef np.ndarray _ylim
    cdef double _maxtemp
    cdef double _min
    cdef double _max

    ## constructor
    def __init__(self, char* geofile, char* datfile):
        self._patch = np.array([])
        self._dset  = np.array([])
        self._geof = XGeoFileLoad(geofile)
        self._matf = XMatFileLoad(datfile)
        self._phi = 1
        self._direct = kZ
        self._info = kConductivity
        self._name = ""
        self._save = ""
        self._xlim = np.array([0., 0.])
        self._ylim = np.array([0., 0.])
        self._maxtemp = 0.
        self._max = 0.
        self._min = 4.5


    ## set the range of x
    cpdef void setxrange(self, double xmin, double xmax):
        self._xlim[0] = xmin
        self._xlim[1] = xmax

    ## set the range of y
    cpdef void setyrange(self, double ymin, double ymax):
        self._ylim[0] = ymin
        self._ylim[1] = ymax

    ## set the max temperature
    cpdef void setmaxtemp(self, double maxtemp):
        self._maxtemp = maxtemp

    ## set color bar limits
    cpdef void setcolorlim(self, double minT, double maxT):
        self._min = minT
        self._max = maxT

    ## fill the geometry info into patches
    cpdef void fillpatch(self, int phi=1, MatType opt=kCapacity, PosType dim=kZ):
        cdef np.ndarray geoinfo = self._geof.GetCollection()
        cdef np.ndarray matinfo = self._matf.GetCollection()
        cdef double z
        cdef double r
        cdef double lz
        cdef double lr

        if len(geoinfo)!=len(matinfo):
            print "size of geometry info is not equal to the size of material info"
            raise

        for i in range(len(geoinfo)): 
            if geoinfo[i].GetId(kPhi)==phi:
                z  = geoinfo[i].GetPosition(kZ)
                r  = geoinfo[i].GetPosition(kR)
                lz = geoinfo[i].GetWidth(kZ)
                lr = geoinfo[i].GetWidth(kR)
                self._patch = np.append( self._patch, Rectangle((z,r), lz, lr, linewidth=0.1) )
                if opt==kCapacity:
                    self._dset = np.append( self._dset, matinfo[i].GetCapacity() )
                elif opt==kRRR:
                    self._dset = np.append( self._dset, matinfo[i].GetRRR() )
                elif opt==kConductivity:
                    self._dset = np.append( self._dset, matinfo[i].GetConductivity(dim) )
                elif opt==kField:
                    self._dset = np.append( self._dset, matinfo[i].GetField() )
                elif opt==kTemperature:
                    self._dset = np.append( self._dset, matinfo[i].GetTemperature() )
                elif opt==kDose:
                    self._dset = np.append( self._dset, matinfo[i].GetDose() )
                else:
                    raise


    ## set one piece along the phi direction
    cpdef SetPhi(self, int phi=1):
        self._phi = phi


    ## set the direction
    cpdef SetDirection(self, PosType dim=kZ):
        self._direct = dim


    ## set plot info
    def SetMatInfo(self, MatType info=kConductivity):
        self._info = info
        if info==kTemperature:
            self._name = "Temperature [K]"
            self._save = "temperature"
        elif info==kConductivity:
            self._name = "Thermal Conductivity [W/m/K]"
            self._save = "conductivity"
        elif info==kCapacity:
            self._name = "Heat Capacity [J/kg/K]"
            self._save = "capacity"
        elif info==kDose:
            self._name = "Heat Generation [W/m$^{3}$]"
            self._save = "dose"
        elif info==kField:
            self._name = "Magnetic Field [Tesla]"
            self._save = "field"
        elif info==kRRR:
            self._name = "RRR"
            self._save = "rrr"
        else:
            print "Warning: there is no this kind of material property."


    ## plot geometry with given info
    def Draw(self, path=""):
        self.fillpatch(self._phi, self._info, self._direct)
        p = PatchCollection(self._patch, cmap=matplotlib.cm.jet, linewidth=0.1)
        #p = PatchCollection(self._patch, cmap=matplotlib.cm.inferno, linewidth=0.1)
        p.set_array(self._dset)
        
        fig, ax = plt.subplots(1, 1, figsize=(14,5.8))
        ax.add_collection(p)
        if self._xlim[0]==0. and self._xlim[1]==0.:
            ax.set_xlim( self._geof.CheckLimits(kZ) )
        else:
            ax.set_xlim( [self._xlim[0], self._xlim[1]] )

        if self._ylim[0]==0. and self._ylim[1]==0.:
            ax.set_ylim( self._geof.CheckLimits(kR) )
        else:
            ax.set_ylim( [self._ylim[0], self._ylim[1]] )

        if self._maxtemp!=0. and self._info==kTemperature:
            ax.set_title("Max. temperature = %.3f [K]" %self._maxtemp)

        if self._max!=0.:
            p.set_clim(self._min, self._max)
        ax.set_xlabel( "Z [mm]" , fontsize=16)
        ax.set_ylabel( "R [mm]" , fontsize=16)
        ax.tick_params(axis="x", labelsize=16)
        ax.tick_params(axis="y", labelsize=16)
        cbar = plt.colorbar(p)
        cbar.ax.tick_params(labelsize=16)
        cbar.set_label( self._name, fontsize=16)
        plt.tight_layout()
        fig.savefig(path+self._save+".pdf", bbox_inches="tight")
        #plt.show()

