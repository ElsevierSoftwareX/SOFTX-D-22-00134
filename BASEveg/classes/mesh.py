# Implementation of classes
# F. Caponi, 2020 (modified from S.J. Peter, 2015)

import numpy as np

class NODE:
    def __init__ (self, x=0.0, y=0.0, z=0.0):
        self.__coords = [x, y, z]
        self.__index = None
        self.__attribute = None

    def __getitem__ (self,ind):
        return self.__coords[ind]
    
    def setID (self, index):
        self.__index = index
    
    def getID (self):
        return self.__index

    def getX (self):
        return self.__coords[0]

    def getY (self):
        return self.__coords[1]

    def getZ (self):
        return self.__coords[2]
    
    def getPointXY(self):
        return [self.__coords[0], self.__coords[1]]
        
    def getPointXYZ(self):
        return [self.__coords[0], self.__coords[1], self.__coords[2]]
    
    def setAttribute(self, attribute):
        self.__attribute = attribute

    def getAttribute(self):
        return self.__attribute

class CELL:

    def __init__ (self, node1, node2, node3):
        self.__nodes = [node1,node2,node3]
        self.__matID = None
        self._center = [0.0,0.0,0.0]
        self._wse_list = []
        self._z_list = 0.0
        self._time_list = []
        self._pol = np.poly1d([1,2,3])
        self._freq_wse = []
        self._freq_bin_wse = []
        self._fit = True
        self._pdf_parameters = [-1.,-1.]
        self._perr = -1.
        self._Bc = 0.0
        self._Br = 0.0
        self._B = 0.0
        self._D = 0.0
        self._H = 0.0
        self._upr = self._D
        self._burial = self._H
        self._root_dist = []
        self._sigmaB = 0.0
        self._WSE = 0.0
        self._waterTable = []
        self._isgrowing = False
        self._netgrowth = 0
        self._netdecay = 0
        self._B0 = 0.0
        self._D0 = 0.0
        #print("Cell obj created at nodes")
        #print(self.__nodes)

# set methods
    def setWseList(self,value):
        self._wse_list = value
    def setTimeList(self,value):
        self._time_list = value
    def setPol(self,value):
        self._pol = value
    def setFit(self,value):
        self._fit = value
    def setParameters(self,values):
        self._pdf_parameters = values
    def setB(self,value):
        self._B = value
    def setBc(self,value):
        self._Bc = value
    def setBr(self,value):
        self._Br = value
    def setD(self,value):
        self._D = value
    def setH(self,value):
        self._H = value
    def setUpr(self,value):
        self._upr = value
    def setBurial(self,value):
        self._burial = value
    def setRootDist(self,value):
        self._root_dist = value
    def setSigmaB(self,value):
        self._sigmaB = value
    def setPerr(self,value):
        self._perr = value
    def setWSE(self,value):
        self._WSE = value 
    def setWaterTable(self,value):
        self._waterTable = value 
    def setIsgrowing(self,value):
        self._isgrowing = value
    def setNETgrowth(self,value):
        self._netgrowth = value 
    def setNETdecay(self,value):
        self._netdecay = value 
    def setB0(self,value):
        self._B0 = value 
    def setD0(self,value):
        self._D0 = value

# get methods
    def getTimeList(self):
        return self._time_list
    def getWseList(self):
        return self._wse_list
    def getPol(self):
        return self._pol 
    def getFit(self):
        return self._fit
    def getParameters(self):
        return self._pdf_parameters
    def getBc(self):
        return self._Bc
    def getBr(self):
        return self._Br
    def getB(self):
        return self._B
    def getD(self):
        return self._D
    def getH(self):
        return self._H
    def getUpr(self):
        return self._upr
    def getBurial(self):
        return self._burial
    def getRootDist(self):
        return self._root_dist
    def getSigmaB(self):
        return self._sigmaB
    def getPerr(self):
        return self._perr
    def getWaterTable(self):
        return self._waterTable
    def getIsgrowing(self):
        return self._isgrowing
    def getNETgrowth(self):
        return self._netgrowth
    def getNETdecay(self):
        return self._netdecay        
    def getB0(self):
        return self._B0
    def getD0(self):
        return self._D0 
# add methods
    def addTime(self,value):
        self._time_list.append(value)
    def addWSE(self,value):
        self._wse_list.append(value)

    def setFrequency(self,value,bin_edg):
        self._freq_wse = value
        self._freq_bin_wse = bin_edg

    def getFrequency(self):
        return [self._freq_wse, self._freq_bin_wse]

    def _init_veg(self,B,D,Bc,Br):
        self.setB( B )
        self.setD( D )
        self.setBc( Bc )
        self.setBr( Br )
        #print "Init of elements done!"

    def calcProperties (self,center): 
        nodes = self.getNodes()
        # calculate the element center coordinates
        x1 = 0.0
        y1 = 0.0
        N = len(nodes)
        if center=='average':
            # average of triangle nodes coordinates as interpolation point
            for nn in nodes:
                x1 += nn.getX()
                y1 += nn.getY()
            x = x1 / N
            y = y1 / N
        elif center=='outer':
            # outer circle center as interpolation point
            d =  0.0
            for ii in range(N):
                n1 = nodes[(ii-1)%N]
                n2 = nodes[ii]
                n3 = nodes[(ii+1)%N]
                dx = n3.getX()-n1.getX()
                dy = n3.getY()-n1.getY()
                d = d + n2.getX()*dy
                v = n2.getX()**2 + n2.getY()**2
                x = x + v*dy
                y = y - v*dx
            x = x / (2.0*d)
            y = y / (2.0*d)
        elif center=='inner':
            # inner circle cener as interpolation point
            lsum = 0.0
            for ii in range(N):
                n1 = nodes[(ii+1)%N]
                n2 = nodes[(ii+2)%N]
                l = math.sqrt((n2.getX()-n1.getX())**2 + (n2.getY()-n1.getY())**2)
                lsum += l
                x = x + l*nodes[ii].getX()
                y = y + l*nodes[ii].getY()
            x = x/lsum
            y = y/lsum
        elif center=='mid':
            # mid of xmax and xmin (ymax and ymin)
            xmax = max(nodes[0].getX(),max(nodes[1].getX(),nodes[2].getX()))
            xmin = min(nodes[0].getX(),min(nodes[1].getX(),nodes[2].getX()))
            ymax = max(nodes[0].getY(),max(nodes[1].getY(),nodes[2].getY()))
            ymin = min(nodes[0].getY(),min(nodes[1].getY(),nodes[2].getY()))
            x = 0.5*(xmin+xmax)
            y = 0.5*(ymin+ymax)

        self._center = [x, y, 0.0]

    def getCenter (self):
        return self._center

    def setElevation (self, z):
        self._center[2] = z

    def getElevation (self):
        return self._center[2]
    
    #def getElementCenterCoords (self, elID):
     #   return self._elements[elID].getCenter()		
    
    def getNodes (self):
        return self.__nodes

    def getNodeIDs (self):
        nids = []
        for ni in self.__nodes:
            nids.append(ni.getID())
        return nids

    def setMaterial (self, value):
        self.__matID = value

    def getMaterial (self):
        return self.__matID

    def calcArea (self):
        n1,n2,n3 = self.getNodes()
        return abs((n1.getX()*(n2.getY()-n3.getY()) + n2.getX()*(n3.getY()-n1.getY())+ n3.getX()*(n1.getY()-n2.getY()))/2.0)

    def isInside (self, x, y):
        whichSide = lambda pt, n1, n2: True if (pt[0] - n2.getX()) * (n1.getY() - n2.getY()) - (pt[1] - n2.getY())*(n1.getX() - n2.getX()) >= 0.0 else False
        side1 = whichSide((x, y), self.__nodes[0], self.__nodes[1])
        N = len(self.__nodes)
        for ii in range(1, N):
            side2 = whichSide((x, y), self.__nodes[ii], self.__nodes[(ii+1)%N])
            if side2 != side1:
                return False
            else:
                side1 = side2
        return True

    def interpolate (self, x, y):
        P = self.__nodes[0]
        Q = self.__nodes[1]
        R = self.__nodes[2]
        a =  -(R.getY()*Q.getZ() - P.getY()*Q.getZ() - R.getY()*P.getZ() + Q.getY()*P.getZ() + P.getY()*R.getZ() - Q.getY()*R.getZ())
        b =   P.getY()*R.getX() + Q.getY()*P.getX() + R.getY()*Q.getX() - Q.getY()*R.getX() - P.getY()*Q.getX() - R.getY()*P.getX()
        c =   Q.getZ()*R.getX() + P.getZ()*Q.getX() + R.getZ()*P.getX() - P.getZ()*R.getX() - Q.getZ()*P.getX() - Q.getX()*R.getZ()
        d = -a*P.getX() - b*P.getZ() - c*P.getY()
        if b == 0.0:
            if ((x == P.getX()) and (y == P.getY)):
                return P.getZ()
            if (x == Q.getX() and y == Q.getY()):
                return Q.getZ()
            if (x == R.getX() and y == R.getY()):
                return R.getZ()
            return 0.0
        return -(a*x+c*y+d) / b