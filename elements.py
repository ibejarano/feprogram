#!/usr/bin/python3.6
import scipy.sparse as sp
from gmshtools import readGmshFile
import numpy as np
import math

def Constructor(entity,inpGmsh,cNodes,bc=False):
    '''
    entity: Entidad para extraer datos de malla, puede ser $Nodes o $Elements
    inpGmsh: Input archivo .gmsh para extraer datos
    cNodes: Lista con Nodos objetos
    nEntity: Cantidad de entidades encontradas
    entityList: Lista con entidades
    '''
    nEntity , entityFile = readGmshFile(entity,inpGmsh)
    entityList = []
    if entity == '$Nodes':
        for items in range(nEntity):
            line = entityFile.readline()  
            intLine = tuple(map(float,line.split()))
            entityList.append(Node2D(intLine[1:4]))
    else:
        line = entityFile.readline()  
        intLine = tuple(map(int,line.split()))
        counter = 0
        if bc:
            while (len(intLine) < 9):
                nod = intLine[5:len(intLine)]
                physGroup = intLine[3]
                for intNode in nod:
                    objNode = cNodes[intNode-1]
                    objNode.setBG(physGroup)
                line = entityFile.readline()
                intLine = tuple(map(int,line.split()))
                entityList.append(nod)
                counter +=1
            return nEntity , entityList
        else:
            if len(intLine) == 8:
                entityList = readElements(entityFile,entityList,'Q9',cNodes)

            else:
                entityList = readElements(entityFile,entityList,'Q4',cNodes)
            nEntity -=1
    return nEntity  , entityList

def readElements(fileGmsh,conectivity,elemType,coords):
    '''
    Once it identifies the elemType it gets here and iterate
    '''
    line = ' '
    intLine = []
    while (len(intLine) < 9):
        line = fileGmsh.readline() 
        intLine = tuple(map(int,line.split()))

    while (line != '$EndElements\n'):
        intLine = tuple(map(int,line.split()))
        objNodes = nodesAsign(intLine[5:len(intLine)],coords)
        conectivity.append(Elem(objNodes,elemType))
        line = fileGmsh.readline() 

    return conectivity

def nodesAsign(tNodes,coordNodes):
    '''
    tNodes: Tupla con los nodos enteros
    tobj: Tupla con los objetos nodos
    '''
    objNod = []
    for i in tNodes:
        objNod.append(coordNodes[i-1])
    return objNod

class Elem:
    def __init__(self,nodes,elemType):
        self.elemType = elemType
        self.nloc = nodes

    def localNodes(self):
        initList = []
        for i in range(len(self.nloc)):
            initList.append(self.nloc[i].nglob)
        return initList

    def funB(self,dHrs,J):
        dof = len(self.nloc)*2
        B = np.matrix(np.zeros((3,dof)))
        Jinv = np.linalg.inv(J)
        Dh = Jinv * dHrs
        for i in range(2):
            B[i,i::2] = Dh[i]
            B[2,i::2] = Dh[i-1]
        return B

    def getpos(self):
        auxlist = []
        nNodeslocal = len(self.nloc)
        for i in range(nNodeslocal):
            auxlist.append([self.nloc[i].x,self.nloc[i].y])
        return np.matrix(auxlist)   


    def calcJacobian(self,Hrs):
        pos = self.getpos()
        return Hrs * pos
    # @profile
    def getKe(self,H,Hrs,C,gpWei):
        dof = len(self.nloc)*2
        J = np.matrix(np.zeros((2,2)))
        Ke = np.matrix(np.zeros((dof,dof)))
        Be = np.matrix(np.zeros((3,dof)))
        for i in range(int(dof/2)):
            J = self.calcJacobian(Hrs[i])
            detJ = np.linalg.det(J)
            B = self.funB(Hrs[i],J)
            Ke += B.T * C * B * detJ*gpWei[i]*10
            Be += B * detJ *10
        self.Bstress = Be
        return Ke

    def getLocalDisplacement(self):
        numNodos = len(self.nloc)
        U = np.matrix(np.zeros((numNodos,2)))
        count = 0
        for node in self.nloc:
            U[count , 0] = node.xValue
            U[count,1] = node.yValue
            count+=1
        U = U.reshape((numNodos*2,1))
        return U

    def computeStress(self,C,Hrs):
        dH = Hrs[-1]
        J = dH * self.getpos()
        B = self.funB(dH,J)
        U = self.getLocalDisplacement()
        Stress = C * B * U
        return Stress

    def elemStressToNodes(self,stress):
        for localNode in self.nloc:
            localNode.storeStress(stress)
        return 0

class FemProblem:
    def __init__(self,nelem,nnode,elemType,conectivity,bcNodes):
        '''Instance the FEM problem, the code get around this class
        
        Arguments:
            meshName {string} -- name of mesh, DEPRECATED
            nelem {integer} -- total element number of problem
            nnode {integer} -- total node number of problem
            elemType {string} -- declares if type is quad9 or quad4
            conectivity {array} -- contains objects class Elem
            bcNodes {array} -- contains a list with nodes with boundary conditions
        '''

        self.nelem = nelem
        self.nnode = nnode
        self.elemType = elemType
        self.conectivity = conectivity
        self.bcNodes = bcNodes

    def setMatrix(self):
        """[summary]
            This function setups Empty K , brhs in all cases
            H & Hrs depends on ElemType

        Raises:
            Exception -- if no elemType defined nor implemented yet

        Returns:
            K[lil_matrix] -- Empty Global rigidity matrix
            brhs[lil_matrix] -- Empty Right hand side vector
            H[array] -- Array with form functions evaluated in gps
            Hrs[array] -- It haves the derivatives in gps
        """
        self.C = self.Cmat()
        ncols = (self.nnode)*2
        nrows = (self.nnode)*2
        self.K = sp.lil_matrix((nrows,ncols))
        self.brhs = sp.lil_matrix((nrows,1))

        if self.elemType == 'Quad9':
            gpsList = [-0.774,0, 0.774]
            gpsWei = [0.55555 , 0.88888 , 0.55555]
            self.H = self.quadH(gpsList,self.funH9,gpsWei)
            self.Hrs = self.quadHrs(gpsList,self.funHrs9,gpsWei)
        
        elif self.elemType == 'Quad4':
            gpsList = [-0.5773,0.5773]
            gpsWei = [1,1]
            self.H = self.quadH(gpsList,self.funH4,gpsWei)
            self.Hrs = self.quadHrs(gpsList,self.funHrs4,gpsWei)

        else:
            raise Exception('Invalid elemType defined you must use Quad4 or Quad9')

    def quadH(self,gpoints,formFunction,gpweights=None):
        '''Build the H matrix with the formfunctions and evaluates at gauss points and weights if its needed
        
        Arguments:
            gpoints {list} -- list containing iterable points
            formFunction {callback} -- function that returns the value of form functions at gauss point
        
        Keyword Arguments:
            gpweights {list} -- weights of numerical integration, only needed when elemType=Quad9 (default: {None})
        
        Returns:
            array -- The columns are hi form function and row are gauss points number
        '''

        H = []
        for gpr in gpoints:
            for gps in gpoints:
                dH = formFunction(gpr,gps)
                H.append(dH)
        return H

    def quadHrs(self,gpoints,formFunction,gpweights=None):
        '''Creates an array that has the derivatives of form
        function evaluated at gauss points
        
        Arguments:
            gpoints {list[float]} -- List with gauss points
            formFunction {callback} -- returns the gp evaluated
        
        Keyword Arguments:
            gpweights {list[float]} -- it weight numerical integration (default: {None})
        
        Returns:
            He [array] -- contains the functions evaluated at gausspoints and ordered
        '''

        He = []
        gpWei = []
        for weir , gpr in enumerate(gpoints):
            for weis, gps in enumerate(gpoints):
                gprWei = gpweights[weir]
                gpsWei = gpweights[weis]
                dHrs = formFunction(gpr,gps)
                gpWei.append(gprWei*gpsWei)
                He.append(dHrs)
        self.gpWei = gpWei
        return He

    def funH4(self,r,s):
        r_plus = 1+r
        r_minus = 1-r
        s_plus = 1+s
        s_minus = 1-s

        h4 = 0.25*r_plus*s_minus
        h3 = 0.25*r_minus*s_minus
        h2 = 0.25*r_minus*s_plus
        h1 = 0.25*r_plus*s_plus

        H=np.matrix([h1,h2,h3,h4])
        return H

    def funH9(self,r,s):
        r_power = 1-r**2
        s_power = 1-s**2
        r_plus = 1+r
        r_minus = 1-r
        s_plus = 1+s
        s_minus = 1-s

        h9 = r_power*s_power
        h8 = 0.5*s_power*r_plus - 0.5*h9
        h7 = 0.5*r_power*s_minus - 0.5*h9
        h6 = 0.5*s_power*r_minus - 0.5*h9
        h5 = 0.5*r_power*s_plus - 0.5*h9
        h4 = 0.25*r_plus*s_minus - 0.5*h7 - 0.5*h8 -0.25*h9
        h3 = 0.25*r_minus*s_minus - 0.5*h6 - 0.5*h7 - 0.25*h9
        h2 = 0.25*r_minus*s_plus - 0.5*h5 - 0.5*h6 - 0.25*h9
        h1 = 0.25*r_plus*s_plus - 0.5*h5 -0.5*h8 -0.25*h9

        H=np.matrix([h1,h2,h3,h4,h5,h6,h7,h8,h9])

        return H

    def funHrs4(self,r,s):
        hrquad4 = [1+s , -1-s , -1+s ,1-s]
        hsquad4 = [1+r,1-r , -1 +r , -1 -r]

        hr = [0]*4
        hs = [0]*4
        for i in range(4):
            hr[i] = 0.25*hrquad4[i]
            hs[i] = 0.25*hsquad4[i]

        dhrs = np.matrix(np.zeros((2,4)))

        dhrs[0] = hr
        dhrs[1] = hs
        return dhrs    

    def funHrs9(self,r,s):
        hrquad4 = [1+s , -1-s , -1 +s ,1-s]
        hrquad8 = [-2*r*(1+s),-1*(1-s**2),-2*r*(1-s),(1-s**2)]
        hrquad9 = -2*r*(1-s**2)
        hsquad4 = [1+r,1-r ,-1+r ,-1-r]
        hsquad8 = [(1-r**2),-2*s*(1-r),-1*(1-r**2),-2*s*(1+r)]
        hsquad9 = -2*s*(1-r**2)

        hr = [0]*9
        hs = [0]*9

        for i in range(4):
            hr[i+4] = 0.5*hrquad8[i] - 0.5*hrquad9
            hs[i+4] = 0.5*hsquad8[i] - 0.5*hsquad9

        hr[8] = hrquad9
        hs[8] = hsquad9
        
        for i in range(1,4):
            hr[i] = 0.25*hrquad4[i] - 0.5*hr[i+4] - 0.5*hr[i+3] - 0.25*hrquad9
            hs[i] = 0.25*hsquad4[i] - 0.5*hs[i+4] - 0.5*hs[i+3] - 0.25*hsquad9
            
        hr[0] = 0.25*hrquad4[0] - 0.5*hr[4] - 0.5*hr[7] - 0.25*hrquad9
        hs[0] = 0.25*hsquad4[0] - 0.5*hs[4] - 0.5*hs[7] - 0.25*hsquad9
        
        dhrs = np.matrix(np.zeros((2,9)))

        dhrs[0] = hr
        dhrs[1] = hs
        return dhrs    

    def getRowData(self,elem,Ke):
        '''
        Optimized one
        '''
        dof = len(elem.nloc)*2
        listedNodes = [(x-1)*2+y for x in elem.localNodes() for y in [0,1]]
        col = listedNodes*dof
        row = [ind for ind in listedNodes for i in range(dof)]
        Kelem = Ke.ravel()
        return row , col , Kelem

    def assemble(self):
        #Armado de matrices elementales y ensamblaje
        row = []
        col = []
        Kelem = []
        forceX = 1
        forceY = 0
        for elem in self.conectivity:
            Ke = elem.getKe(self.H,self.Hrs,self.C,self.gpWei)
            #self.K , self.brhs = Assemble(elem, Ke , self.K , self.brhs)
            rowap , colap , Kelemap = self.getRowData(elem,Ke)
            row.extend(rowap)
            col.extend(colap)
            Kelem.append(Kelemap)

        #Setup de matrices esparsas
        Kelem = np.array(Kelem).ravel()
        self.K = sp.coo_matrix((Kelem,(row,col)),shape=(self.nnode*2,self.nnode*2)).tolil()

        indexList = []

        for dirInd in self.nodesDIR:
            i = (dirInd-1)*2
            k = i+1
            _, J = self.K[i,:].nonzero()
            _, T = self.K[k,:].nonzero()
            for j in J:
                if j != i:
                    indexList.append(j)
                    self.K[i,j] = 0.0
            for t in T:
                if t != k:
                    self.K[k,t] = 0.0
            
            self.K[i,i] = 1.0
            self.K[k,k] = 1.0

        for nodeForce in self.nodesForce:
            glx = (nodeForce - 1)*2
            gly = glx+1
            self.brhs[glx] = forceX
            self.brhs[gly] = forceY

        self.K = self.K.tocsc()
        self.brhs = self.brhs.tocsc()

    def setBoundaryConditions(self, coords):
            '''
            it only receives a tuple with nodes and assign only forces and construct brhs from element
            '''
            #the first One is Dirichlet and its meant to be K = 1 in that index
            #if NEU brhs equals force
            nodesDirichlet = []
            nodesForce = []
            for node in self.bcNodes:
                if coords[node-1].bcGroup == 1:
                    nodesDirichlet.append(node)
                
                else:
                    nodesForce.append(node)

            self.nodesDIR = nodesDirichlet
            self.nodesForce = nodesForce

    def Cmat(self):
        '''Calculates relation between stress-strains
        
        Returns:
            C (array) -- matrix with relation
        '''

        nu = 0.3
        E = 200
        mat = np.zeros((3,3))
        auxlist = [[1,nu,0],[nu,1,0],[0,0,(1-nu)*0.5]]
        ind = 0
        for i in auxlist:
            mat[ind] = i
            ind +=1
        const = E / (1- nu**2)
        return mat*const

    def getStress(self):
        listStress = []
        C = self.C
        for elem in self.conectivity:
            localStress = elem.computeStress(C,self.Hrs)
            vonMises = math.sqrt(localStress[0]**2 + localStress[1]**2 - localStress[0]*localStress[1] + 3*localStress[2]**2)
            listStress.append(vonMises)
        return listStress

class Node2D:
    nglob = 1
    def __init__(self,coords):
        self.x = coords[0]
        self.y = coords[1]
        self.DIRx = False
        self.DIRy = False
        self.NEU = False
        self.nglob = self.nglob
        self.BG = False
        Node2D.nglob += 1
        self.stress = [0,0,0]
        self.markStress = 1

    def setBG(self, group):
        self.BG = True
        self.bcGroup = group

    def storeCalcValue(self,xValue,yValue):
        self.xValue = xValue
        self.yValue = yValue

    def storeStress(self,stressVector):
        for i in range(3):
            self.stress[i] = (stressVector[0,i] + self.stress[i]) / (self.markStress+1) - self.stress[i]/(self.markStress)
        self.markStress += 1

# @profile