#!/usr/bin/python3.6
import scipy.sparse as sp
from gmshtools import readGmshFile
import numpy as np

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
                while (line != '$EndElements\n'):
                    intLine = tuple(map(int,line.split()))
                    objNodes = nodesAsign(intLine[5:len(intLine)],cNodes)
                    entityList.append(Elem(objNodes,'Q4'))
                    line = entityFile.readline()  
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
        conectivity.append(Elem(objNodes,'Q9')) #FIXME : corregir que no asigne solo a Q4
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

class ElemQ4:
    def __init__(self,nodes):
        self.nloc = nodes
        self.H = self.buildMatrix()

    def localNodes(self):
        initList = []
        for i in range(len(self.nloc)):
            initList.append(self.nloc[i].nglob)
        return initList

    def Cmat(self):
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

    def fundHrs(self,r,s):
        hr=[1+s , -1-s , -1 +s ,1-s]
        hs=[1+r,1-r , -1 +r , -1 -r]
        dhrs = np.matrix(np.zeros((2,4)))
        dhrs[0] = hr
        dhrs[1] = hs
        return dhrs*0.25

    def funB(self,dHrs,J):
        B = np.matrix(np.zeros((3,8)))
        Jinv = np.linalg.inv(J)
        DH = Jinv * dHrs
        for i in range(2):
            B[i,i::2] = DH[i]
            B[2,i::2] = DH[i-1]
        return B

    def getpos(self):
        auxlist = []
        for i in range(4):
            auxlist.append([self.nloc[i].x,self.nloc[i].y])
        return np.matrix(auxlist)   

    def buildMatrix(self):
        gpoints = [-0.5773, 0.5773]
        He = []
        for gpr in gpoints:
            for gps in gpoints:
                dHrs = self.fundHrs(gpr,gps)
                He.append(dHrs)
        return He
    
    def calcJacobian(self,gps):
        return self.H[gps] * self.getpos()

    def getKe(self,C):
        J = np.matrix(np.zeros((2,2)))
        Ke = np.matrix(np.zeros((8,8)))
        Be = np.matrix(np.zeros((3,8)))
        for i in range(4):
            J = self.calcJacobian(i)
            detJ = np.linalg.det(J)
            B = self.funB(self.H[i],J)
            Ke += B.T * C * B * detJ *10
            Be += B * detJ *10
        self.Bstress = Be
        return Ke

    def getLocalDisplacement(self):
        U = np.matrix(np.zeros((4,2)))
        count = 0
        for node in self.nloc:
            U[count , 0] = node.xValue
            U[count,1] = node.yValue
            count+=1
        U = U.reshape((8,1))
        return U

    def computeStress(self,C):
        dH = self.fundHrs(0,0)
        J = dH * self.getpos()
        B = self.funB(dH,J)
        U = self.getLocalDisplacement()
        Stress = C * B * U
        return Stress

    def elemStressToNodes(self,stress):
        for localNode in self.nloc:
            localNode.storeStress(stress)
        return 0

class Elem:
    def __init__(self,nodes,elemType):
        self.elemType = elemType
        self.nloc = nodes

    def localNodes(self):
        initList = []
        for i in range(len(self.nloc)):
            initList.append(self.nloc[i].nglob)
        return initList

    def Cmat(self):
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


    def funB(self,dHrs,J):
        B = np.matrix(np.zeros((3,18)))
        Jinv = np.linalg.inv(J)
        Dh = Jinv * dHrs
        for i in range(2):
            B[i,i::2] = Dh[i]
            B[2,i::2] = Dh[i-1]
        return B

    def getpos(self):
        auxlist = []
        for i in range(9):
            auxlist.append([self.nloc[i].x,self.nloc[i].y])
        return np.matrix(auxlist)   


    def calcJacobian(self,gps,Hrs):
        pos = self.getpos()
        return Hrs[gps] * pos
    # @profile
    def getKe(self,H,Hrs,C):
        J = np.matrix(np.zeros((2,2)))
        Ke = np.matrix(np.zeros((18,18)))
        Be = np.matrix(np.zeros((3,18)))
        for i in range(9):
            J = self.calcJacobian(i,Hrs)
            detJ = np.linalg.det(J)
            B = self.funB(Hrs[i],J)
            Ke += B.T * C * B * detJ *10
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
        dH = Hrs
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
    def __init__(self,meshName,nelem,nnode,elemType,conectivity,bcNodes):
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
            Exception -- [description]

        Returns:
            K[lil_matrix] -- Empty Global rigidity matrix
            brhs[lil_matrix] -- Empty Right hand side vector
            H[array] -- Array with form functions evaluated in gps
            Hrs[array] -- It haves the derivatives in gps
        """

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
            gpsList = [-0.5777,0.5777] 
            self.H = self.quadH(gpsList,self.funH4)
            self.Hrs = self.quadHrs(gpsList,self.funHrs4)

        else:
            raise Exception('Invalid elemType defined you must use Quad4 or Quad9')


    def quadH(self,gpoints,formFunction,gpweights=None):
        H = []
        for gpr in gpoints:
            for gps in gpoints:
                dH = formFunction(gpr,gps)
                H.append(dH)
        return H

    def quadHrs(self,gpoints,formFunction,gpweights=None):
        He = []
        for weir , gpr in enumerate(gpoints):
            for weis, gps in enumerate(gpoints):
                gprWei = gpweights[weir]
                gpsWei = gpweights[weis]
                dHrs = formFunction(gpr,gps)*gprWei*gpsWei
                He.append(dHrs)
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
        r_power = 1- r**2
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
        hrquad8 = [-r*(1+s),-0.5*(1-s**2),-r*(1-s),0.5*(1-s**2)]
        hrquad9 = -2*r*(1-s**2)
        hsquad4 = [1+r,1-r , -1 +r , -1 -r]
        hsquad8 = [0.5*(1-r**2),-s*(1-r),-0.5*(1-r**2),-s*(1+r)]
        hsquad9 = -2*r*(1-r**2)

        hr = [0]*9
        hs = [0]*9
        for i in range(4):
            hr[i] = 0.25*hrquad4[i] - 0.25*hrquad8[i] - 0.25*hrquad8[i-1] - 0.25*hrquad9
            hs[i] = 0.25*hsquad4[i] - 0.25*hsquad8[i] - 0.25*hsquad8[i-1] - 0.25*hsquad9

        for i in range(4):
            hr[i+4] = 0.5*hrquad8[i] - 0.5*hrquad9
            hs[i+4] = 0.5*hsquad8[i] - 0.5*hsquad9

        hr[8] = hrquad9
        hs[8] = hsquad9

        dhrs = np.matrix(np.zeros((2,9)))

        dhrs[0] = hr
        dhrs[1] = hs
        return dhrs    

    def getRowData(self,elem,Ke):
        '''
        Optimized one
        '''
        listedNodes = [(x-1)*2+y for x in elem.localNodes() for y in [0,1]]
        col = listedNodes*18
        row = [ind for ind in listedNodes for i in range(18)]
        Kelem = Ke.ravel()
        return row , col , Kelem

    def assemble(self, C):
        #Armado de matrices elementales y ensamblaje
        row = []
        col = []
        Kelem = []
        forceX = 50
        forceY = 0
        for elem in self.conectivity:
            Ke = elem.getKe(self.H,self.Hrs,C)
            #self.K , self.brhs = Assemble(elem, Ke , self.K , self.brhs)
            rowap , colap , Kelemap = self.getRowData(elem,Ke)
            row.extend(rowap)
            col.extend(colap)
            Kelem.append(Kelemap)

        #Setup de matrices esparsas
        Kelem = np.array(Kelem).ravel()
        self.K = sp.coo_matrix((Kelem,(row,col)),shape=(self.nnode*2,self.nnode*2)).tolil()

        indexList = []

        #FIXME : This has to be done only in DIR nodes!
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
        return None

    def setBoundaryConditions(self, coords):
            '''
            it only receives a tuple with nodes and assing only forces and construct brhs from element
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

            return None

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

    def setNEUx(self,bc):
        self.NEU = True
        self.xForce = bc

    def setNEUy(self,bc):
        self.yForce = bc

    def setDIRxValue(self,value):
        self.DIRx = True
        self.xValue = value

    def setDIRyValue(self,value):
        self.DIRy = True
        self.yValue = value

    def storeCalcValue(self,xValue,yValue):
        self.xValue = xValue
        self.yValue = yValue

    def storeStress(self,stressVector):
        for i in range(3):
            self.stress[i] = (stressVector[0,i] + self.stress[i]) / (self.markStress+1) - self.stress[i]/(self.markStress)
        self.markStress += 1

    def physGrouptoValue(self,bclist):
        '''
        Pasar los numeros del phys group a un valor
        '''
        group = self.bcGroup - 1
        pairValue = bclist[group]
        aux = 0
        for bc in pairValue:
            if group == 0:
                if aux ==0:
                    self.setDIRxValue(bc)
                else:
                    self.setDIRyValue(bc)
            else:
                if aux ==0:
                    self.setNEUx(bc)                  
                else:
                    self.setNEUy(bc)
            aux += 1

def funHrs(r,s):
    #FIXME: FUNCION NO IMPLEMENTADA, CONTIENE LAS FUNCIONES DE INTERPOLACION
    pass

# @profile


