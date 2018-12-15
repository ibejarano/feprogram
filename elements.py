#!/usr/bin/python3.6
import scipy.sparse as sp
from gmshtools import readGmshFile
import numpy as np

def searchElem(line,nElemLocal,listElements):

    pass

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
        if len(intLine) == 8:
            '''
            Si los elementos definidos en los bordes es len 8 entonces Q9
            '''
            elemType = 'Q9'

        else:
            elemType = 'Q4'
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
            if elemType == 'Q9':
                while (line != '$EndElements\n'):
                    intLine = tuple(map(int,line.split()))
                    objNodes = nodesAsign(intLine[5:len(intLine)],cNodes)
                    entityList.append(ElemQ9(objNodes)) #FIXME : corregir que no asigne solo a Q4
                    line = entityFile.readline()  
            else:
                while (line != '$EndElements'):
                    intLine = tuple(map(int,line.split()))
                    objNodes = nodesAsign(intLine[5:len(intLine)],cNodes)
                    entityList.append(ElemQ4(objNodes))
                    line = entityFile.readline()  
            nEntity -=1
    return nEntity  , entityList

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

class ElemQ9:
    def __init__(self,nodes):
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
        if self.elemType == 'Quad9':
            self.H = self.quad9H()
            self.Hrs = self.quad9Hrs()
            ncols = (self.nnode)*2
            nrows = (self.nnode)*2 + len(self.bcNodes)*2
            self.K = sp.lil_matrix((nrows,ncols))
            self.brhs = sp.lil_matrix((nrows,1))
        else:
            raise Exception('No se definió ningún tipo de elemento')

    def quad9H(self):
        gpoints = [-0.774,0, 0.774]
        H = []
        for gpr in gpoints:
            for gps in gpoints:
                dH = self.funH(gpr,gps)
                H.append(dH)
        return H

    def quad9Hrs(self):
        gpoints = [-0.774,0, 0.774]
        He = []
        for gpr in gpoints:
            for gps in gpoints:
                dHrs = self.fundHrs(gpr,gps)
                He.append(dHrs)
        return He

    def funH(self,r,s):
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

    def fundHrs(self,r,s):
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
        dirList = [1,1]
        for elem in self.conectivity:
            Ke = elem.getKe(self.H,self.Hrs,C)
            #self.K , self.brhs = Assemble(elem, Ke , self.K , self.brhs)
            rowap , colap , Kelemap = self.getRowData(elem,Ke)
            row.extend(rowap)
            col.extend(colap)
            Kelem.append(Kelemap)

        for dirNodes in self.bcNodes:
            bcrow = [(dirNodes-1)*2 , (dirNodes-1)*2 + 1 ]
            row.extend(bcrow)
            col.extend(bcrow)
            Kelem.append(dirList)
        #Setup de matrices esparsas
        Kelem = np.array(Kelem).ravel()
        bancaaK = sp.coo_matrix((Kelem,(row,col)),shape=((self.nnode+len(self.bcNodes))*2,self.nnode*2)).tocsc()
        self.brhs = self.brhs.tocsc()
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
    h1 = 0.25*(1+r)*(1+s)
    h2 = 0.25*(1-r)*(1+s)
    h3 = 0.25*(1-r)*(1-s)
    h4 = 0.25*(1+r)*(1-s)
    hrs = np.matrix([h1,h2,h3,h4])
    return hrs

# @profile
def Assemble(elem,Ke,K,brhs):
    '''
    elem: Objeto elemento
    Ke: Matriz de rigidez local a ensamblar en global
    K: Matriz global
    gl: grado de libertad = nnodos x 2, es el indice de la K global
    '''

    row = [(x-1)*2+y for x in elem.localNodes() for y in [0,1]]
    col = row
    for i in range(9):
        ni = elem.nloc[i].nglob-1
        gl = ni*2
        helprow = [ni,ni+1]
        row[i*2:i*2+2]=[ni,ni+1]
        if elem.nloc[i].DIRx and elem.nloc[i].DIRy :
            K[gl,gl] = 1
            K[gl+1,gl+1] = 1
            brhs[gl] = elem.nloc[i].xValue
            brhs[gl+1] = elem.nloc[i].yValue
        else :
            if elem.nloc[i].NEU:
                brhs[gl] = elem.nloc[i].xForce
                brhs[gl+1] = elem.nloc[i].yForce               
        for j in range(9):
            nj = elem.nloc[j].nglob-1
            glj= nj*2
            K[gl:gl+2,glj:glj+2] += Ke[i*2:i*2+2,j*2:j*2+2]
    return K , brhs

