#!/usr/bin/python3.6

from gmshtools import readGmshFile
import numpy as np

def Constructor(entity,inpGmsh,cNodes):
    '''
    entity: Entidad para extraer datos de malla, puede ser $Nodes o $Elements
    inpGmsh: Input archivo .gmsh para extraer datos
    cNodes: Lista con Nodos objetos
    nEntity: Cantidad de entidades encontradas
    entityList: Lista con entidades
    '''
    nEntity , entityFile = readGmshFile(entity,inpGmsh)
    entityList = []
    for items in range(nEntity):
        line = entityFile.readline()
        if entity == '$Nodes':
            intLine = tuple(map(float,line.split()))
            entityList.append(Node2D(intLine[1:4]))
        else:
            intLine = tuple(map(int,line.split()))
            if len(intLine) == 9:
                objNodes = nodesAsign(intLine[5:9],cNodes)
                entityList.append(ElemQ4(objNodes))
            else:
                physGroup = intLine[3]
                for nod in intLine[5:7]:
                    intNode = nod
                    objNode = cNodes[intNode-1]
                    objNode.setBG(physGroup)
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

    def Stress(self,U,C):
        '''
        Deprecated
        '''
        Pos = self.getpos()
        J = np.matrix(np.zeros((2,2)))
        Stre = np.matrix(np.zeros((3,1)))
        Hoo = self.fundHrs(0,0)
        J = Hoo*Pos
        B = self.funB(Hoo,J)
        Stre =  C * B *U *10
        return Stre.T

    def StressPost(self,Bel):
        self.CB = Bel

    def calcBaricenter(self,xList):
        xn = len(xList)
        baricenter = 0
        for x in xList:
            baricenter += x/xn
        return baricenter

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

def Assemble(elem,Ke,K,brhs):
    '''
    elem: Objeto elemento
    Ke: Matriz de rigidez local a ensamblar en global
    K: Matriz global
    gl: grado de libertad = nnodos x 2, es el indice de la K global
    '''
    for i in range(4):
        ni = elem.nloc[i].nglob-1
        gl = ni*2
        if elem.nloc[i].DIRx and elem.nloc[i].DIRy :
            K[gl,gl] = 1
            K[gl+1,gl+1] = 1
            brhs[gl] = elem.nloc[i].xValue
            brhs[gl+1] = elem.nloc[i].yValue
        else :
            if elem.nloc[i].NEU:
                brhs[gl] = elem.nloc[i].xForce
                brhs[gl+1] = elem.nloc[i].yForce               
            for j in range(4):
                nj = elem.nloc[j].nglob-1
                glj= nj*2
                K[gl:gl+2,glj:glj+2] += Ke[i*2:i*2+2,j*2:j*2+2]
    return K , brhs