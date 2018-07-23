import numpy as np
import scipy.sparse as sp
import sys
import logging
from gmshtools import readGmshFile
from scipy.sparse.linalg import spsolve

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
        intLine = tuple(map(int,line.split()))
        if entity == '$Nodes':
            entityList.append(Node2D(intLine[1:4]))
        elif len(intLine) == 9:
            objNodes = nodesAsign(intLine[5:9],cNodes)
            entityList.append(ElemQ4(objNodes))
        else:
            physGroup = intLine[3]
            intNode = intLine[6]
            objNode = cNodes[intNode-1]
            objNode.setBG(physGroup)
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

    def getpos(self):
        auxlist = []
        for i in range(4):
            auxlist.append([self.nloc[i].x,self.nloc[i].y])
        return np.matrix(auxlist)

    def fundHrs(self,r,s):
        hr=[0.25*y*x*(1+x*s) for x in [1,-1] for y in [1,-1]]
        hs=[0.25*x*(1+y*x*s) for x in [1,-1] for y in [1,-1]]
        dhrs = np.matrix(np.zeros((2,4)))
        dhrs[0] = hr
        dhrs[1] = hs
        return dhrs

    def funB(self,dHrs,J):
        B = np.matrix(np.zeros((3,8)))
        Jinv = np.linalg.inv(J)
        DH = Jinv * dHrs
        for i in range(2):
            B[i,i::2] = DH[i]
            B[1-i::2,::2] = DH[i]
        return B

    def getKe(self):
        Pos = self.getpos()
        gpoints = [-0.5773, 0.5773]
        J = np.matrix(np.zeros((2,2)))
        Ke = np.matrix(np.zeros((8,8)))
        for gpr in gpoints:
            for gps in gpoints:
                dHrs = self.fundHrs(gpr,gps)
                J = dHrs*Pos
                B = self.funB(dHrs,J)
                detJ = np.linalg.det(J)
                Ke += B.T * B * detJ
        return Ke

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

    def setBG(self, group):
        self.BG = True
        self.bcGroup = group

    def setNEU(self):
        self.NEU = True

    def setDIRxValue(self,value):
        self.DIRx = True
        self.xValue = value

    def setDIRyValue(self,value):
        self.DIRy = True
        self.yValue = value

    def storeCalcValue(self,xValue,yValue):
        self.xValue = xValue
        self.yValue = yValue

    def physGrouptoValue(self,bclist):
        '''
        Pasar los numeros del phys group a un valor
        '''
        group = self.bcGroup - 1
        pairValue = bclist[group]
        aux = 0
        for bc in pairValue:
            if bc != None:
                if aux == 0:
                    self.setDIRxValue(bc)
                else:
                    self.setDIRyValue(bc)
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
            K[gl,gl:gl+2] = 1
            brhs[gl:gl+2] = elem.nloc[i].xValue
            brhs[gl+1] = elem.nloc[i].yValue
            test = [elem.nloc[i].xValue , elem.nloc[i].yValue]
            print(K[gl,:])
        elif elem.nloc[i].DIRx:
            K[gl,gl] = 1
            brhs[gl] = elem.nloc[i].xValue
        elif elem.nloc[i].DIRy:
            K[gl+1,gl+1] = 1
            brhs[gl+1] = elem.nloc[i].yValue
        else :
            for j in range(4):
                nj = elem.nloc[j].nglob
                glj= nj*2
                K[gl,glj:glj+2] += Ke[i*2,j*2:j*2+1]
    return K , brhs

fileGmsh = sys.argv[1]
print(fileGmsh)
n, coord = Constructor('$Nodes',fileGmsh,None)
nelem, conect = Constructor('$Elements',fileGmsh,coord)

#Inicio de matriz global
K = sp.lil_matrix((n*2,n*2))
brhs = sp.lil_matrix((n*2,1))
#Formato de cond de borde
bcList = [[0,0],[1,None],[0,0],[None,1]]

for nodin in coord:
    if nodin.BG:
        nodin.physGrouptoValue(bcList)

for elemtest in conect:
    Ke = elemtest.getKe()
    K , bhrs = Assemble(elemtest, Ke , K , brhs)

print(K)

#K = K.tocsc
#brhs = brhs.tocsc

#U = spsolve(K,brhs)
#print(U)