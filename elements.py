import numpy as np

class elementQ4:
    def __init__(self,num, nodes):
        self.num = num
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

    def funJ(self,r,s):
        #FIXME: FUNCION NO IMPLEMENTADA, BORRAR SI NO SE VA A USAR
        #J = np.matrix(np.zeros((2,2)))
        #J = dHrs*self.Pos
        return None

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

class node:
    def __init__(self,num,coords):
        self.ng = num
        self.x = coords[0]
        self.y = coords[1]
        self.DIR = False
        self.NEU = False
    def setDIR(self):
        self.DIR = True
    def setNEU(self):
        self.NEU = True

def funHrs(r,s):
    #FIXME: FUNCION NO IMPLEMENTADA, CONTIENE LAS FUNCIONES DE INTERPOLACION
    h1 = 0.25*(1+r)*(1+s)
    h2 = 0.25*(1-r)*(1+s)
    h3 = 0.25*(1-r)*(1-s)
    h4 = 0.25*(1+r)*(1-s)
    hrs = np.matrix([h1,h2,h3,h4])
    return hrs

#LISTA PARA TESTEOS DE MATRICES, EL ELEMENTO 1 ESTA DISTORSIONADO
matcon = [[0,1,4,3],[1,2,5,4],[3,4,7,6],[4,5,8,7]]
coordenadas = [[3,2],[1,2],[0,2],[2,1],[1,1],[0,1],[2,0],[1,0],[0,0]]

nodelist = []
for numglobal , coord in enumerate(coordenadas):
    nodelist.append(node(numglobal, coord))

#CONSTRUCTOR DE LISTA CON OBJETOS ELEMENTOS EN CADA UNO
elemlist = []
for num,elem in enumerate(matcon):
    locallist = []
    for node in elem:
        locallist.append(nodelist[node])
    elemlist.append(elementQ4(num,locallist))

#coordstest = [[3,2],[-3,2],[-3,-2],[3,-2]]
#test = Q4(coordstest)
print(elemlist[1].getKe())