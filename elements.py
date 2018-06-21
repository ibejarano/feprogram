import numpy as np

class elementQ4:
    def __init__(self, coords):
        self.coords = np.matrix(coords)

def Hrs(r,s):
    h1 = 0.25*(1+r)*(1+s)
    h2 = 0.25*(1-r)*(1+s)
    h3 = 0.25*(1-r)*(1-s)
    h4 = 0.25*(1+r)*(1-s)
    funHrs = np.matrix([h1,h2,h3,h4])
    return funHrs

def dHrs(r,s):
    dhr1 = 0.25*(1+s)
    dhr2 = -0.25*(1+s)
    dhr3 = -0.25*(1-s)
    dhr4 = 0.25*(1-s)
    dhs1 = 0.25*(1+r)
    dhs2 = 0.25*(1-r)
    dhs3 = -0.25*(1-r)
    dhs4 = -0.25*(1+r)
    funDHrs = np.matrix(np.zeros((2,4)))
    funDHrs[0] = [dhr1,dhr2,dhr3,dhr4]
    funDHrs[1] = [dhs1,dhs2,dhs3,dhs4]
    return funDHrs

def funJ(r,s,Coords):
    #FALTA DEFINIR LAS COORD, SEGURAMENTE SERAN EN PUNTOS DE GAUSS
    return dHrs(r,s)*Coords

def funB(r,s,J):
    B = np.matrix(np.zeros((3,8)))
    Jinv = np.linalg.inv(J)
    DH = (dHrs(r,s).T * Jinv).T
    for i in range(2):
        B[i,i::2] = DH[i]
        B[1-i::2,::2] = DH[i]  
    return B

def Q4(NodesVec):
    Coords = np.matrix(NodesVec)
    gpoints = [-0.5773, 0.5773]
    J = np.matrix(np.zeros((2,2)))
    Ke = np.matrix(np.zeros((8,8)))
    for gpr in gpoints:
        for gps in gpoints:
            J = funJ(gpr,gps,Coords)
            B = funB(gpr,gps,J)
            detJ =  np.linalg.det(J)
            Kel += B.T * B
            print(Kel)
    return J

#LISTA PARA TESTEOS DE MATRICES, EL ELEMENTO 1 ESTA DISTORSIONADO
matcon = [[0,1,4,3],[1,2,5,4],[3,4,7,6],[4,5,8,7]]
coordenadas = [[3,2],[1,2],[0,2],[2,1],[1,1],[0,1],[2,0],[1,0],[0,0]]

#CONSTRUCTOR DE LISTA CON OBJETOS ELEMENTOS EN CADA UNO
lista = []
for elem in matcon:
    nodeloc = []
    for node in elem:
        nodeloc.append(coordenadas[node])
    lista.append(elementQ4(nodeloc))

coordstest = [[3,2],[-3,2],[-3,-2],[3,-2]]
test = Q4(coordstest)
