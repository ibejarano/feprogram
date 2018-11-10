#!/usr/bin/python3.6

import numpy as np
import scipy.sparse as sp
import sys
import logging
from gmshtools import readGmshFile
from scipy.sparse.linalg import spsolve
from elements import Node2D , ElemQ4 , Assemble , Constructor ,Ulocal
from datetime import datetime
from python_to_xml import writeXML


t1 = datetime.now()
logging.basicConfig(level='INFO')

fileGmsh = sys.argv[1]
nNodos, coord = Constructor('$Nodes',fileGmsh,None)
nelem, conect = Constructor('$Elements',fileGmsh,coord)

He = conect[0].buildMatrix()
C = conect[0].Cmat()

#Inicio de matriz global
K = sp.lil_matrix((nNodos*2,nNodos*2))
brhs = sp.lil_matrix((nNodos*2,1))
#Formato de cond de borde
bcList = [[0,0],[100,0]]

for nodin in coord:
    if nodin.BG:
        nodin.physGrouptoValue(bcList)

# # TestH = conect[0].fundHrs(-1,1)

for elemtest in conect:
    Ke = elemtest.getKe(He,C)
    K , brhs = Assemble(elemtest, Ke , K , brhs)
    # elemtest.StressPost(CBe)

K = K.tocsc()
brhs = brhs.tocsc()

U = spsolve(K,brhs)

Stress = np.matrix(np.zeros((nelem,3)))

for count,elem in enumerate(conect):
    B = elem.Bstress
    Uloc = Ulocal(U,elem.nloc)
    Stress[count] = elem.Stress(Uloc,C,He)
    for nodalStress in elem.localNodes(): #Lista con nodos globales del elemento
            coord[nodalStress-1].storeStress(Stress[count])


nelem = len(conect)
outCoords = np.zeros((nNodos,3))
outNodesStress = np.zeros((nNodos,3))
outConect = np.zeros((nelem,4))

for i in range(nNodos):
    outCoords[i] = [coord[i].x , coord[i].y , 0.0]
    outNodesStress[i] = coord[i].stress

for i in range(nelem):
    outConect[i] = [conect[i].nloc[k].nglob -1  for k in range(4)]

writeXML(outCoords, outConect ,U,Stress , sys.argv[1],outNodesStress)


# t2 = datetime.now()

# logging.info('Tiempo total: %f sec', (t2 - t1).total_seconds())

