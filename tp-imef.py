#!/usr/bin/python3.6

import numpy as np
import scipy.sparse as sp
import sys
import logging
from gmshtools import readGmshFile
from scipy.sparse.linalg import spsolve
from elements import Node2D , ElemQ4 , Assemble , Constructor
from datetime import datetime
from python_to_xml import writeXML

def storeValuesToNodes(nodeArray,values):
        '''
        Solo almacena los desplazamientos calculados en los nodos
        '''
        nNodos = len(nodeArray)
        values = U.reshape(nNodos,2)
        for i in range(nNodos):
                nodeArray[i].storeCalcValue(values[i,0],values[i,1])
        return None



t1 = datetime.now()
logging.basicConfig(level='INFO')

fileGmsh = sys.argv[1]
nNodos, coord = Constructor('$Nodes',fileGmsh,None)
nelem, conect = Constructor('$Elements',fileGmsh,coord)

C = conect[0].Cmat()

#Inicio de matriz global
K = sp.lil_matrix((nNodos*2,nNodos*2))
brhs = sp.lil_matrix((nNodos*2,1))

#Formato de cond de borde
bcList = [[0,0],[100,0]]

for nodin in coord:
    if nodin.BG:
        nodin.physGrouptoValue(bcList)

for elem in conect:
    Ke = elem.getKe(C)
    K , brhs = Assemble(elem, Ke , K , brhs)
    # elemtest.StressPost(CBe)

K = K.tocsc()
brhs = brhs.tocsc()

U = spsolve(K,brhs)

storeValuesToNodes(coord,U)

Stress = np.matrix(np.zeros((nelem,3)))
nelem = len(conect)
outCoords = np.zeros((nNodos,3))
outNodesStress = np.zeros((nNodos,3))
outConect = np.zeros((nelem,4))

for i in range(nelem):
    outConect[i] = [conect[i].nloc[k].nglob -1  for k in range(4)]
    elemStress = conect[i].computeStress(C).T
    Stress[i] = elemStress
    conect[i].elemStressToNodes(elemStress)
    
for i in range(nNodos):
    outCoords[i] = [coord[i].x , coord[i].y , 0.0]
    outNodesStress[i] = coord[i].stress

writeXML(outCoords, outConect ,U,Stress, sys.argv[1],outNodesStress)


# t2 = datetime.now()

# logging.info('Tiempo total: %f sec', (t2 - t1).total_seconds())

