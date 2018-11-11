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

def computeDomainStress(elemMat):
        nelem = len(elemMat)
        outConect = np.zeros((nelem,4))
        Stress = np.matrix(np.zeros((nelem,3)))

        for elemRow , elem in enumerate(elemMat):
                outConect[elemRow] = [elem.nloc[k].nglob -1  for k in range(4)]
                elemStress = elem.computeStress(C).T
                Stress[elemRow] = elemStress
                elemMat[elemRow].elemStressToNodes(elemStress)
        return outConect


def createNodeData(nodeMat):
        outCoords = np.zeros((nNodos,3))
        outNodesStress = np.zeros((nNodos,3))
        for nodeRow, node in enumerate(nodeMat):
                outCoords[nodeRow] = [node.x , node.y , 0.0]
                outNodesStress[nodeRow] = node.stress
        return outCoords , outNodesStress

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
bcList = [[0,0],[400,200]]

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

conectivity = computeDomainStress(conect)

nodeCoordinates , nodeStress = createNodeData(coord)

writeXML(nodeCoordinates, conectivity , U, sys.argv[1], nodeStress)


# t2 = datetime.now()

# logging.info('Tiempo total: %f sec', (t2 - t1).total_seconds())

