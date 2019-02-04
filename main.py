#!/usr/bin/env python3
import numpy as np
import scipy.sparse as sp
import sys
import logging
from gmshtools import readGmshFile , writeGmshOut
from scipy.sparse.linalg import spsolve
from elements import Node2D , Constructor , FemProblem
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

def createNodeData(nodeMat):
        outCoords = np.zeros((nNodos,3))
        outNodesStress = np.zeros((nNodos,3))
        for nodeRow, node in enumerate(nodeMat):
                outCoords[nodeRow] = [node.x , node.y , 0.0]
                outNodesStress[nodeRow] = node.stress
        return outCoords , outNodesStress

def removeDuplicates(tupledNodes):
        listNodes = []
        for i in tupledNodes:
                for j in i:
                        listNodes.append(j)
        return listNodes

t1 = datetime.now()
logging.basicConfig(level='INFO')

fileGmsh = sys.argv[1]
nNodos, coord = Constructor('$Nodes',fileGmsh,None)
nNodesBc , bcTupledNodes = Constructor('$Elements',fileGmsh,coord,bc=True)
bcNodes = removeDuplicates(bcTupledNodes)
nElem, conect = Constructor('$Elements',fileGmsh,coord)

t2 = datetime.now()
logging.info('Leer de Gmsh: %f sec', (t2 - t1).total_seconds())

if len(conect[0].nloc) == 9:
        elemType = 'Quad9'
        elemNodes = 9
else:
        elemType = 'Quad4'
        elemNodes = 4

fem = FemProblem(nElem,nNodos,elemType,conect,bcNodes)

t3 = datetime.now()
logging.info('Crear elemento fem: %f sec', (t3 - t2).total_seconds())

fem.setMatrix()

t4 = datetime.now()
logging.info('Seteo inicial de matrices: %f sec', (t4 - t3).total_seconds())

#Seteo de condiciones de borde sobre los nodos
fem.setBoundaryConditions(coord)

t5 = datetime.now()
logging.info('Seteo de BC en nodos: %f sec', (t5 - t4).total_seconds())

fem.assemble()

t6 = datetime.now()
logging.info('Ensamblaje de matrices: %f sec', (t6 - t5).total_seconds())

U = spsolve(fem.K,fem.brhs)

t7 = datetime.now()
logging.info('Calculo de desplazamientos: %f sec', (t7 - t6).total_seconds())

storeValuesToNodes(coord,U)
vonMisesStress = fem.getStress()
nodeCoordinates , nodeStress = createNodeData(coord)
U = U.reshape((nNodos,2))


conectivity_xml = np.zeros((len(conect),elemNodes))

for ind , local in enumerate(conect):
        conectivity_xml[ind] = local.localNodes()

#Escribir archivo .vtu para ver en Paraview
#writeXML(nodeCoordinates, conectivity_xml , U, sys.argv[1], nodeStress)
elemSkip = len(bcTupledNodes)+1
writeGmshOut(fileGmsh,U,vonMisesStress,elemSkip)

t8 = datetime.now()
logging.info('Tiempo total: %f sec', (t8 - t1).total_seconds())

