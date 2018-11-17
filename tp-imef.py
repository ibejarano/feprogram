#!/usr/bin/python3.6

import numpy as np
import scipy.sparse as sp
import sys
import logging
from gmshtools import readGmshFile
from scipy.sparse.linalg import spsolve
from elements import Node2D , ElemQ4 , Assemble , Constructor , FemProblem
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

def computeDomainStress(elemMat,Hrs):
        nelem = len(elemMat)
        outConect = np.zeros((nelem,4))
        Stress = np.matrix(np.zeros((nelem,3)))

        for elemRow , elem in enumerate(elemMat):
                outConect[elemRow] = [elem.nloc[k].nglob -1  for k in range(4)]
                elemStress = elem.computeStress(C,Hrs).T
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

def setBC(coordenadas , bcList):
        #Seteo de BC en nodos
        for nodin in coordenadas:
                if nodin.BG:
                        nodin.physGrouptoValue(bcList)
        pass

t1 = datetime.now()
logging.basicConfig(level='INFO')

fileGmsh = sys.argv[1]
nNodos, coord = Constructor('$Nodes',fileGmsh,None)
nElem, conect = Constructor('$Elements',fileGmsh,coord)

t2 = datetime.now()
logging.info('Leer de Gmsh: %f sec', (t2 - t1).total_seconds())

C = conect[0].Cmat()


if len(conect[0].nloc) == 9:
        elemType = 'Quad9'
else:
        elemType = 'Quad4'


        #Inicio de matriz global
        K = sp.lil_matrix((nNodos*2,nNodos*2))
        brhs = sp.lil_matrix((nNodos*2,1))

        #Formato de cond de borde
        bcList = [[0,0],[400,200]]

        #Seteo de BC en nodos
        for nodin in coord:
                if nodin.BG:
                        nodin.physGrouptoValue(bcList)

        #Armado de matrices elementales y ensamblaje
        for elem in conect:
                Ke = elem.getKe(C)
                K , brhs = Assemble(elem, Ke , K , brhs)

        #Setup de matrices esparsas
        K = K.tocsc()
        brhs = brhs.tocsc()

        #Resolucion
        logging.info('Resolviendo...')
        U = spsolve(K,brhs)
        logging.info('Resuelto')
        storeValuesToNodes(coord,U)

        #Preparacion de Arrays para postproceso
        conectivity = computeDomainStress(conect)
        nodeCoordinates , nodeStress = createNodeData(coord)
        U = U.reshape((nNodos,2))

        #Escribir archivo .vtu para ver en Paraview
        writeXML(nodeCoordinates, conectivity , U, sys.argv[1], nodeStress)

fem = FemProblem(fileGmsh,nElem,nNodos,elemType,conect)

t3 = datetime.now()
logging.info('Crear elemento fem: %f sec', (t3 - t2).total_seconds())

fem.setMatrix()

t4 = datetime.now()
logging.info('Seteo inicial de matrices: %f sec', (t4 - t3).total_seconds())

#Seteo de condiciones de borde sobre los nodos
bcList = [[0,0],[50,0]]

setBC(coord, bcList)

t5 = datetime.now()
logging.info('Seteo de BC en nodos: %f sec', (t5 - t4).total_seconds())

fem.assemble(C)

t6 = datetime.now()
logging.info('Ensamblaje de matrices: %f sec', (t6 - t5).total_seconds())

U = spsolve(fem.K,fem.brhs)

t7 = datetime.now()
logging.info('Calculo de desplazamientos: %f sec', (t7 - t6).total_seconds())

storeValuesToNodes(coord,U)
conectivity = computeDomainStress(conect,fem.Hrs[8])
nodeCoordinates , nodeStress = createNodeData(coord)
U = U.reshape((nNodos,2))

#Escribir archivo .vtu para ver en Paraview
writeXML(nodeCoordinates, conectivity , U, sys.argv[1], nodeStress)


t3 = datetime.now()
logging.info('Tiempo total: %f sec', (t3 - t1).total_seconds())

