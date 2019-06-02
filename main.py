#!/usr/bin/env python3
import numpy as np
import scipy.sparse as sp
import sys
import logging
from scipy.sparse.linalg import spsolve
from elements import FemProblem
from datetime import datetime
from python_to_xml import writeXML , writeXML_nodeData
from gmsh_api import getMeshInfo

t1 = datetime.now()
logging.basicConfig(level='INFO')
conectivity ,coordinates, bcNodeTags = getMeshInfo()

t2 = datetime.now()
logging.info('Crear malla: %f sec', (t2 - t1).total_seconds())

fem = FemProblem(conectivity, coordinates)

t3 = datetime.now()
logging.info('Crear elemento fem: %f sec', (t3 - t2).total_seconds())

fem.setMatrix()

t4 = datetime.now()
logging.info('Seteo inicial de matrices: %f sec', (t4 - t3).total_seconds())

#Seteo de condiciones de borde sobre los nodos
fem.setBoundaryConditions(bcNodeTags)

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


#conectivity_xml = np.zeros((len(conectivity),elemNodes))

#for ind , local in enumerate(conectivity):
        #conectivity_xml[ind] = local.localNodes()

#Escribir archivo .vtu para ver en Paraview
#writeXML(nodeCoordinates, conectivity_xml , U, sys.argv[1], nodeStress)
elemSkip = len(bcTupledNodes)+1
writeGmshOut(fileGmsh,U,vonMisesStress,elemSkip)

t8 = datetime.now()
logging.info('Tiempo total: %f sec', (t8 - t1).total_seconds())