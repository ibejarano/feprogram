#!/usr/bin/python3.6
import scipy.sparse as sp
import numpy as np
import math
from scipy.sparse.linalg import spsolve
from tools.interpFunctions import *


class Elem:
    def __init__(self,nodes,elemType):
        self.elemType = elemType
        self.nnodesloc = nodes
        #self.nloc = nodes (This are ojects)

    def funB(self,dHrs,J):
        dof = self.nnodesloc*2
        B = np.matrix(np.zeros((4,dof)))
        Jinv = np.linalg.inv(J)
        Dh = Jinv * dHrs
        for i in range(2):
            B[i,i::2] = Dh[i]
            B[2,i::2] = Dh[i-1]
        return B

    def getHgrad(self,dHrs, J):
        dof = self.nnodesloc*2
        H_grad = np.matrix(np.zeros((4,dof)))
        Jinv = np.linalg.inv(J)
        Dh = Jinv * dHrs
        for i in range(2):
            H_grad[i*2:i*2+2,i::2] = Dh
        return H_grad

    def getpos(self, nodesTagList, coordinates):
        auxlist = []
        for node in nodesTagList:
            auxlist.append([coordinates[node,0], coordinates[node,1]])
        return np.matrix(auxlist)   

class FemProblem:
    def __init__(self,conectivity, coordinates):
        '''Instance the FEM problem, the code get around this class
        
        Arguments:
            meshName {string} -- name of mesh, DEPRECATED
            nelem {integer} -- total element number of problem
            nnode {integer} -- total node number of problem
            elemType {string} -- declares if type is quad9 or quad4
            conectivity {array} -- contains objects class Elem
            bcNodes {array} -- contains a list with nodes with boundary conditions
        '''
        self.nnode = coordinates.shape[0]
        self.nelem = len(conectivity)

        self.elem = Elem(4,'Quad4') if len(conectivity[0])==4 else Elem(9,'Quad9')

        self.conectivity = conectivity
        self.coordinates = coordinates

    def setMatrix(self):
        """[summary]
            This function setups Empty K , brhs in all cases
            H & Hrs depends on ElemType

        Raises:
            Exception -- if no elemType defined nor implemented yet

        Returns:
            K[lil_matrix] -- Empty Global rigidity matrix
            brhs[lil_matrix] -- Empty Right hand side vector
            H[array] -- Array with form functions evaluated in gps
            Hrs[array] -- It haves the derivatives in gps
        """
        ncols = (self.nnode)*2
        nrows = (self.nnode)*2
        self.K = sp.lil_matrix((nrows,ncols))
        self.brhs = sp.lil_matrix((nrows,1))

        if self.elem.elemType == 'Quad9':
            gpsList = [-0.774,0, 0.774]
            gpsWei = [0.55555 , 0.88888 , 0.55555]
            self.H = self.quadH(gpsList,funH9,gpsWei)
            self.Hrs = self.quadHrs(gpsList,funHrs9,gpsWei)
        
        elif self.elem.elemType == 'Quad4':
            gpsList = [-0.5773,0.5773]
            gpsWei = [1,1]
            self.H = self.quadH(gpsList,funH4,gpsWei)
            self.Hrs = self.quadHrs(gpsList,funHrs4,gpsWei)

        else:
            raise Exception('Invalid elemType defined you must use Quad4 or Quad9')

    def quadH(self,gpoints,formFunction,gpweights=None):
        '''Build the H matrix with the formfunctions and evaluates at gauss points and weights if its needed
        
        Arguments:
            gpoints {list} -- list containing iterable points
            formFunction {callback} -- function that returns the value of form functions at gauss point
        
        Keyword Arguments:
            gpweights {list} -- weights of numerical integration, only needed when elemType=Quad9 (default: {None})
        
        Returns:
            array -- The columns are hi form function and row are gauss points number
        '''

        H = []
        for gpr in gpoints:
            for gps in gpoints:
                dH = formFunction(gpr,gps)
                H.append(dH)
        return H

    def quadHrs(self,gpoints,formFunction,gpweights=None):
        '''Creates an array that has the derivatives of form
        function evaluated at gauss points
        
        Arguments:
            gpoints {list[float]} -- List with gauss points
            formFunction {callback} -- returns the gp evaluated
        
        Keyword Arguments:
            gpweights {list[float]} -- it weight numerical integration (default: {None})
        
        Returns:
            He [array] -- contains the functions evaluated at gausspoints and ordered
        '''

        He = []
        gpWei = []
        for weir , gpr in enumerate(gpoints):
            for weis, gps in enumerate(gpoints):
                gprWei = gpweights[weir]
                gpsWei = gpweights[weis]
                dHrs = formFunction(gpr,gps)
                gpWei.append(gprWei*gpsWei)
                He.append(dHrs)
        self.gpWei = gpWei
        return He


    def getRowData(self,elemNodesTag):
        '''
        Optimized one
        '''
        dof = len(elemNodesTag)*2
        listedNodes = [(x-1)*2+y for x in elemNodesTag for y in [0,1]]
        col = listedNodes*dof
        row = [ind for ind in listedNodes for i in range(dof)]
        return row , col

    def getKe(self, elemNodeTags):
        dof = self.elem.nnodesloc*2
        J = np.matrix(np.zeros((2,2)))
        Ke = np.matrix(np.zeros((dof,dof)))
        elemCorners = self.elem.getpos(elemNodeTags, self.coordinates)
        for i in range(int(dof/2)):
            J = self.Hrs[i] * elemCorners
            detJ = np.linalg.det(J)
            B = self.elem.funB(self.Hrs[i],J)
            Ke += B.T * B * detJ
        return Ke

    def assemble(self):
        #Armado de matrices elementales y ensamblaje
        row = []
        col = []
        Kelem = []

        for elem in self.conectivity:
            Ke = self.getKe(elem-1)
            #self.K , self.brhs = Assemble(elem, Ke , self.K , self.brhs)
            rowap , colap = self.getRowData(elem)
            Kelemap = Ke.ravel()
            row.extend(rowap)
            col.extend(colap)
            Kelem.append(Kelemap)

        #Setup de matrices esparsas
        Kelem = np.array(Kelem).ravel()

        self.K = sp.coo_matrix((Kelem,(row,col)),shape=(self.nnode*2,self.nnode*2)).tolil()

        for dof2Set in self.dofDir:
            dof2Set = int(dof2Set)
            _, J = self.K[dof2Set,:].nonzero()
            for j in J:
                if j != dof2Set:
                    self.K[dof2Set,j] = 0.0
            
            self.K[dof2Set,dof2Set] = 1.0
            self.brhs[dof2Set] = 0

        for dofNeu in self.dofNeu:
            
            self.brhs[dofNeu] = 1 #FIXME : Hardcoded force value

        self.K = self.K.tocsc()
        self.brhs = self.brhs.tocsc()

    def setBoundaryConditions(self, bcNodesList):
            '''
            it only receives a tuple with nodes and assign only forces and construct brhs from element
            '''
            #the first One is Dirichlet and its meant to be K = 1 in that index
            #if NEU brhs equals force
            #bcNodeTags[0,1] = Borde inferior x y
            #bcNodeTags[2,3] = Borde derecho  x y
            #bcNodeTags[4,5] = Borde superior x y
            #bcNodeTags[6,7] = Borde izquierdo x y

            dofSetDir = [[0,1],[],[],[]]
            dofSetNeu = [[],[],[],[0,1]]
            dirDof = set()
            neuDof = set()

            for enu, nodes in enumerate(bcNodesList):
                dof = [int((node-1)*2 +j) for node in nodes for j in range(2)]
                for k in dofSetDir[enu]:
                    dirDof.update(dof[k::2])

                for m in dofSetNeu[enu]:
                    neuDof.update(dof[m::2])

            print(dirDof)
            print(neuDof)
            self.dofDir = list(dirDof)
            self.dofNeu = list(neuDof)

    def getLocalVec(self, ind, Vec):
        """Dados unos indices, devuelve un vector con esos grados de libertad"""
        return Vec[ind].reshape((2,4))

    def initialVector(self):
        return np.ones(self.nnode*2)

    def assembleMatrix(self, vec):
        return 0

    def solveProblem(self):
        V = spsolve(self.K , self.brhs)
        return V