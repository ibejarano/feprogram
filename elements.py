#!/usr/bin/python3.6
import scipy.sparse as sp
import numpy as np
import math
from scipy.sparse.linalg import spsolve


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
        self.M = sp.lil_matrix((nrows,ncols))
        self.Kp = sp.lil_matrix((nrows,ncols))
        self.brhs = sp.lil_matrix((nrows,1))

        if self.elem.elemType == 'Quad9':
            gpsList = [-0.774,0, 0.774]
            gpsWei = [0.55555 , 0.88888 , 0.55555]
            self.H = self.quadH(gpsList,self.funH9,gpsWei)
            self.Hrs = self.quadHrs(gpsList,self.funHrs9,gpsWei)
        
        elif self.elem.elemType == 'Quad4':
            gpsList = [-0.5773,0.5773]
            gpsWei = [1,1]
            self.H = self.quadH(gpsList,self.funH4,gpsWei)
            self.Hrs = self.quadHrs(gpsList,self.funHrs4,gpsWei)

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

    def funH4(self,r,s):
        r_plus = 1+r
        r_minus = 1-r
        s_plus = 1+s
        s_minus = 1-s

        h4 = 0.25*r_plus*s_minus
        h3 = 0.25*r_minus*s_minus
        h2 = 0.25*r_minus*s_plus
        h1 = 0.25*r_plus*s_plus

        H=np.matrix([h1,h2,h3,h4])
        return H

    def funH9(self,r,s):
        r_power = 1-r**2
        s_power = 1-s**2
        r_plus = 1+r
        r_minus = 1-r
        s_plus = 1+s
        s_minus = 1-s

        h9 = r_power*s_power
        h8 = 0.5*s_power*r_plus - 0.5*h9
        h7 = 0.5*r_power*s_minus - 0.5*h9
        h6 = 0.5*s_power*r_minus - 0.5*h9
        h5 = 0.5*r_power*s_plus - 0.5*h9
        h4 = 0.25*r_plus*s_minus - 0.5*h7 - 0.5*h8 -0.25*h9
        h3 = 0.25*r_minus*s_minus - 0.5*h6 - 0.5*h7 - 0.25*h9
        h2 = 0.25*r_minus*s_plus - 0.5*h5 - 0.5*h6 - 0.25*h9
        h1 = 0.25*r_plus*s_plus - 0.5*h5 -0.5*h8 -0.25*h9

        H=np.matrix([h1,h2,h3,h4,h5,h6,h7,h8,h9])

        return H

    def funHrs4(self,r,s):
        hrquad4 = [1+s , -1-s , -1+s ,1-s]
        hsquad4 = [1+r,1-r , -1 +r , -1 -r]

        hr = [0]*4
        hs = [0]*4
        for i in range(4):
            hr[i] = 0.25*hrquad4[i]
            hs[i] = 0.25*hsquad4[i]

        dhrs = np.matrix(np.zeros((2,4)))

        dhrs[0] = hr
        dhrs[1] = hs
        return dhrs    

    def funHrs9(self,r,s):
        hrquad4 = [1+s , -1-s , -1 +s ,1-s]
        hrquad8 = [-2*r*(1+s),-1*(1-s**2),-2*r*(1-s),(1-s**2)]
        hrquad9 = -2*r*(1-s**2)
        hsquad4 = [1+r,1-r ,-1+r ,-1-r]
        hsquad8 = [(1-r**2),-2*s*(1-r),-1*(1-r**2),-2*s*(1+r)]
        hsquad9 = -2*s*(1-r**2)

        hr = [0]*9
        hs = [0]*9

        for i in range(4):
            hr[i+4] = 0.5*hrquad8[i] - 0.5*hrquad9
            hs[i+4] = 0.5*hsquad8[i] - 0.5*hsquad9

        hr[8] = hrquad9
        hs[8] = hsquad9
        
        for i in range(1,4):
            hr[i] = 0.25*hrquad4[i] - 0.5*hr[i+4] - 0.5*hr[i+3] - 0.25*hrquad9
            hs[i] = 0.25*hsquad4[i] - 0.5*hs[i+4] - 0.5*hs[i+3] - 0.25*hsquad9
            
        hr[0] = 0.25*hrquad4[0] - 0.5*hr[4] - 0.5*hr[7] - 0.25*hrquad9
        hs[0] = 0.25*hsquad4[0] - 0.5*hs[4] - 0.5*hs[7] - 0.25*hsquad9
        
        dhrs = np.matrix(np.zeros((2,9)))

        dhrs[0] = hr
        dhrs[1] = hs
        return dhrs    

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
        mu = 1e-2
        rho = 1e3
        dof = self.elem.nnodesloc*2
        J = np.matrix(np.zeros((2,2)))
        Ke = np.matrix(np.zeros((dof,dof)))
        Me = np.matrix(np.zeros((dof,dof)))
        Kpe = np.matrix(np.zeros((dof,dof)))
        m = np.matrix(np.array([1,1,0,1]).reshape((4,1)))
        I = np.matrix(np.eye(4))
        elemCorners = self.elem.getpos(elemNodeTags, self.coordinates)
        Hv = np.matrix(np.zeros((2,8)))
        for i in range(int(dof/2)):
            J = self.Hrs[i] * elemCorners
            detJ = np.linalg.det(J)
            B = self.elem.funB(self.Hrs[i],J)
            Hv[0,::2] = self.H[i]
            Hv[1,1::2] = self.H[i]
            Ke += 2*mu* B.T * (I - 0.333 * m * m.T).T * (I - 0.333 * m * m.T) * B * detJ
            Me += rho * Hv.T * Hv * detJ
            Kpe += B.T * m * m.T * B * detJ
        return Ke, Me, Kpe

    def assemble(self):
        #Armado de matrices elementales y ensamblaje
        row = []
        col = []
        Kelem = []
        Melem = []
        Kpelem = []
        deltaP = 1e-9

        for elem in self.conectivity:
            Ke, Me , Kpe = self.getKe(elem-1)
            #self.K , self.brhs = Assemble(elem, Ke , self.K , self.brhs)
            rowap , colap = self.getRowData(elem)
            Kelemap = Ke.ravel()
            Melemap = Me.ravel()
            Kpelemap = Kpe.ravel()
            row.extend(rowap)
            col.extend(colap)
            Kelem.append(Kelemap)
            Melem.append(Melemap)
            Kpelem.append(Kpelemap)

        #Setup de matrices esparsas
        Kelem = np.array(Kelem).ravel()
        Melem = np.array(Melem).ravel()
        Kpelem = np.array(Kpelem).ravel()

        self.K = sp.coo_matrix((Kelem,(row,col)),shape=(self.nnode*2,self.nnode*2)).tolil()
        self.M = sp.coo_matrix((Melem,(row,col)),shape=(self.nnode*2,self.nnode*2)).tolil()
        self.Kp = sp.coo_matrix((Kpelem,(row,col)),shape=(self.nnode*2,self.nnode*2)).tolil()

        for dof2Set in self.dofDir:
            dof2Set = int(dof2Set)
            _, J = self.K[dof2Set,:].nonzero()
            for j in J:
                if j != dof2Set:
                    self.K[dof2Set,j] = 0.0
            
            self.K[dof2Set,dof2Set] = 1.0
            self.brhs[dof2Set] = 0

        for dofNeu in self.dofNeu:
            self.brhs[dofNeu] = deltaP

        self.K = self.K.tocsc()
        self.M = self.M.tocsc()
        penalizacion = 1e2
        self.Kp = (self.Kp * penalizacion).tocsc()
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

            dofSetDir = [[0,1],[1],[1],[1]]
            dirDof = set()
            neuDof = set()

            for enu, nodes in enumerate(bcNodesList):
                dof = [(node-1)*2 +j for node in nodes for j in range(2)]
                for k in dofSetDir[enu]:
                    dirDof.update(dof[k::2])
                    if enu == 3:
                        neuDof.update(set(dof) - set(dof[k::2]))

            self.dofDir = list(dirDof)
            self.dofNeu = list(neuDof)

    def getLocalVec(self, ind, Vec):
        """Dados unos indices, devuelve un vector con esos grados de libertad"""
        return Vec[ind].reshape((2,4))

    def getNe(self, elemNodeTags, velocity):
        rho = 1e3
        dof = self.elem.nnodesloc*2
        J = np.matrix(np.zeros((2,2)))
        Ne = np.matrix(np.zeros((dof,dof)))
        elemCorners = self.elem.getpos(elemNodeTags, self.coordinates)
        Hv = np.matrix(np.zeros((2,8)))
        dofLocal= [ int(i*2+j) for i in elemNodeTags for j in range(2)]
        velLocal = self.getLocalVec(dofLocal, velocity)
        for i in range(int(dof/2)):
            J = self.Hrs[i] * elemCorners
            detJ = np.linalg.det(J)
            Hv[0,::2] = self.H[i]
            Hv[1,1::2] = self.H[i]
            H_grad = self.elem.getHgrad(self.Hrs[i],J)
            Ne += Hv.T * rho * velLocal * H_grad * detJ
        return Ne

    def assemble_N(self, velGlobal):
        #Armado de matrices elementales y ensamblaje
        row = []
        col = []
        Nelem = []

        for elem in self.conectivity:
            Ne = self.getNe(elem-1, velGlobal)
            #self.K , self.brhs = Assemble(elem, Ke , self.K , self.brhs)
            rowap , colap = self.getRowData(elem)
            Nelemap = Ne.ravel()
            row.extend(rowap)
            col.extend(colap)
            Nelem.append(Nelemap)

        #Setup de matrices esparsas
        Nelem = np.array(Nelem).ravel()

        return sp.coo_matrix((Nelem,(row,col)),shape=(self.nnode*2,self.nnode*2)).tocsc()

    def initialVector(self):
        return np.ones(self.nnode*2)

    def solveProblem(self):
        V_prev = self.initialVector()
        N = self.assemble_N(V_prev)
        tol = 1e-11
        for ite in range(50):
            V = spsolve(self.K + N , self.brhs)
            err = np.linalg.norm(V - V_prev , np.inf)
            if err < tol:
                print("INFO: Problem Converged!")
                break
            else:
                print("INFO: error : {}".format(err))
                print("INFO: Iterating number: {}".format(ite))
            V_prev = V
            N = self.assemble_N(V_prev)

        return V
