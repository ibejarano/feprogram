import shutil

def searchline(palabra, f ):
    print("Buscando: ",palabra)
    linea = f.readline()
    condicion = (palabra in linea)
    while not condicion:
        linea = f.readline()
        condicion = (palabra in linea)
    if palabra in linea:
        # found linea
        return linea, f
    else:
        print (palabra, ' not found')
        return 0
    return 0

def readGmshFile(word,inpGmsh):
    '''
    INPUT
    word: String con palabra a buscar
    inpgmsh: Archivo en formato .gmsh donde busca word
    OUTPUT
    nitems: Cantidad de elementos hallados
    Parser del file de Gmsh: Lista de listas con datos sobre el item necesario
    '''
    fileGmsh = open(inpGmsh,"r")
    lines, fileGmsh = searchline(word,fileGmsh)
    lines = fileGmsh.readline()
    nitems = int(lines.split()[0])
    return nitems , fileGmsh

def writeGmshOut(inputGmsh,U):
    """A Function to write the output proccessing file
    
    Arguments:
        inputGmsh {string} -- Only a string with the name of gmsh
        U {array} --          Displacement array vector
    """
    dim = 2
    nNodes = len(U)
    outputGmsh = 'Out'+inputGmsh
    shutil.copyfile(inputGmsh,outputGmsh)
    fOut = open(outputGmsh,'a')
    fOut.write( '$NodeData\n1\n"Displacement [mm]"\n1\n0\n3\n0\n3\n')
    fOut.write('{} \n'.format(nNodes))
    counter = 1
    for displacement in U:
        fOut.write('{} {} {} 0\n'.format(counter,displacement[0],displacement[1]))
        counter+=1
    fOut.write('$EndNodeData')