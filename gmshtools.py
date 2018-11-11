
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