def writeXML(coords, conect, nodedata, nameFile, nodeStress):
    '''
    coords : Array with coordinates x & y (i have to add 0.0 below)
    conect : Array with elements nodes
    nodedata: Array with node displacements
    nodeStress: Array with stress in Nodes
    '''
    nnode = coords.shape[0]
    nelem = conect.shape[0]
    fOut = open(nameFile.split(".")[0]+".vtu","w")
    # fOut = open("Testingtesting.vtu","w")
    #HEADER
    fOut.write( '<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
    fOut.write(' <UnstructuredGrid>\n')
    fOut.write('  <Piece NumberOfPoints="{}" NumberOfCells="{}">\n'.format(nnode,nelem))
    #ESCRIBO LOS DATOS DE LOS NODOS STRESS
    fOut.write('   <PointData Vectors="Tensiones">\n')
    fOut.write('    <DataArray Name="Tensiones" NumberOfComponents="3" type="Float32"  format="ascii">\n    ')
    for data in nodeStress:
        fOut.write(' {} {} {}\n'.format(data[0],data[1],data[2]))
    fOut.write('\n    </DataArray>')
    #ESCRIBO LOS DATOS DE LOS NODOS DISPLACEMENTS
    fOut.write('    <DataArray Name="Displacement" NumberOfComponents="2" type="Float32" format="ascii">\n    ')
    for data in nodedata:
        fOut.write(' {} {}\n'.format(data[0],data[1]))
    fOut.write('\n    </DataArray>')
    fOut.write('\n   </PointData>\n')
    #DATA ELEMENTOS
    # fOut.write('   <CellData Tensors="Stress">\n')
    # fOut.write('    <DataArray NumberOfComponents="3" type="Float32" Name="Stress" format="ascii">\n    ')
    # for data in elemdata:
    #     fOut.write(' ' + str(data[0,0]) + ' ' + str(data[0,1]) + ' '+str(data[0,2]))
    # fOut.write('\n    </DataArray>')
    # fOut.write('\n   </CellData>\n')
    #DEFINO NODOS
    fOut.write('   <Points>\n')
    fOut.write('    <DataArray NumberOfComponents="3" type="Float32"  format="ascii">\n    ')
    for data in coords:
        fOut.write(' ' + str(data[0]) + ' ' + str(data[1]) + ' 0.0')
    fOut.write('\n    </DataArray>')
    fOut.write('\n   </Points>\n')
    #DEFINO ELEMENTOS
    fOut.write('   <Cells>\n')
    fOut.write('    <DataArray type="Int32" Name="connectivity" format="ascii" RangeMin="{}" RangeMax="{}">\n    '.format(1,nnode))
    for data in conect:
        fOut.write(' ' + str(int(data[0])) + ' ' + str(int(data[1]))  + ' '+ str(int(data[2]))  + ' '+ str(int(data[3])))
    fOut.write('\n    </DataArray>\n')

    fOut.write('    <DataArray type="Int32" Name="offsets" format="ascii" RangeMin="{}" RangeMax="{}">\n    '.format(4,4*nelem))
    aux = 4
    for data in conect:
        fOut.write(' ' + str(aux))
        aux +=4
    fOut.write('\n    </DataArray>\n')

    fOut.write('    <DataArray type="Int32" Name="types" format="ascii" RangeMin="{}" RangeMax="{}">\n    '.format(9,9))
    for data in conect:
        fOut.write(' ' + str(9))
    fOut.write('\n    </DataArray>')
    fOut.write('\n   </Cells>\n')

    fOut.write('  </Piece>\n')
    fOut.write(' </UnstructuredGrid>\n')
    fOut.write('</VTKFile>\n')

    fOut.close()
if __name__ == '__main__':
    writeXML(1,1,1,1,0)


def writeVTK(coords, conect, nodedata, nameFile, nodeStress):
    '''
    coords : Array with coordinates x & y (i have to add 0.0 below)
    conect : Array with elements nodes
    nodedata: Array with node displacements
    nodeStress: Array with stress in Nodes
    '''
    nnode = coords.shape[0]
    nelem = conect.shape[0]
    fOut = open(nameFile.split(".")[0]+".vtk","w")
    #HEADER
    fOut.write("# vtk DataFile Version 2.0\n")
    fOut.write("ASCII\n")
    fOut.write("DATASET UNSTRUCTURED_GRID\n \n")
    fOut.write("POINT {} float\n".format(nnode))
    for coordinate in coords:
        fOut.write("{} {} {}  \n".format(coordinate[0],coordinate[1],coordinate[2]))
    fOut.write("\n")

    fOut.write("CELLS {} {} \n".format(nelem,nelem*5))
    for nloc in conect:
        fOut.write("4 {} {} {} {} \n".format(int(nloc[0]),int(nloc[1]),int(nloc[2]),int(nloc[3])))
    fOut.write("\n")
    
    fOut.write("CELL_TYPES {}\n".format(nelem))
    for celltype in range(nelem):
        fOut.write("{}\n".format(4))
    fOut.write("\n")

    fOut.write("POINT_DATA {}\n".format(nnode))
    fOut.write("SCALARS Displacement float\n")
    fOut.write("LOOKUP_TABLE default\n")
    for scalar in nodedata:
        fOut.write("{}\n".format(scalar[0]))
    fOut.write("\n")
    fOut.close()