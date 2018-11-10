def writeXML(coords,conect,nodedata,elemdata , nameFile,nodeStress):
    '''
    coords : Array with coordinates z,y
    coenct : Array with elements nodes
    nodedata: Array with node displacements
    elemdata: Array with Stressxx Stressyy, Shearxy
    '''
    nnode = coords.shape[0]
    nelem = elemdata.shape[0]
    fOut = open(nameFile.split(".")[0]+".vtu","w")
    #HEADER
    fOut.write( '<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
    fOut.write(' <UnstructuredGrid>\n')
    fOut.write('  <Piece NumberOfPoints="{}" NumberOfCells="{}">\n'.format(nnode,nelem))
    #ESCRIBO LOS DATOS DE LOS NODOS
    fOut.write('   <PointData Vectors="Tensiones">\n')
    fOut.write('    <DataArray NumberOfComponents="2" type="Float32" Name="Tensiones" format="ascii">\n    ')
    for data in nodeStress:
        fOut.write(' ' + str(data[0]) + ' ' + str(data[1]))
    fOut.write('\n    </DataArray>')
    fOut.write('\n   </PointData>\n')
    #ESCRIBO LOS DATOS DE LOS NODOS STRESS
    # fOut.write('   <PointData Vectors="Displacement">\n')
    # fOut.write('    <DataArray NumberOfComponents="2" type="Float32" Name="Displacement" format="ascii">\n    ')
    # for data in nodedata:
    #     fOut.write(' ' + str(data))
    # fOut.write('\n    </DataArray>')
    # fOut.write('\n   </PointData>\n')
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
    writeXML(1,1,1,1)