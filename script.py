# coding: utf-8
import sys
import numpy as np

from vtk import *
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

import integration_point_pb2 as ip


def readElementData(ip_data, offsets):

    if offsets.GetNumberOfComponents() != 1:
        print('Error: integration point offsets array has wrong number of components.')
        print('Expected one component.')
        sys.exit(2)
    if offsets.GetDataType() != VTK_UNSIGNED_LONG:
        print('Error: integration point offsets array has wrong data type.')
        print('Expected VTK_UNSIGNED_LONG.')
        sys.exit(2)

    if ip_data.GetNumberOfComponents() != 1:
        print('Error: integration point data array has wrong number of components.')
        print('Expected one component.')
        sys.exit(2)
    if ip_data.GetDataType() != VTK_CHAR:
        print('Error: integration point data array has wrong data type.')
        print('Expected VTK_CHAR.')
        sys.exit(2)

    ip_data_array = vtk_to_numpy(ip_data)

    element_data = []       # Resulting array

    # Process all but the last one.
    for i in range(offsets.GetNumberOfTuples() - 1):
        start = int(offsets.GetTuple1(i))
        end = int(offsets.GetTuple1(i+1))
        ed = ip.ElementData()
        ed.ParseFromString(ip_data_array[start:end].tobytes())
        element_data.append(ed)
    # Process the last element
    start = int(offsets.GetTuple1(i+1))
    end = ip_data.GetNumberOfTuples()
    if start >= end:
        print('Error: integration point offsets and integration point data \
               array sizes do not match.')
        sys.exit(2);
    ed = ip.ElementData()
    ed.ParseFromString(ip_data_array[start:end].tobytes())
    element_data.append(ed)

    return element_data

def averageOverIntegrationPoints(element_data):
    if element_data.n_integration_points == 0:
        print("Error: Got 0 integration points.")
        sys.exit(2)

    sigma = np.array([0.0, 0.0, 0.0, 0.0])
    for i in range(element_data.n_integration_points):
        s = element_data.sigma[i]
        for j in range(4):
            sigma[j] += s.value[j]
    return sigma / element_data.n_integration_points

def readMesh(filename):
    r = vtkXMLUnstructuredGridReader()
    r.SetFileName(filename)
    r.Update()
    return r.GetOutput()

def writeMesh(mesh, filename):
    w = vtkXMLUnstructuredGridWriter()
    w.SetFileName(filename)
    w.SetInputData(mesh)
    w.SetDataModeToAscii()
    w.Update()

def getIntegrationPointArrays(mesh):
    offsets = mesh.GetCellData().GetArray('integration_point_offsets')
    if offsets == None:
        print("Error: integration_point_offsets cell data array not found.")
        sys.exit(2)
    ip_data = mesh.GetFieldData().GetArray('integration_point_data')
    if ip_data == None:
        print("Error: integration_point_data field data array not found.")
        sys.exit(2)
    return [offsets, ip_data]

def processIntegrationPointData(filename):
    m = readMesh(filename)

    [offsets, ip_data] = getIntegrationPointArrays(m)

    element_data = readElementData(ip_data, offsets)

    np_sigma = np.empty([m.GetNumberOfCells(), 4])

    for ed in element_data:
        np_sigma[ed.id] = averageOverIntegrationPoints(ed)

    sigma = numpy_to_vtk(np_sigma);
    sigma.SetName("sigma")
    m.GetCellData().AddArray(sigma)

    writeMesh(m, filename)

def main():
    if len(sys.argv) != 2:
        print('Expected input file name.')
        sys.exit(1)

    processIntegrationPointData(sys.argv[1])

if __name__ == '__main__':
    main()


def garbage():
    ### Create a protobuf object.
    a = ip.KelvinVector()
    a.value.append(5)
    a.value.append(3)

    ### Create an object from a string.
    ed = ip.ElementData()
    ed.ParseFromString('aeouaoeuaoeu')


