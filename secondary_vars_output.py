# coding: utf-8
import sys
import numpy as np

from vtk import *
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

import integration_point_pb2 as ip
import common_pb2
import material_pb2


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

def averageOverIntegrationPointsSigma(element_data):
    if element_data.n_integration_points == 0:
        print("Error: Got 0 integration points.")
        sys.exit(2)

    sd = element_data.small_deformation
    result = np.array([0.0, 0.0, 0.0, 0.0])
    for i in range(element_data.n_integration_points):
        result += sd.sigma[i].value
    return result / element_data.n_integration_points

def averageOverIntegrationPointsEpsPD(element_data):
    if element_data.n_integration_points == 0:
        print("Error: Got 0 integration points.")
        sys.exit(2)

    sd = element_data.small_deformation
    result = np.array([0.0, 0.0, 0.0, 0.0])
    for i in range(element_data.n_integration_points):
        result += sd.material_state[i].ehlers.eps_p_D.value

    return result / element_data.n_integration_points

def averageOverIntegrationPointsEpsPV(element_data):
    if element_data.n_integration_points == 0:
        print("Error: Got 0 integration points.")
        sys.exit(2)

    sd = element_data.small_deformation
    result = 0;
    for i in range(element_data.n_integration_points):
        result += sd.material_state[i].ehlers.eps_p_V
    return result / element_data.n_integration_points

def averageOverIntegrationPointsEpsPEff(element_data):
    if element_data.n_integration_points == 0:
        print("Error: Got 0 integration points.")
        sys.exit(2)

    sd = element_data.small_deformation
    result = 0;
    for i in range(element_data.n_integration_points):
        result += sd.material_state[i].ehlers.eps_p_eff
    return result / element_data.n_integration_points

def averageOverIntegrationPointsKappaD(element_data):
    if element_data.n_integration_points == 0:
        print("Error: Got 0 integration points.")
        sys.exit(2)

    sd = element_data.small_deformation
    result = 0;
    for i in range(element_data.n_integration_points):
        result += sd.material_state[i].ehlers.kappa_d
    return result / element_data.n_integration_points

def averageOverIntegrationPointsDamage(element_data):
    if element_data.n_integration_points == 0:
        print("Error: Got 0 integration points.")
        sys.exit(2)

    sd = element_data.small_deformation
    result = 0;
    for i in range(element_data.n_integration_points):
        result += sd.material_state[i].ehlers.damage
    return result / element_data.n_integration_points

def processIntegrationPointData(filename):
    m = readMesh(filename)

    [offsets, ip_data] = getIntegrationPointArrays(m)

    element_data = readElementData(ip_data, offsets)

    np_sigma = np.empty([m.GetNumberOfCells(), 4])
    np_eps_p_D = np.empty([m.GetNumberOfCells(), 4])
    np_eps_p_V = np.empty([m.GetNumberOfCells(), 1])
    np_eps_p_eff = np.empty([m.GetNumberOfCells(), 1])
    np_kappa_d = np.empty([m.GetNumberOfCells(), 1])
    np_damage = np.empty([m.GetNumberOfCells(), 1])

    for ed in element_data:
        np_sigma[ed.element_id] = averageOverIntegrationPointsSigma(ed)
        np_eps_p_D[ed.element_id] = averageOverIntegrationPointsEpsPD(ed)
        np_eps_p_V[ed.element_id] = averageOverIntegrationPointsEpsPV(ed)
        np_eps_p_eff[ed.element_id] = averageOverIntegrationPointsEpsPEff(ed)
        np_kappa_d[ed.element_id] = averageOverIntegrationPointsKappaD(ed)
        np_damage[ed.element_id] = averageOverIntegrationPointsDamage(ed)

    sigma = numpy_to_vtk(np_sigma);
    sigma.SetName("sigma")
    m.GetCellData().AddArray(sigma)

    eps_p_D = numpy_to_vtk(np_eps_p_D);
    eps_p_D.SetName("eps_p_D")
    m.GetCellData().AddArray(eps_p_D)

    eps_p_V = numpy_to_vtk(np_eps_p_V);
    eps_p_V.SetName("eps_p_V")
    m.GetCellData().AddArray(eps_p_V)

    eps_p_eff = numpy_to_vtk(np_eps_p_eff);
    eps_p_eff.SetName("eps_p_eff")
    m.GetCellData().AddArray(eps_p_eff)

    kappa_d = numpy_to_vtk(np_kappa_d);
    kappa_d.SetName("kappa_d")
    m.GetCellData().AddArray(kappa_d)

    damage = numpy_to_vtk(np_damage);
    damage.SetName("damage")
    m.GetCellData().AddArray(damage)

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

def generateIntegrationPoints(mesh):
    qp_dict_gen = vtkQuadratureSchemeDictionaryGenerator()
    qp_dict_gen.SetInputData(m)
    qp_dict_gen.Update()
    x = qp_dict_gen.GetOutput()

    qp_gen = vtkQuadraturePointsGenerator()
    qp_gen.SetInputConnection(qp_dict_gen.GetOutputPort())
    qp_gen.SetInputArrayToProcess(0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, "QuadratureOffset");
    qp_gen.Update()
    qp_gen.GetOutput()

    w = vtkXMLPolyDataWriter()
    w.SetFileName("/tmp/x.vtp")
    w.SetInputConnection(qp_gen.GetOutputPort())
    w.Update()

