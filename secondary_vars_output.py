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
        print(
            'Error: integration point offsets array has wrong number of components.'
        )
        print('Expected one component.')
        sys.exit(2)
    if offsets.GetDataType() != VTK_UNSIGNED_LONG:
        print('Error: integration point offsets array has wrong data type.')
        print('Expected VTK_UNSIGNED_LONG.')
        sys.exit(2)

    if ip_data.GetNumberOfComponents() != 1:
        print(
            'Error: integration point data array has wrong number of components.'
        )
        print('Expected one component.')
        sys.exit(2)
    if ip_data.GetDataType() != VTK_CHAR:
        print('Error: integration point data array has wrong data type.')
        print('Expected VTK_CHAR.')
        sys.exit(2)

    ip_data_array = vtk_to_numpy(ip_data)

    element_data = []  # Resulting array

    # Process all but the last one.
    for i in range(offsets.GetNumberOfTuples() - 1):
        start = int(offsets.GetTuple1(i))
        end = int(offsets.GetTuple1(i + 1))
        ed = ip.ElementData()
        ed.ParseFromString(ip_data_array[start:end].tobytes())
        element_data.append(ed)
    # Process the last element
    start = int(offsets.GetTuple1(offsets.GetNumberOfTuples() - 1))
    end = ip_data.GetNumberOfTuples()
    if start >= end:
        print('Error: integration point offsets and integration point data \
               array sizes do not match.\nstart', start, ', end', end)
        sys.exit(2)
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


def averageOverIntegrationPointsSigma(sd_common, n_ip):
    assert len(sd_common.sigma) == n_ip
    result = np.zeros(len(sd_common.sigma[0].value))
    for i in range(n_ip):
        result += sd_common.sigma[i].value
    return result / n_ip


def averageOverIntegrationPointsEpsPD(sd_common, n_ip):
    assert len(sd_common.sigma) == n_ip
    result = np.zeros(len(sd_common.sigma[0].value))
    for i in range(n_ip):
        result += sd_common.material_state[i].ehlers.eps_p_D.value

    return result / n_ip


def averageOverIntegrationPointsEpsPV(sd_common, n_ip):
    result = 0
    for i in range(n_ip):
        result += sd_common.material_state[i].ehlers.eps_p_V
    return result / n_ip


def averageOverIntegrationPointsEpsPEff(sd_common, n_ip):
    result = 0
    for i in range(n_ip):
        result += sd_common.material_state[i].ehlers.eps_p_eff
    return result / n_ip


def averageOverIntegrationPointsKappaD(sd_common, n_ip):
    result = 0
    for i in range(n_ip):
        result += sd_common.material_state[i].ehlers.kappa_d
    return result / n_ip


def averageOverIntegrationPointsDamage(sd_common, n_ip):
    result = 0
    for i in range(n_ip):
        result += sd_common.material_state[i].ehlers.damage
    return result / n_ip


def averageOverIntegrationPointsNonlocalDamage(sdn):
    return np.average(sdn.nonlocal_damage)


def averageOverIntegrationPointsNonlocalLength(sdn):
    return np.average(sdn.nonlocal_length)


def processIntegrationPointData(filename):
    m = readMesh(filename)

    [offsets, ip_data] = getIntegrationPointArrays(m)

    element_data = readElementData(ip_data, offsets)

    assert m.GetNumberOfCells() > 0
    mesh_dimension = m.GetCell(0).GetCellDimension()
    kelvin_vector_size = [0, 0, 4, 6][mesh_dimension]
    assert kelvin_vector_size > 0

    np_sigma = np.empty([m.GetNumberOfCells(), kelvin_vector_size])
    np_eps_p_D = np.empty([m.GetNumberOfCells(), kelvin_vector_size])
    np_eps_p_V = np.empty([m.GetNumberOfCells(), 1])
    np_eps_p_eff = np.empty([m.GetNumberOfCells(), 1])
    np_kappa_d = np.empty([m.GetNumberOfCells(), 1])
    np_damage = np.empty([m.GetNumberOfCells(), 1])

    # only used if SDN is present
    has_sdn = False
    np_nonlocal_damage = np.empty([m.GetNumberOfCells(), 1])
    np_nonlocal_length = np.empty([m.GetNumberOfCells(), 1])

    for ed in element_data:
        n_ip = ed.n_integration_points
        if not n_ip > 0:
            print("Error: Got", n_ip,
                  "integration points, but positive number expected")
            sys.exit(2)

        ipdata = getattr(ed, ed.WhichOneof('ipdata'))  # one of SD or SDN
        c = ipdata.common

        np_sigma[ed.element_id] = averageOverIntegrationPointsSigma(c, n_ip)
        np_eps_p_D[ed.element_id] = averageOverIntegrationPointsEpsPD(c, n_ip)
        np_eps_p_V[ed.element_id] = averageOverIntegrationPointsEpsPV(c, n_ip)
        np_eps_p_eff[ed.element_id] = averageOverIntegrationPointsEpsPEff(
            c, n_ip)
        np_kappa_d[ed.element_id] = averageOverIntegrationPointsKappaD(c, n_ip)
        np_damage[ed.element_id] = averageOverIntegrationPointsDamage(c, n_ip)

        if ipdata.DESCRIPTOR.full_name == 'OGS.SmallDeformationNonlocal':
            has_sdn = True
            np_nonlocal_damage[
                ed.element_id] = averageOverIntegrationPointsNonlocalDamage(
                    ipdata)
            np_nonlocal_length[
                ed.element_id] = averageOverIntegrationPointsNonlocalLength(
                    ipdata)

    sigma = numpy_to_vtk(np_sigma)
    sigma.SetName("sigma")
    m.GetCellData().AddArray(sigma)

    eps_p_D = numpy_to_vtk(np_eps_p_D)
    eps_p_D.SetName("eps_p_D")
    m.GetCellData().AddArray(eps_p_D)

    eps_p_V = numpy_to_vtk(np_eps_p_V)
    eps_p_V.SetName("eps_p_V")
    m.GetCellData().AddArray(eps_p_V)

    eps_p_eff = numpy_to_vtk(np_eps_p_eff)
    eps_p_eff.SetName("eps_p_eff")
    m.GetCellData().AddArray(eps_p_eff)

    kappa_d = numpy_to_vtk(np_kappa_d)
    kappa_d.SetName("kappa_d")
    m.GetCellData().AddArray(kappa_d)

    damage = numpy_to_vtk(np_damage)
    damage.SetName("damage")
    m.GetCellData().AddArray(damage)

    if has_sdn:
        nonlocal_damage = numpy_to_vtk(np_nonlocal_damage)
        nonlocal_damage.SetName("nonlocal_damage")
        m.GetCellData().AddArray(nonlocal_damage)

        nonlocal_length = numpy_to_vtk(np_nonlocal_length)
        nonlocal_length.SetName("nonlocal_length")
        m.GetCellData().AddArray(nonlocal_length)

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
    qp_gen.SetInputArrayToProcess(
        0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, "QuadratureOffset")
    qp_gen.Update()
    qp_gen.GetOutput()

    w = vtkXMLPolyDataWriter()
    w.SetFileName("/tmp/x.vtp")
    w.SetInputConnection(qp_gen.GetOutputPort())
    w.Update()
