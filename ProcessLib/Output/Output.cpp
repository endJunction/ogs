/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Output.h"

#include <cassert>
#include <fstream>
#include <vector>

#include <logog/include/logog.hpp>

#include "Applications/InSituLib/Adaptor.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "ProcessLib/Process.h"

namespace
{
//! Converts a vtkXMLWriter's data mode string to an int. See
/// Output::_output_file_data_mode.
int convertVtkDataMode(std::string const& data_mode)
{
    if (data_mode == "Ascii")
    {
        return 0;
    }
    if (data_mode == "Binary")
    {
        return 1;
    }
    if (data_mode == "Appended")
    {
        return 2;
    }
    OGS_FATAL(
        "Unsupported vtk output file data mode \"%s\". Expected Ascii, "
        "Binary, or Appended.",
        data_mode.c_str());
}
}  // namespace

namespace ProcessLib
{
bool Output::shallDoOutput(unsigned timestep, double const t)
{
    unsigned each_steps = 1;

    for (auto const& pair : _repeats_each_steps)
    {
        each_steps = pair.each_steps;

        if (timestep > pair.repeat * each_steps)
        {
            timestep -= pair.repeat * each_steps;
        }
        else
        {
            break;
        }
    }

    bool make_output = timestep % each_steps == 0;

    if (_fixed_output_times.empty())
        return make_output;

    const double specific_time = _fixed_output_times.back();
    const double zero_threshold = std::numeric_limits<double>::min();
    if (std::fabs(specific_time - t) < zero_threshold)
    {
        _fixed_output_times.pop_back();
        make_output = true;
    }

    return make_output;
}

Output::Output(std::string output_directory, std::string prefix,
               bool const compress_output, std::string const& data_mode,
               bool const output_nonlinear_iteration_results,
               std::vector<PairRepeatEachSteps> repeats_each_steps,
               std::vector<double>&& fixed_output_times,
               std::vector<MathLib::Point3d>&& timeseries_output_points)
    : _output_directory(std::move(output_directory)),
      _output_file_prefix(std::move(prefix)),
      _output_file_compression(compress_output),
      _output_file_data_mode(convertVtkDataMode(data_mode)),
      _output_nonlinear_iteration_results(output_nonlinear_iteration_results),
      _repeats_each_steps(std::move(repeats_each_steps)),
      _fixed_output_times(std::move(fixed_output_times)),
      _timeseries_output_points(std::move(timeseries_output_points))
{
    if (!_timeseries_output_points.empty())
    {
        // Create output file.
        auto const filename =
            BaseLib::joinPaths(_output_directory, _output_file_prefix + ".csv");
        _timeseries_output_stream.open(filename.c_str(),
                                       std::ios_base::out);
        _timeseries_output_stream.precision(
            std::numeric_limits<double>::digits10);
    }
}

void Output::addProcess(ProcessLib::Process const& process,
                        const int process_id)
{
    auto const filename = BaseLib::joinPaths(
        _output_directory,
        _output_file_prefix + "_pcs_" + std::to_string(process_id) + ".pvd");
    _process_to_process_data.emplace(std::piecewise_construct,
                                     std::forward_as_tuple(&process),
                                     std::forward_as_tuple(filename));
}

// TODO return a reference.
Output::ProcessData* Output::findProcessData(Process const& process,
                                             const int process_id)
{
    auto range = _process_to_process_data.equal_range(&process);
    int counter = 0;
    ProcessData* process_data = nullptr;
    for (auto spd_it = range.first; spd_it != range.second; ++spd_it)
    {
        if (counter == process_id)
        {
            process_data = &spd_it->second;
            break;
        }
        counter++;
    }
    if (process_data == nullptr)
    {
        OGS_FATAL(
            "The given process is not contained in the output"
            " configuration. Aborting.");
    }

    return process_data;
}

void Output::writeTimeseriesData(
    GlobalVector const& x,
    double const t,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh)
{
    if (t == 0)
    {
        //
        // Find points in the mesh to get node_ids
        //
        _timeseries_output_node_ids.reserve(_timeseries_output_points.size());
        {
            auto search_length_algorithm =
                std::make_unique<MeshGeoToolsLib::SearchLength>(1e-8);

            MeshGeoToolsLib::MeshNodeSearcher const& mesh_node_searcher =
                MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
                    mesh, std::move(search_length_algorithm));

            std::vector<GeoLib::Point> search_points;
            search_points.reserve(_timeseries_output_points.size());
            std::transform(begin(_timeseries_output_points),
                           end(_timeseries_output_points),
                           std::back_inserter(search_points),
                           [](MathLib::Point3d const& p) {
                               return GeoLib::Point{
                                   p, std::numeric_limits<std::size_t>::max()};
                           });

            for (auto const& p : search_points)
            {
                std::vector<std::size_t> id =
                    mesh_node_searcher.getMeshNodeIDsForPoint(p);
                if (id.empty())
                    OGS_FATAL("Output point not found in the mesh.");
                if (id.size() > 1)
                    OGS_FATAL(
                        "Multiple output points found in the mesh but only one "
                        "is allowed.");
                DBUG("Found node %d for point %g, %g, %g.", id[0], p[0], p[1],
                     p[2]);
                _timeseries_output_node_ids.push_back(id[0]);
            }
        }
    }
    _timeseries_output_stream << t << "; ";

    //
    // Extract primary variables from the solution vector
    //
    {
        for (auto const& node_id : _timeseries_output_node_ids)
        {
            _timeseries_output_stream << node_id << " ";
            MeshLib::Location const l(mesh.getID(), MeshLib::MeshItemType::Node,
                                      node_id);
            auto global_indices = dof_table.getGlobalIndices(l);
            for (auto const& i : global_indices)
            {
                _timeseries_output_stream << x.get(i) << " ";
            }
            _timeseries_output_stream << ";";
        }
        _timeseries_output_stream << "\n";
    }
}

void Output::doOutputAlways(Process const& process,
                            const int process_id,
                            ProcessOutput const& process_output,
                            unsigned timestep,
                            const double t,
                            GlobalVector const& x)
{
    BaseLib::RunTime time_output;
    time_output.start();

    // Need to add variables of process to vtu even no output takes place.
    processOutputData(t, x, process.getMesh(), process.getDOFTable(process_id),
                      process.getProcessVariables(process_id),
                      process.getSecondaryVariables(),
                      process.getIntegrationPointWriter(),
                      process_output);

    // For the staggered scheme for the coupling, only the last process, which
    // gives the latest solution within a coupling loop, is allowed to make
    // output.
    if (!(process_id == static_cast<int>(_process_to_process_data.size()) - 1 ||
          process.isMonolithicSchemeUsed()))
        return;

    writeTimeseriesData(x, t, process.getDOFTable(process_id),
                        process.getMesh());

    std::string const output_file_name =
        _output_file_prefix + "_pcs_" + std::to_string(process_id) + "_ts_" +
        std::to_string(timestep) + "_t_" + std::to_string(t) + ".vtu";
    std::string const output_file_path =
        BaseLib::joinPaths(_output_directory, output_file_name);

    DBUG("output to %s", output_file_path.c_str());

    ProcessData* process_data = findProcessData(process, process_id);
    process_data->pvd_file.addVTUFile(output_file_name, t);
    INFO("[time] Output of timestep %d took %g s.", timestep,
         time_output.elapsed());

    makeOutput(output_file_path, process.getMesh(), _output_file_compression,
               _output_file_data_mode);
}

void Output::doOutput(Process const& process,
                      const int process_id,
                      ProcessOutput const& process_output,
                      unsigned timestep,
                      const double t,
                      GlobalVector const& x)
{
    if (shallDoOutput(timestep, t))
    {
        doOutputAlways(process, process_id, process_output, timestep, t, x);
    }
#ifdef USE_INSITU
    // Note: last time step may be output twice: here and in
    // doOutputLastTimestep() which throws a warning.
    InSituLib::CoProcess(process.getMesh(), t, timestep, false);
#endif
}

void Output::doOutputLastTimestep(Process const& process,
                                  const int process_id,
                                  ProcessOutput const& process_output,
                                  unsigned timestep,
                                  const double t,
                                  GlobalVector const& x)
{
    if (!shallDoOutput(timestep, t))
    {
        doOutputAlways(process, process_id, process_output, timestep, t, x);
    }
#ifdef USE_INSITU
    InSituLib::CoProcess(process.getMesh(), t, timestep, true);
#endif
}

void Output::doOutputNonlinearIteration(Process const& process,
                                        const int process_id,
                                        ProcessOutput const& process_output,
                                        const unsigned timestep, const double t,
                                        GlobalVector const& x,
                                        const unsigned iteration)
{
    if (!_output_nonlinear_iteration_results)
    {
        return;
    }

    BaseLib::RunTime time_output;
    time_output.start();

    processOutputData(t, x, process.getMesh(),
                      process.getDOFTable(process_id),
                      process.getProcessVariables(process_id),
                      process.getSecondaryVariables(),
                      process.getIntegrationPointWriter(),
                      process_output);

    // For the staggered scheme for the coupling, only the last process, which
    // gives the latest solution within a coupling loop, is allowed to make
    // output.
    if (!(process_id == static_cast<int>(_process_to_process_data.size()) - 1 ||
          process.isMonolithicSchemeUsed()))
        return;

    // Only check whether a process data is available for output.
    findProcessData(process, process_id);

    std::string const output_file_name =
        _output_file_prefix + "_pcs_" + std::to_string(process_id) + "_ts_" +
        std::to_string(timestep) + "_t_" + std::to_string(t) + "_nliter_" +
        std::to_string(iteration) + ".vtu";
    std::string const output_file_path =
        BaseLib::joinPaths(_output_directory, output_file_name);

    DBUG("output iteration results to %s", output_file_path.c_str());

    INFO("[time] Output took %g s.", time_output.elapsed());

    makeOutput(output_file_path, process.getMesh(), _output_file_compression,
               _output_file_data_mode);
}
}  // namespace ProcessLib
