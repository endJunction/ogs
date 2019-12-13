/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PK.h"

#include <iostream>

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
void PK::initializePhreeqcGeneralSettings()
{
    do_initialize();
}

void PK::loadDatabase(std::string const& database)
{
    std::ifstream in(database);
    if (!in)
    {
        std::cerr << "Unable to open database file " << database;
        std::terminate();
    }
    assert(phrq_io->get_istream() == nullptr);
    phrq_io->push_istream(&in, false);
    read_database();
}

void PK::setConvergenceTolerance()
{
    convergence_tolerance = 1e-12;
}

void PK::configureOutputSettings()
{
    pr.all = false;
}

void PK::addRxnSolution(int const chemical_system_id,
                        cxxSolution const* const solution)
{
    Rxn_solution_map.emplace(chemical_system_id, *solution);
}

}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
