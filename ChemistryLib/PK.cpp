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

}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
