/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ThirdParty/iphreeqc/src/src/phreeqcpp/Phreeqc.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
struct PK : private Phreeqc
{
    void initializePhreeqcGeneralSettings();
    void loadDatabase(std::string const& database);
    void setConvergenceTolerance();
    void configureOutputSettings();
    void addRxnSolution(int const chemical_system_id,
                        cxxSolution const* const solution);
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
