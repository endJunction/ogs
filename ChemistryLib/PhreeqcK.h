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

#include <string>

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
struct PK;
class PhreeqcK
{
public:
    void initializePhreeqcGeneralSettings();
    void loadDatabase(std::string const& database);
    void setConvergenceTolerance();
    void configureOutputSettings();

private:
    PK* impl;
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
