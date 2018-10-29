/**
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "mpProperty.h"
#include "Properties/properties.h"
#include "mpComponent.h"
#include "mpMedium.h"
#include "mpPhase.h"

#include <sstream>
#include <string>
#include <iostream>
#include <vector>

namespace MaterialPropertyLib
{

PropertyDataType Property::value() const
{
    return _value;
}
/// The default implementation of this method only returns the
/// property value without altering it.
PropertyDataType Property::value(VariableArray const& /*unused*/)
{
    return _value;
}

/// Default implementation: derivative of any constant property is zero.
PropertyDataType Property::dvalue(VariableArray const&, Variables const)
{
    return _dvalue;
}

/// Default implementation: 2nd derivative of any constant property is zero.
PropertyDataType Property::ddvalue(VariableArray const&, Variables const pv1,
        Variables const pv2)
{
    return 0.0;
}

void Property::notImplemented(std::string property, std::string material)
{
    OGS_FATAL(
        "The property \'%s\' is not "
        "available on the \'%s\' scale",
        property.c_str(), material.c_str());
}

std::unique_ptr<Property> newProperty(BaseLib::ConfigTree const& config,
                                      Medium* m)
{
    return selectProperty(config, m);
}

std::unique_ptr<Property> newProperty(BaseLib::ConfigTree const& config,
                                      Phase* p)
{
    return selectProperty(config, p);
}

std::unique_ptr<Property> newProperty(BaseLib::ConfigTree const& config,
                                      Component* c)
{
    return selectProperty(config, c);
}

/**
 * This method creates a new Property. It uses the information stored
 * inside the configuration tree to select the property that is specified
 * by the user. It further passes a pointer to the material that requests
 * the property.
*/
template <typename MaterialType>
std::unique_ptr<Property> selectProperty(BaseLib::ConfigTree const& config,
                                         MaterialType M)
{
    // Parsing the property type:
    auto const property_type = config.getConfigParameter<std::string>("type");

    // If (and only if) the given property type is 'constant', a
    // corresponding value is needed.
    if (boost::iequals(property_type, "constant"))
    {
        /// \todo: Creating constant properties from prj-file is
        /// currently restricted to scalar values. The Constant
        /// constructor, however, can handle any datatype defined
        /// by PropertyDataType. This could be enhanced in order
        /// to define vectors or even tensors as constant properties.

        auto const sValue = config.getConfigParameter<std::string>("value");

        std::stringstream issValue(sValue);

        std::vector<double> values;
        double dummy(0);

        while (issValue >> dummy)
            values.push_back(dummy);

        switch (values.size())
        {
        case 1:
        {
            // scalar
            PropertyDataType property_value = values[0];
            std::cout << "---create a constant scalar value property.\n";

            return std::make_unique<Constant>(Constant(property_value));
        }
        case 2:
        {
            // Pair
            PropertyDataType property_value = (Pair){values[0], values[1]};
            std::cout << "---create a constant paired value property.\n";
            return std::make_unique<Constant>(Constant(property_value));
        }
        case 3:
        {
            // Vector
            PropertyDataType property_value =
                    (Vector){values[0], values[1], values[2]};
            std::cout << "---create a constant vector value property.\n";
            return std::make_unique<Constant>(Constant(property_value));
        }
        case 6:
        {
            // Symmetric Tensor - xx, yy, zz, xy, xz, yz
            PropertyDataType property_value =
                    (SymmTensor){values[0], values[1], values[2],
                                 values[3], values[4], values[5]};
            std::cout << "---create a constant symmetric tensor value property.\n";
            return std::make_unique<Constant>(Constant(property_value));
        }
        case 9:
        {
            // Tensor
            PropertyDataType property_value =
                    (Tensor){values[0], values[1], values[2],
                             values[3], values[4], values[5],
                             values[6], values[7], values[8]};
            std::cout << "---create a constant tensor value property.\n";
            return std::make_unique<Constant>(Constant(property_value));
        }

        default:
        {
            OGS_FATAL ("Given number of components for constant property. /%i", values.size());
        }
        }

        PropertyDataType property_value;
        return std::make_unique<Constant>(Constant(property_value));

       // return std::make_unique<Constant>(Constant(property_value));
    }
    /// Properties can be medium, phase, or component properties.
    /// Some of them require information about the respective material.
    /// Thus, we pass a pointer to the material that requests the property.
    /// In this method, this pointer is realized via typename MaterialType,
    /// which replaces either Medium*, Phase*, or Component*.
    /// Note that most property constructors (only those that request material
    /// pointers) must be overloaded for any type of material.
    if (boost::iequals(property_type, "BilinearTemperaturePressure"))
    {
        return createBilinearTemperaturePressure(config, M);
    }
    if (boost::iequals(property_type, "LinearTemperature"))
    {
        return std::make_unique<LinearTemperature>(M);
    }
    if (boost::iequals(property_type, "LinearEpsilon"))
    {
        OGS_FATAL("TODO: Implementation of %s property!!", property_type.c_str());
    }
    if (boost::iequals(property_type, "AverageMolarMass"))
    {
        OGS_FATAL("TODO: Implementation of %s property!!", property_type.c_str());
    }
    if (boost::iequals(property_type, "AverageVolumeFraction"))
    {
        return std::make_unique<AverageMoleFraction>(M);
    }
    if (boost::iequals(property_type, "Duan_2012"))
    {
        OGS_FATAL("TODO: Implementation of %s property!!", property_type.c_str());
    }
    if (boost::iequals(property_type, "IAPWS_2008"))
    {
        return std::make_unique<ViscosityWaterIAPWS>(M);
    }
    if (boost::iequals(property_type, "ideal_gas_law"))
    {
        return std::make_unique<IdealGasLaw>(M);
    }
    if (boost::iequals(property_type, "Islam_Carlson_2012"))
    {
        OGS_FATAL("TODO: Implementation of %s property!!", property_type.c_str());
    }
    if (boost::iequals(property_type, "Fenghour_1998"))
    {
        return std::make_unique<ViscosityCO2Fenghour>(M);
    }
    if (boost::iequals(property_type, "Mualem_1978"))
    {
        OGS_FATAL("TODO: Implementation of %s property!!", property_type.c_str());
    }
    if (boost::iequals(property_type, "Buddenberg_Wilke_1949"))
    {
        OGS_FATAL("TODO: Implementation of %s property!!", property_type.c_str());
    }
    if (boost::iequals(property_type, "Peng_Robinson_1976"))
    {
        return std::make_unique<PengRobinson>(M);
    }
    if (boost::iequals(property_type, "Brooks_Corey_1964_saturation"))
    {
        return std::make_unique<BrooksCoreySaturation>(M);
    }
    if (boost::iequals(property_type, "RelPermBrooksCorey"))
    {
        return createRelPermBrooksCorey(config, M);
    }
    if (boost::iequals(property_type, "RelPermLiakopoulos"))
    {
        return std::make_unique<RelPermLiakopoulos>(M);
    }
    if (boost::iequals(property_type, "SaturationLiakopoulos"))
    {
        return std::make_unique<SaturationLiakopoulos>(M);
    }
    if (boost::iequals(property_type, "SaturationFredlund"))
    {
        return createSaturationFredlund(config, M);
    }
    // If none of the above property types are found, OGS throws an error.
    OGS_FATAL(
        "The specified component property type \"%s\" was not "
        "recognized",
        property_type.c_str());
    //return nullptr;  // to avoid 'no return' warnings
}

}  // namespace MaterialPropertyLib
