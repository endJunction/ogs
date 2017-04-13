/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>
#include <boost/math/special_functions/pow.hpp>
#include <cmath>

#include "MathLib/Integration/GaussLegendre.h"
#include "MathLib/TemplateWeightedPoint.h"


namespace NumLib
{

/// Gauss quadrature rule for regular shape elements: line, quad and hex.
///
/// \tparam N_DIM    Spatial dimension
template <unsigned N_DIM>
class IntegrationGaussRegular
{
    using WeightedPoint =
        typename MathLib::TemplateWeightedPoint<double, double, N_DIM>;

public:
    /// Create IntegrationGaussRegular of the given Gauss-Legendre integration
    /// order.
    ///
    /// @param order     integration order (default 2)
    explicit IntegrationGaussRegular(unsigned order = 2)
    : _order(order), _n_sampl_pt(0)
    {
        this->setIntegrationOrder(order);
    }

    /// Change the integration order.
    void setIntegrationOrder(unsigned order)
    {
        this->_n_sampl_pt =
            static_cast<unsigned>(boost::math::pow<N_DIM>(order));
        this->_order = order;
    }

    /// return current integration order.
    unsigned getIntegrationOrder() const {return _order;}

    /// return the number of sampling points
    unsigned getNumberOfPoints() const {return _n_sampl_pt;}

    /// Get coordinates of the integration point.
    ///
    /// @param igp       The integration point index
    /// @return a weighted point
    WeightedPoint
    getWeightedPoint(unsigned igp) const
    {
        return getWeightedPoint(getIntegrationOrder(), igp);
    }

    /// Get position indexes in r-s-t axis.
    ///
    /// @param order     The number of integration points
    /// @param igp       The integration point index
    /// @return  a tuple of position indexes
    static std::array<unsigned, N_DIM> getPositionIndices(unsigned order, unsigned igp);

    /// Get coordinates of the integration point.
    ///
    /// @param order     The number of integration points
    /// @param igp       The integration point index
    /// @return a weighted point
    static WeightedPoint
    getWeightedPoint(unsigned order, unsigned igp);

private:
    /// Computes weighted point using given integration method.
    ///
    /// \tparam Method  Integration method to use.
    /// \param  pos     Point indices computed by getPositionIndices.
    template <typename Method>
    static
    WeightedPoint
    getWeightedPoint(std::array<unsigned, N_DIM> const& pos);

private:
    unsigned _order;
    unsigned _n_sampl_pt;
};

} // NumLib

#include "IntegrationGaussRegular-impl.h"
