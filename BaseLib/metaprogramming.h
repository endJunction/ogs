/**
 * \copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef BASELIB_METAPROGRAMMING_H
#define BASELIB_METAPROGRAMMING_H


namespace BaseLib
{

template <typename T>
inline constexpr T pow(T const x, unsigned const y)
{
	return (y == 0) ? 1 : x * pow(x, y - 1);
};

}	// namespace BaseLib

#endif	// BASELIB_METAPROGRAMMING_H
