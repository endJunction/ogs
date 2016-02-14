/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROPERTYVECTOR_H_
#define PROPERTYVECTOR_H_

#include <iterator>
#include <ostream>
#include <string>
#include <vector>

#include "BaseLib/excludeObjectCopy.h"
#include "Location.h"

namespace MeshLib
{

//template <typename T> class PropertyVector;

class PropertyVectorBase
{
public:
	virtual PropertyVectorBase* clone(
		std::vector<std::size_t> const& exclude_positions
	) const = 0;
	virtual ~PropertyVectorBase() = default;
};

/// Class template PropertyVector is a std::vector with template parameter
/// PROP_VAL_TYPE. The reason for the derivation of std::vector is
/// the template specialisation for pointer types below.
/// \tparam PROP_VAL_TYPE typical this is a scalar, a vector or a matrix
template <typename PROP_VAL_TYPE>
class PropertyVector : public std::vector<PROP_VAL_TYPE>,
	public PropertyVectorBase
{
friend class Properties;

public:
	std::size_t getNumberOfComponents() const { return _n_components; }
	std::size_t getNumberOfTuples() const
	{
		return std::vector<PROP_VAL_TYPE>::size() / _n_components;
	}
	MeshItemType getMeshItemType() const { return _mesh_item_type; }
	std::string const& getPropertyName() const { return _property_name; }

	PropertyVectorBase* clone(std::vector<std::size_t> const& exclude_positions) const
	{
		PropertyVector<PROP_VAL_TYPE> *t(new PropertyVector<PROP_VAL_TYPE>(_property_name,
			_mesh_item_type, _n_components));
		BaseLib::excludeObjectCopy(*this, exclude_positions, *t);
		return t;
	}

	/// Method returns the number of tuples times the number of tuple components.
	std::size_t size() const
	{
		return std::vector<PROP_VAL_TYPE>::size();
	}

protected:
	/// @brief The constructor taking meta information for the data.
	/// @param property_name a string describing the property
	/// @param mesh_item_type the values of the property are either assigned to
	/// nodes or cells (see enumeration MeshItemType)
	/// @param n_components the number of elements of a tuple
	explicit PropertyVector(std::string const& property_name,
	                        MeshItemType mesh_item_type,
	                        std::size_t n_components)
	    : std::vector<PROP_VAL_TYPE>(),
	      _n_components(n_components),
	      _mesh_item_type(mesh_item_type),
	      _property_name(property_name)
	{}

	/// @brief The constructor taking meta information for the data.
	/// @param n_property_values number of property values (value can be a tuple
	/// with several entries)
	/// @param property_name a string describing the property
	/// @param mesh_item_type the values of the property are either assigned to
	/// nodes or cells (see enumeration MeshItemType)
	/// @param n_components the number of elements of a tuple
	PropertyVector(std::size_t n_property_values,
	               std::string const& property_name,
	               MeshItemType mesh_item_type,
	               std::size_t n_components)
	    : std::vector<PROP_VAL_TYPE>(n_property_values * n_components),
	      _mesh_item_type(mesh_item_type),
	      _property_name(property_name)
	{}

	std::size_t const _n_components;
	MeshItemType const _mesh_item_type;
	std::string const _property_name;
};

/// Class template PropertyVector is a std::vector with template parameter
/// T, where T is a pointer type.
/// The behaviour has changed for the constructor, destructor and the operator[].
/// The user has to provide the size and an item to group mapping for construction.
/// The destructor takes care to delete the entries of the vector.
/// The operator[] uses an item-to-group property map to access the
/// correct property.
/// \tparam T pointer type, the type the type points to is typical a scalar,
/// a vector or a matrix type
template <typename T>
class PropertyVector<T*> : public std::vector<std::size_t>,
	public PropertyVectorBase
{
friend class Properties;
public:
	/// Destructor ensures the deletion of the heap-constructed objects.
	~PropertyVector()
	{
		for (auto v : _values)
			delete v;
	}

	/// The operator[] uses the item to group property map to access to the
	/// correct property value/object.
	T* const& operator[](std::size_t id) const
	{
		return _values[std::vector<std::size_t>::operator[](id)];
	}

	T* & operator[](std::size_t id)
	{
		return _values[std::vector<std::size_t>::operator[](id)];
	}

	void initPropertyValue(std::size_t group_id, T const& value)
	{
		_values[group_id] = new T(value);
	}

	std::size_t getNumberOfComponents() const { return _n_components; }
	std::size_t getNumberOfTuples() const
	{
		return std::vector<std::size_t>::size();
	}
	/// Method returns the number of tuples times the number of tuple components.
	std::size_t size() const
	{
		return _n_components * std::vector<std::size_t>::size();
	}
	MeshItemType getMeshItemType() const { return _mesh_item_type; }
	std::string const& getPropertyName() const { return _property_name; }

	PropertyVectorBase* clone(std::vector<std::size_t> const& exclude_positions) const
	{
		// create new PropertyVector with modified mapping
		PropertyVector<T*> *t(new PropertyVector<T*>
			(
				_values.size()/_n_components,
				BaseLib::excludeObjectCopy(*this, exclude_positions),
				_property_name, _mesh_item_type, _n_components
			)
		);
		// copy pointers to property values
		for (std::size_t j(0); j<_values.size(); j++) {
			t->initPropertyValue(j, *(_values[j]));
		}
		return t;
	}

#ifndef NDEBUG
	std::ostream& print(std::ostream &os) const
	{
		os << "\nPropertyVector<T*> at address: " << this << ":\n";
		os << "\tmapping (" << size() <<"):\n";
		std::copy(this->cbegin(), this->cend(),
			std::ostream_iterator<std::size_t>(os, " "));
		os << "\n\tvalues (" << _values.size() << "):\n";
		for (std::size_t k(0); k<_values.size(); k++) {
			os << "val: " << *(_values[k]) << ", address: " << _values[k] << "\n";
		}
		return os;
	}
#endif

protected:
	/// @brief The constructor taking meta information for the data.
	/// @param n_prop_groups number of different property values
	/// @param item2group_mapping Class Mesh has a mapping from the mesh items
	/// (Node or Element) to an index (position in the data structure).
	/// The vector item2group_mapping must have the same number of entries as
	/// the above mapping and the values have to be in the range
	/// \f$[0, \text{n_prop_groups})\f$.
	/// @param property_name a string describing the property
	/// @param mesh_item_type the values of the property are either assigned to
	/// nodes or cells (see enumeration MeshItemType)
	/// @param n_components the number of elements of a tuple
	PropertyVector(std::size_t n_prop_groups,
	               std::vector<std::size_t> const& item2group_mapping,
	               std::string const& property_name,
	               MeshItemType mesh_item_type,
	               std::size_t n_components)
	    : std::vector<std::size_t>(item2group_mapping),
	      _n_components(n_components),
	      _mesh_item_type(mesh_item_type),
	      _property_name(property_name),
	      _values(n_prop_groups * n_components)
	{}

protected:
	std::size_t const _n_components;
	MeshItemType const _mesh_item_type;
	std::string const _property_name;

private:
	std::vector<T*> _values;
	// hide method
	T* at(std::size_t);
};

} // end namespace MeshLib

#endif
