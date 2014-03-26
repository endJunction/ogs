/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_SIMPLEASSEMBLER_H_
#define ASSEMBLERLIB_SIMPLEASSEMBLER_H_

namespace AssemblerLib
{

template<
    typename MESH_ITEM_,
    typename ASSEMBLER_>
class SimpleAssembler
{
public:
    SimpleAssembler(ASSEMBLER_ &assembler)
    : _assembler(assembler) {}

    ~SimpleAssembler() {}

    /// Executes given assembler for the given mesh item and a data item.
    /// \attention The index \c id is not necesserily the mesh item's id.
    template <typename ITEM_DATA>
    void operator()(const MESH_ITEM_* item, std::size_t id, ITEM_DATA& item_data) const
    {
        _assembler(*item, item_data);
    }

protected:
    ASSEMBLER_ &_assembler;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_SIMPLEASSEMBLER_H_
