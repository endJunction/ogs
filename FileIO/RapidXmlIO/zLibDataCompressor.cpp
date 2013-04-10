/**
 * \file
 * \author Karsten Rink
 * \date   2012-11-26
 * \brief  Implementation of the zLibDataCompressor class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Based on the vtkZLibDataCompressor-class in VTK 5.6
 */

#include <cstddef>
#include <iostream>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "zlib/zlib.h"

#include "zLibDataCompressor.h"

void printZlibErrorMessage(int const error)
{
	ERR("zLibDataCompressor::UncompressBuffer(): Zlib error while uncompressing data:");

	switch (error)
	{
		case Z_MEM_ERROR:    ERR("Not enought memory."); break;
		case Z_BUF_ERROR:    ERR("Not enough room in the output buffer."); break;
		case Z_DATA_ERROR:   ERR("Input data corrupted or incomplete."); break;
		case Z_STREAM_ERROR: ERR("Some stream state is inconsistent. zlib returned Z_STREAM_ERROR."); break;
		default:             ERR("Unknown zlib result: %d.", error);
	}
}

unsigned long zLibDataCompressor::CompressBuffer(const unsigned char* uncompressedData,
                                                 unsigned long uncompressedSize,
                                                 unsigned char* compressedData,
                                                 unsigned long compressionSpace)
{
	int CompressionLevel = Z_DEFAULT_COMPRESSION;
	unsigned long compressedSize = compressionSpace;
	Bytef* cd = reinterpret_cast<Bytef*>(compressedData);
	const Bytef* ud = reinterpret_cast<const Bytef*>(uncompressedData);

	// zlib call
	int const zlibResult = compress2(cd, &compressedSize, ud, uncompressedSize, CompressionLevel);
	if (zlibResult != Z_OK)
	{
		printZlibErrorMessage(zlibResult);
		return 0;
	}

	return compressedSize;
}

unsigned long zLibDataCompressor::UncompressBuffer(const unsigned char* compressedData,
                                                   unsigned long compressedSize,
                                                   unsigned char* uncompressedData,
                                                   unsigned long uncompressedSize)
{
	unsigned long decSize = uncompressedSize;
	Bytef* ud = reinterpret_cast<Bytef*>(uncompressedData);
	const Bytef* cd = reinterpret_cast<const Bytef*>(compressedData);

	// zlib call
	int const zlibResult = uncompress(ud, &decSize, cd, compressedSize);
	if (zlibResult != Z_OK)
	{
		printZlibErrorMessage(zlibResult);
		return 0;
	}

	// Make sure the output size matched that expected.
	if(decSize != uncompressedSize)
	{
		WARN("Decompression produced incorrect size. Expected %d and got %d.",
		     uncompressedSize, decSize);
		return 0;
	}

	return decSize;
}
