/*****************************************************************************
 * Copyright (c) 2015 Adam Heins                                             *
 *                                                                           *
 * This file is part of the mtrx project, which is distributed under the MIT *
 * license. For the full terms, see the included LICENSE file.               *
 *****************************************************************************/

#include <stddef.h>

#include "sclr.h"


#ifndef MTRX_UTIL_H_
#define MTRX_UTIL_H_


/*****************************************************************************
 * Converts a string of ascii characters to a scalar value.                  *
 *                                                                           *
 * Fields:                                                                   *
 *     str - The string that contains the characters to convert.             *
 *     start - The index of the first character in the substring to be       *
 *             converted (inclusive).                                        *
 *     end - The index of the last character in the substring to be          *
 *           converted (exclusive).                                          *
 *                                                                           *
 * Returns: A pointer to the new vector.                                     *
 *****************************************************************************/
scalar_t atos(const char *str, size_t start, size_t end);

#endif

