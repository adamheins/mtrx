#include <stddef.h>

#include "sclr.h"
#include "mtrx_util.h"

// Convert ascii to a scalar.
scalar_t atos(const char *str, size_t start, size_t end) {
  scalar_t value = 0;
  int sign = 1;
  size_t index = start;

  // Check for a negative sign.
  if (str[start] == '-') {
    sign = -1;
    ++index;
  }

  while (str[index] != '.' && index < end) {
    value = value * 10 + str[index] - '0';
    ++index;
  }

  scalar_t fraction = 0;
  size_t len = 0;

  // Get the value of the fractional part.
  while (++index < end) {
    fraction = fraction * 10 + str[index] - '0';
    ++len;
  }

  for (size_t i = 0; i < len; ++i)
    fraction /= 10;

  return sign * (value + fraction);
}
