/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <stdlib.h> // for malloc
#include <string.h> // for memset
#include <string>
using std::cout;
using std::endl;
using std::ptrdiff_t;
using std::size_t;

#include "FciCompat.h"
using boost::format;

#include "CxMemoryStack.h"
#include "CxOpenMpProxy.h"
#include "CxPodArray.h"
using namespace ct;

extern bool g_BosonicSigns;
// clang-format off
#include "CxTypes.h"
#include "CxAlgebra.h"
#include "CxStorageDevice.h"

#include "FciDump.h"
#include "FciVec.h"
#include "FciSubspace.h"
#include "FciOps.h"
#include "FciCsf.h"
// clang-format on
