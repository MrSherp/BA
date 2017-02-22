#ifndef __INCLUDEFROMCERES_H
#define __INCLUDEFROMCERES_H

#include <general.h>

#ifdef __GNUC__
#pragma GCC system_header
#endif




#ifdef USE_CERES
#include "ceres/ceres.h"
#include "glog/logging.h"
//#else
//     throw std::invalid_argument ( aol::strprintf ( "Ceres is not installed. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
#endif



#endif // __INCLUDEFROMCERES_H
