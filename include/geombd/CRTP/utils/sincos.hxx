/**
 *    \file include/geombd/CRTP/utils/sincos.hxx
 *    \author Alvaro Paz, Gustavo Arechavaleta
 *    \version 1.0
 *    \date 2021
 *
 *    Class to implement template-based sin cos functions
 *    Copyright (c) 2021 Cinvestav
 *    This library is distributed under the MIT License.
 */

#ifndef GEOMBD_SINCOS_HXX
#define GEOMBD_SINCOS_HXX

#pragma GCC optimize("O3")
//#pragma GCC flag("fno-builtin lm")

#include <cmath>

namespace geo{

  //! SINCOS
  //!------------------------------------------------------------------------------!//

  // Forward declaration
  template<typename ScalarT> struct SINCOSAlgo;

  template<typename ScalarT>
  inline static void SINCOS(const ScalarT & qi, ScalarT* sqi, ScalarT* cqi)
  { SINCOSAlgo<ScalarT>::run(qi, sqi, cqi); }

  template<typename ScalarT>
  struct SINCOSAlgo
  {
    inline static void run(const ScalarT & qi, ScalarT* sqi, ScalarT* cqi)
    {  (*sqi) = std::sin(qi); (*cqi) = std::cos(qi);  }
  };

  template<>
  struct SINCOSAlgo<float>
  {
    inline static void run(const float & qi, float* sqi, float* cqi)
    {  sincosf(qi, sqi, cqi);  }
  };

  template<>
  struct SINCOSAlgo<double>
  {
    inline static void run(const double & qi, double* sqi, double* cqi)
    {  sincos(qi, sqi, cqi);  }
  };

  template<>
  struct SINCOSAlgo<long double>
  {
    inline static void run(const long double & qi, long double* sqi, long double* cqi)
    {  sincosl(qi, sqi, cqi);  }
  };

}

#endif // GEOMBD_SINCOS_HXX
