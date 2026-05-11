/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#ifndef CXTYPES_H
#define CXTYPES_H

#define _for_each(it, con)                                                     \
  for ((it) = (con).begin(); (it) != (con).end(); ++(it))

#include "FciCompat.h"

using boost::int16_t;
using boost::int32_t;
using boost::int64_t;
using boost::int8_t;
using boost::uint16_t;
using boost::uint32_t;
using boost::uint64_t;
using boost::uint8_t;
using std::ptrdiff_t;
using std::size_t;

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned int uint;

#include "CxDefs.h"
#define RESTRICT AIC_RP

namespace ct {
struct FIntrusivePtrDest;
}

void intrusive_ptr_add_ref(ct::FIntrusivePtrDest const *pExpr);
void intrusive_ptr_release(ct::FIntrusivePtrDest const *pExpr);

namespace ct {
/// A base class for reference counted objects. Classes derived from this can
/// be used as target for boost::intrusive_ptr.
struct FIntrusivePtrDest {
  FIntrusivePtrDest() : m_RefCount(0){};
  inline virtual ~FIntrusivePtrDest() = 0;

  mutable int m_RefCount;
  friend void ::intrusive_ptr_add_ref(FIntrusivePtrDest const *Expr);
  friend void ::intrusive_ptr_release(FIntrusivePtrDest const *Expr);
};

inline FIntrusivePtrDest::~FIntrusivePtrDest(){};
} // namespace ct

inline void intrusive_ptr_add_ref(ct::FIntrusivePtrDest const *pExpr) {
  pExpr->m_RefCount += 1;
}

inline void intrusive_ptr_release(ct::FIntrusivePtrDest const *pExpr) {
  assert(pExpr->m_RefCount > 0);
  pExpr->m_RefCount -= 1;
  if (pExpr->m_RefCount == 0)
    delete pExpr;
}

#endif // CXTYPES_H

// kate: space-indent on; tab-indent on; backspace-indent on; tab-width 4;
// indent-width 4; mixedindent off; indent-mode normal;
