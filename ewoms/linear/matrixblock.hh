// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#ifndef EWOMS_MATRIX_BLOCK_HH
#define EWOMS_MATRIX_BLOCK_HH

#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/amg.hh>

#include <dune/common/fmatrix.hh>

namespace Ewoms {
namespace MatrixBlockHelp {

template <typename K, int m, int n>
static inline void invertMatrix(Dune::FieldMatrix<K, m, n>& matrix)
{ matrix.invert(); }

template <typename K>
static inline void invertMatrix(Dune::FieldMatrix<K, 1, 1>& matrix)
{
    matrix[0][0] = 1.0/matrix[0][0];
}

template <typename K>
static inline void invertMatrix(Dune::FieldMatrix<K, 2, 2>& matrix)
{
    Dune::FieldMatrix<K, 2, 2> tmp(matrix);
    Dune::FMatrixHelp::invertMatrix(tmp, matrix);
}

template <typename K>
static inline void invertMatrix(Dune::FieldMatrix<K, 3, 3>& matrix)
{
    Dune::FieldMatrix<K, 3, 3> tmp(matrix);
    Dune::FMatrixHelp::invertMatrix(tmp, matrix);
}

template <typename K>
static inline void invertMatrix(Dune::FieldMatrix<K, 4, 4>& matrix)
{
    Dune::FieldMatrix<K, 4, 4> tmp(matrix);

    matrix[0][0] =
        tmp[1][1]*tmp[2][2]*tmp[3][3] -
        tmp[1][1]*tmp[2][3]*tmp[3][2] -
        tmp[2][1]*tmp[1][2]*tmp[3][3] +
        tmp[2][1]*tmp[1][3]*tmp[3][2] +
        tmp[3][1]*tmp[1][2]*tmp[2][3] -
        tmp[3][1]*tmp[1][3]*tmp[2][2];

    matrix[1][0] =
       -tmp[1][0]*tmp[2][2]*tmp[3][3] +
        tmp[1][0]*tmp[2][3]*tmp[3][2] +
        tmp[2][0]*tmp[1][2]*tmp[3][3] -
        tmp[2][0]*tmp[1][3]*tmp[3][2] -
        tmp[3][0]*tmp[1][2]*tmp[2][3] +
        tmp[3][0]*tmp[1][3]*tmp[2][2];

    matrix[2][0] =
        tmp[1][0]*tmp[2][1]*tmp[3][3] -
        tmp[1][0]*tmp[2][3]*tmp[3][1] -
        tmp[2][0]*tmp[1][1]*tmp[3][3] +
        tmp[2][0]*tmp[1][3]*tmp[3][1] +
        tmp[3][0]*tmp[1][1]*tmp[2][3] -
        tmp[3][0]*tmp[1][3]*tmp[2][1];

    matrix[3][0] =
       -tmp[1][0]*tmp[2][1]*tmp[3][2] +
        tmp[1][0]*tmp[2][2]*tmp[3][1] +
        tmp[2][0]*tmp[1][1]*tmp[3][2] -
        tmp[2][0]*tmp[1][2]*tmp[3][1] -
        tmp[3][0]*tmp[1][1]*tmp[2][2] +
        tmp[3][0]*tmp[1][2]*tmp[2][1];

    matrix[0][1] =
       -tmp[0][1]*tmp[2][2]*tmp[3][3] +
        tmp[0][1]*tmp[2][3]*tmp[3][2] +
        tmp[2][1]*tmp[0][2]*tmp[3][3] -
        tmp[2][1]*tmp[0][3]*tmp[3][2] -
        tmp[3][1]*tmp[0][2]*tmp[2][3] +
        tmp[3][1]*tmp[0][3]*tmp[2][2];

    matrix[1][1] =
        tmp[0][0]*tmp[2][2]*tmp[3][3] -
        tmp[0][0]*tmp[2][3]*tmp[3][2] -
        tmp[2][0]*tmp[0][2]*tmp[3][3] +
        tmp[2][0]*tmp[0][3]*tmp[3][2] +
        tmp[3][0]*tmp[0][2]*tmp[2][3] -
        tmp[3][0]*tmp[0][3]*tmp[2][2];

    matrix[2][1] =
       -tmp[0][0]*tmp[2][1]*tmp[3][3] +
        tmp[0][0]*tmp[2][3]*tmp[3][1] +
        tmp[2][0]*tmp[0][1]*tmp[3][3] -
        tmp[2][0]*tmp[0][3]*tmp[3][1] -
        tmp[3][0]*tmp[0][1]*tmp[2][3] +
        tmp[3][0]*tmp[0][3]*tmp[2][1];

    matrix[3][1] =
        tmp[0][0]*tmp[2][1]*tmp[3][2] -
        tmp[0][0]*tmp[2][2]*tmp[3][1] -
        tmp[2][0]*tmp[0][1]*tmp[3][2] +
        tmp[2][0]*tmp[0][2]*tmp[3][1] +
        tmp[3][0]*tmp[0][1]*tmp[2][2] -
        tmp[3][0]*tmp[0][2]*tmp[2][1];

    matrix[0][2] =
        tmp[0][1]*tmp[1][2]*tmp[3][3] -
        tmp[0][1]*tmp[1][3]*tmp[3][2] -
        tmp[1][1]*tmp[0][2]*tmp[3][3] +
        tmp[1][1]*tmp[0][3]*tmp[3][2] +
        tmp[3][1]*tmp[0][2]*tmp[1][3] -
        tmp[3][1]*tmp[0][3]*tmp[1][2];

    matrix[1][2] =
       -tmp[0][0] *tmp[1][2]*tmp[3][3] +
        tmp[0][0]*tmp[1][3]*tmp[3][2] +
        tmp[1][0]*tmp[0][2]*tmp[3][3] -
        tmp[1][0]*tmp[0][3]*tmp[3][2] -
        tmp[3][0]*tmp[0][2]*tmp[1][3] +
        tmp[3][0]*tmp[0][3]*tmp[1][2];

    matrix[2][2] =
        tmp[0][0]*tmp[1][1]*tmp[3][3] -
        tmp[0][0]*tmp[1][3]*tmp[3][1] -
        tmp[1][0]*tmp[0][1]*tmp[3][3] +
        tmp[1][0]*tmp[0][3]*tmp[3][1] +
        tmp[3][0]*tmp[0][1]*tmp[1][3] -
        tmp[3][0]*tmp[0][3]*tmp[1][1];

    matrix[3][2] =
       -tmp[0][0]*tmp[1][1]*tmp[3][2] +
        tmp[0][0]*tmp[1][2]*tmp[3][1] +
        tmp[1][0]*tmp[0][1]*tmp[3][2] -
        tmp[1][0]*tmp[0][2]*tmp[3][1] -
        tmp[3][0]*tmp[0][1]*tmp[1][2] +
        tmp[3][0]*tmp[0][2]*tmp[1][1];

    matrix[0][3] =
       -tmp[0][1]*tmp[1][2]*tmp[2][3] +
        tmp[0][1]*tmp[1][3]*tmp[2][2] +
        tmp[1][1]*tmp[0][2]*tmp[2][3] -
        tmp[1][1]*tmp[0][3]*tmp[2][2] -
        tmp[2][1]*tmp[0][2]*tmp[1][3] +
        tmp[2][1]*tmp[0][3]*tmp[1][2];

    matrix[1][3] =
        tmp[0][0]*tmp[1][2]*tmp[2][3] -
        tmp[0][0]*tmp[1][3]*tmp[2][2] -
        tmp[1][0]*tmp[0][2]*tmp[2][3] +
        tmp[1][0]*tmp[0][3]*tmp[2][2] +
        tmp[2][0]*tmp[0][2]*tmp[1][3] -
        tmp[2][0]*tmp[0][3]*tmp[1][2];

    matrix[2][3] =
       -tmp[0][0]*tmp[1][1]*tmp[2][3] +
        tmp[0][0]*tmp[1][3]*tmp[2][1] +
        tmp[1][0]*tmp[0][1]*tmp[2][3] -
        tmp[1][0]*tmp[0][3]*tmp[2][1] -
        tmp[2][0]*tmp[0][1]*tmp[1][3] +
        tmp[2][0]*tmp[0][3]*tmp[1][1];

    matrix[3][3] =
        tmp[0][0]*tmp[1][1]*tmp[2][2] -
        tmp[0][0]*tmp[1][2]*tmp[2][1] -
        tmp[1][0]*tmp[0][1]*tmp[2][2] +
        tmp[1][0]*tmp[0][2]*tmp[2][1] +
        tmp[2][0]*tmp[0][1]*tmp[1][2] -
        tmp[2][0]*tmp[0][2]*tmp[1][1];

    K det =
        tmp[0][0]*matrix[0][0] +
        tmp[0][1]*matrix[1][0] +
        tmp[0][2]*matrix[2][0] +
        tmp[0][3]*matrix[3][0];

    // return identity for singular or nearly singular matrices.
    if (std::abs(det) < 1e-40)
        matrix = std::numeric_limits<K>::quiet_NaN();
    else
        matrix *= 1.0/det;
}

template <typename K, int m, int n, class X, class Y>
static void mv( const Dune::FieldMatrix<K, n, m>& matrix, const X& x, Y& y )
{
    matrix.mv( x, y );
}

template <typename K, int m, int n, class X, class Y>
static void umv( const Dune::FieldMatrix<K, n, m>& matrix, const X& x, Y& y )
{
    matrix.umv( x, y );
}

template <typename K, int m, int n, class X, class Y>
static void mmv( const Dune::FieldMatrix<K, n, m>& matrix, const X& x, Y& y )
{
    matrix.mmv( x, y );
}

template <typename K, class X, class Y>
static void mv( const Dune::FieldMatrix<K, 3, 3>& matrix, const X& x, Y& y )
{
  y[ 0 ] = matrix[ 0 ][ 0 ] * x[ 0 ] + matrix[ 0 ][ 1 ] * x[ 1 ] + matrix[ 0 ][ 2 ] * x[ 2 ] ;
  y[ 1 ] = matrix[ 1 ][ 0 ] * x[ 0 ] + matrix[ 1 ][ 1 ] * x[ 1 ] + matrix[ 1 ][ 2 ] * x[ 2 ] ;
  y[ 2 ] = matrix[ 2 ][ 0 ] * x[ 0 ] + matrix[ 2 ][ 1 ] * x[ 1 ] + matrix[ 2 ][ 2 ] * x[ 2 ] ;
}

template <typename K, class X, class Y>
static void umv( const Dune::FieldMatrix<K, 3, 3>& matrix, const X& x, Y& y )
{
  y[ 0 ] += matrix[ 0 ][ 0 ] * x[ 0 ] + matrix[ 0 ][ 1 ] * x[ 1 ] + matrix[ 0 ][ 2 ] * x[ 2 ] ;
  y[ 1 ] += matrix[ 1 ][ 0 ] * x[ 0 ] + matrix[ 1 ][ 1 ] * x[ 1 ] + matrix[ 1 ][ 2 ] * x[ 2 ] ;
  y[ 2 ] += matrix[ 2 ][ 0 ] * x[ 0 ] + matrix[ 2 ][ 1 ] * x[ 1 ] + matrix[ 2 ][ 2 ] * x[ 2 ] ;
}

template <typename K, class X, class Y>
static void mmv( const Dune::FieldMatrix<K, 3, 3>& matrix, const X& x, Y& y )
{
  y[ 0 ] -= matrix[ 0 ][ 0 ] * x[ 0 ] + matrix[ 0 ][ 1 ] * x[ 1 ] + matrix[ 0 ][ 2 ] * x[ 2 ] ;
  y[ 1 ] -= matrix[ 1 ][ 0 ] * x[ 0 ] + matrix[ 1 ][ 1 ] * x[ 1 ] + matrix[ 1 ][ 2 ] * x[ 2 ] ;
  y[ 2 ] -= matrix[ 2 ][ 0 ] * x[ 0 ] + matrix[ 2 ][ 1 ] * x[ 1 ] + matrix[ 2 ][ 2 ] * x[ 2 ] ;
}

template<typename K, int m, int n, typename M2>
static Dune::FieldMatrix<K, n, m>&
leftmultiply (Dune::FieldMatrix<K, n, m>& A, const Dune::DenseMatrix<M2>& B)
{
  return A.leftmultiply( B );
}

template <class K>
static inline K dot( const Dune::FieldVector< K, 3 >& a, const Dune::FieldVector< K, 3 >& b)
{
    return a[ 0 ] * b[ 0 ] + a[ 1 ] * b[ 1 ] + a[ 2 ] * b[ 2 ];
}

template <class K>
static inline K dot( const Dune::FieldVector< K, 3 >& a, const Dune::FieldMatrix< K, 3, 3 >& b, const int col)
{
    return a[ 0 ] * b[ 0 ][ col ] + a[ 1 ] * b[ 1 ][ col ] + a[ 2 ] * b[ 2 ][ col ];
}

template <typename K, class X, class Y>
static inline void
leftmultiply (Dune::FieldMatrix<K, 3, 3>& A, const Dune::FieldMatrix< K, 3, 3>& B, const Dune::FieldMatrix< K, 3, 3>& C)
{
    A[ 0 ][ 0 ] = B[ 0 ][ 0 ] * C[ 0 ][ 0 ] + B[ 0 ][ 1 ] * C[ 0 ][ 1 ] + B[ 0 ][ 2 ] * C[ 0 ][ 2 ];
    A[ 0 ][ 1 ] = B[ 0 ][ 0 ] * C[ 1 ][ 0 ] + B[ 0 ][ 1 ] * C[ 1 ][ 1 ] + B[ 0 ][ 2 ] * C[ 1 ][ 2 ];
    A[ 0 ][ 2 ] = B[ 0 ][ 0 ] * C[ 2 ][ 0 ] + B[ 0 ][ 1 ] * C[ 2 ][ 1 ] + B[ 0 ][ 2 ] * C[ 2 ][ 2 ];
    A[ 1 ][ 0 ] = B[ 1 ][ 0 ] * C[ 0 ][ 0 ] + B[ 1 ][ 1 ] * C[ 0 ][ 1 ] + B[ 1 ][ 2 ] * C[ 0 ][ 2 ];
    A[ 1 ][ 1 ] = B[ 1 ][ 0 ] * C[ 1 ][ 0 ] + B[ 1 ][ 1 ] * C[ 1 ][ 1 ] + B[ 1 ][ 2 ] * C[ 1 ][ 2 ];
    A[ 1 ][ 2 ] = B[ 1 ][ 0 ] * C[ 2 ][ 0 ] + B[ 1 ][ 1 ] * C[ 2 ][ 1 ] + B[ 1 ][ 2 ] * C[ 2 ][ 2 ];
    A[ 2 ][ 0 ] = B[ 2 ][ 0 ] * C[ 0 ][ 0 ] + B[ 2 ][ 1 ] * C[ 0 ][ 1 ] + B[ 2 ][ 2 ] * C[ 0 ][ 2 ];
    A[ 2 ][ 1 ] = B[ 2 ][ 0 ] * C[ 1 ][ 0 ] + B[ 2 ][ 1 ] * C[ 1 ][ 1 ] + B[ 2 ][ 2 ] * C[ 1 ][ 2 ];
    A[ 2 ][ 2 ] = B[ 2 ][ 0 ] * C[ 2 ][ 0 ] + B[ 2 ][ 1 ] * C[ 2 ][ 1 ] + B[ 2 ][ 2 ] * C[ 2 ][ 2 ];
}

template <typename K, class X, class Y>
static inline void
rightmultiply (Dune::FieldMatrix<K, 3, 3>& A, const Dune::FieldMatrix< K, 3, 3>& B, const Dune::FieldMatrix< K, 3, 3>& C)
{
    A[ 0 ][ 0 ] = B[ 0 ][ 0 ] * C[ 0 ][ 0 ] + B[ 0 ][ 1 ] * C[ 1 ][ 0 ] + B[ 0 ][ 2 ] * C[ 2 ][ 0 ];
    A[ 0 ][ 1 ] = B[ 0 ][ 0 ] * C[ 0 ][ 1 ] + B[ 0 ][ 1 ] * C[ 1 ][ 1 ] + B[ 0 ][ 2 ] * C[ 2 ][ 1 ];
    A[ 0 ][ 2 ] = B[ 0 ][ 0 ] * C[ 0 ][ 2 ] + B[ 0 ][ 1 ] * C[ 1 ][ 2 ] + B[ 0 ][ 2 ] * C[ 2 ][ 2 ];
    A[ 1 ][ 0 ] = B[ 1 ][ 0 ] * C[ 0 ][ 0 ] + B[ 1 ][ 1 ] * C[ 1 ][ 0 ] + B[ 1 ][ 2 ] * C[ 2 ][ 0 ];
    A[ 1 ][ 1 ] = B[ 1 ][ 0 ] * C[ 0 ][ 1 ] + B[ 1 ][ 1 ] * C[ 1 ][ 1 ] + B[ 1 ][ 2 ] * C[ 2 ][ 1 ];
    A[ 1 ][ 2 ] = B[ 1 ][ 0 ] * C[ 0 ][ 2 ] + B[ 1 ][ 1 ] * C[ 1 ][ 2 ] + B[ 1 ][ 2 ] * C[ 2 ][ 2 ];
    A[ 2 ][ 0 ] = B[ 2 ][ 0 ] * C[ 0 ][ 0 ] + B[ 2 ][ 1 ] * C[ 1 ][ 0 ] + B[ 2 ][ 2 ] * C[ 2 ][ 0 ];
    A[ 2 ][ 1 ] = B[ 2 ][ 0 ] * C[ 0 ][ 1 ] + B[ 2 ][ 1 ] * C[ 1 ][ 1 ] + B[ 2 ][ 2 ] * C[ 2 ][ 1 ];
    A[ 2 ][ 2 ] = B[ 2 ][ 0 ] * C[ 0 ][ 2 ] + B[ 2 ][ 1 ] * C[ 1 ][ 2 ] + B[ 2 ][ 2 ] * C[ 2 ][ 2 ];
}

template <typename K, class X, class Y>
static Dune::FieldMatrix<K, 3, 3>&
leftmultiply (Dune::FieldMatrix<K, 3, 3>& A, const Dune::FieldMatrix< K, 3, 3>& B)
{
    Dune::FieldMatrix< K, 3, 3 > tmp;

    // transpose matrix
    tmp[ 0 ][ 0 ] = A[ 0 ][ 0 ];
    tmp[ 0 ][ 1 ] = A[ 1 ][ 0 ];
    tmp[ 0 ][ 2 ] = A[ 2 ][ 0 ];
    tmp[ 1 ][ 0 ] = A[ 0 ][ 1 ];
    tmp[ 1 ][ 1 ] = A[ 1 ][ 1 ];
    tmp[ 1 ][ 2 ] = A[ 2 ][ 1 ];
    tmp[ 2 ][ 0 ] = A[ 0 ][ 2 ];
    tmp[ 2 ][ 1 ] = A[ 1 ][ 2 ];
    tmp[ 2 ][ 2 ] = A[ 2 ][ 2 ];

    leftMultiply( A, B, tmp );
    return A;
}

template<typename K, int m, int n, typename M2>
static Dune::FieldMatrix<K, n, m>&
rightmultiply (Dune::FieldMatrix<K, n, m>& A, const Dune::DenseMatrix<M2>& B)
{
  return A.rightmultiply( B );
}

template <typename K, class X, class Y>
static Dune::FieldMatrix<K, 3, 3>&
rightmultiply (Dune::FieldMatrix<K, 3, 3>& A, const Dune::FieldMatrix< K, 3, 3>& B)
{
    Dune::FieldMatrix< K, 3, 3 > tmp( A );
    rightmultiply( A, tmp, B );
    return A;
}
} // namespace MatrixBlockHelp

template <class Scalar, int n, int m>
class MatrixBlock : public Dune::FieldMatrix<Scalar, n, m>
{
    typedef MatrixBlock< Scalar, n, m> ThisType;
public:
    typedef Dune::FieldMatrix<Scalar, n, m>  BaseType;

    using BaseType::operator= ;
    using BaseType::rows;
    using BaseType::cols;

    MatrixBlock()
        : BaseType(Scalar(0.0))
    {}

    explicit MatrixBlock(const Scalar value)
        : BaseType(value)
    {}

    void invert()
    { Ewoms::MatrixBlockHelp::invertMatrix(asBase()); }

    template<class X, class Y>
    void mv (const X& x, Y& y) const { Ewoms::MatrixBlockHelp::mv( asBase(), x, y ); }

    template<class X, class Y>
    void umv (const X& x, Y& y) const { Ewoms::MatrixBlockHelp::umv( asBase(), x, y ); }

    template<class X, class Y>
    void mmv (const X& x, Y& y) const { Ewoms::MatrixBlockHelp::mmv( asBase(), x, y ); }

    //! Multiplies M from the left to this matrix
    template<typename M2>
    ThisType& leftmultiply (const Dune::DenseMatrix<M2>& M)
    {
      return static_cast< ThisType& > (Ewoms::MatrixBlockHelp::leftmultiply( asBase(), static_cast< const BaseType& > (M) ));
    }

    //! Multiplies M from the right to this matrix
    template<typename M2>
    ThisType& rightmultiply (const Dune::DenseMatrix<M2>& M)
    {
      return static_cast< ThisType& > (Ewoms::MatrixBlockHelp::rightmultiply( asBase(), static_cast< const BaseType& > (M) ));
    }

    const BaseType& asBase() const
    { return static_cast<const BaseType&>(*this); }

    BaseType& asBase()
    { return static_cast<BaseType&>(*this); }
};



} // namespace Ewoms

namespace Dune {

template<class K, int n, int m>
void print_row(std::ostream& s, const Ewoms::MatrixBlock<K, n, m>& A,
               typename FieldMatrix<K, n, m>::size_type I,
               typename FieldMatrix<K, n, m>::size_type J,
               typename FieldMatrix<K, n, m>::size_type therow,
               int width,
               int precision)
{ print_row(s, A.asBase(), I, J, therow, width, precision); }

template<class K, int n, int m>
K& firstmatrixelement(Ewoms::MatrixBlock<K, n, m>& A)
{ return firstmatrixelement(A.asBase()); }

template <typename Scalar, int n, int m>
struct MatrixDimension<Ewoms::MatrixBlock<Scalar, n, m> >
    : public MatrixDimension<typename Ewoms::MatrixBlock<Scalar, n, m>::BaseType>
{ };


#if HAVE_UMFPACK
/// \brief UMFPack specialization for Ewoms::MatrixBlock to make AMG happy
///
/// Without this the empty default implementation would be used.
template <typename T, typename A, int n, int m>
class UMFPack<BCRSMatrix<Ewoms::MatrixBlock<T, n, m>, A> >
    : public UMFPack<BCRSMatrix<FieldMatrix<T, n, m>, A> >
{
    typedef UMFPack<BCRSMatrix<FieldMatrix<T, n, m>, A> > Base;
    typedef BCRSMatrix<FieldMatrix<T, n, m>, A> Matrix;

public:
    typedef BCRSMatrix<Ewoms::MatrixBlock<T, n, m>, A> RealMatrix;

    UMFPack(const RealMatrix& matrix, int verbose, bool)
        : Base(reinterpret_cast<const Matrix&>(matrix), verbose)
    {}
};
#endif

#if HAVE_SUPERLU
/// \brief SuperLU specialization for Ewoms::MatrixBlock to make AMG happy
///
/// Without this the empty default implementation would be used.
template <typename T, typename A, int n, int m>
class SuperLU<BCRSMatrix<Ewoms::MatrixBlock<T, n, m>, A> >
    : public SuperLU<BCRSMatrix<FieldMatrix<T, n, m>, A> >
{
    typedef SuperLU<BCRSMatrix<FieldMatrix<T, n, m>, A> > Base;
    typedef BCRSMatrix<FieldMatrix<T, n, m>, A> Matrix;

public:
    typedef BCRSMatrix<Ewoms::MatrixBlock<T, n, m>, A> RealMatrix;

    SuperLU(const RealMatrix& matrix, int verbose, bool reuse=true)
        : Base(reinterpret_cast<const Matrix&>(matrix), verbose, reuse)
    {}
};
#endif

} // end namespace Dune


#endif
