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
/*!
 * \file
 * \copydoc Ewoms::Linear::ISTLMatrixBackend
 */
#ifndef EWOMS_ISTL_MATRIX_BACKEND_HH
#define EWOMS_ISTL_MATRIX_BACKEND_HH

#include <dune/istl/bcrsmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

namespace Ewoms {
namespace Linear {

/*!
 * \ingroup Linear
 * \brief A backend for BCRSMatrix from dune-istl.
 */
template<class Block, class A=std::allocator< Block > >
class ISTLMatrixBackend
{
    typedef A AllocatorType;
public:
    // \brief Implementation of matrix
    typedef Dune::BCRSMatrix< Block, AllocatorType >   Matrix;

    // \brief block type forming the matrix entries (same as Block)
    typedef typename Matrix :: block_type          block_type;

    /*!
     * \brief Constructor creating an empty matrix.
     */
    ISTLMatrixBackend( const std::string& name, const size_t rows, const size_t columns )
        : rows_( rows )
        , columns_( columns )
        , matrix_()
    {}

    /*!
     * \brief Constructor creating an empty matrix.
     */
    explicit ISTLMatrixBackend( const std::string& name, const int rows, const int columns )
        : ISTLMatrixBackend( name, size_t(rows), size_t(columns) )
    {}

    /*!
     * \brief Constructor creating an empty matrix.
     */
    template <class DomainSpace, class RangeSpace>
    ISTLMatrixBackend( const std::string& name, const DomainSpace& domainSpace, const RangeSpace& rangeSpace )
        : ISTLMatrixBackend( name, domainSpace.size()/DomainSpace::dimRange, rangeSpace.size()/RangeSpace::dimRange )
    {
    }

    /*!
     * \brief Allocate matrix structure give a sparsity pattern.
     */
    template <class Set>
    inline void reserve( const std::vector< Set >& sparsityPattern )
    {
        // allocate raw matrix
        matrix_.reset( new Matrix(rows_, columns_, Matrix::random) );

        // make sure sparsityPattern is consistent with number of rows
        assert( rows_ == sparsityPattern.size() );

        // allocate space for the rows of the matrix
        for (size_t dofIdx = 0; dofIdx < rows_; ++ dofIdx)
        {
            matrix_->setrowsize(dofIdx, sparsityPattern[dofIdx].size());
        }

        matrix_->endrowsizes();

        // fill the rows with indices. each degree of freedom talks to
        // all of its neighbors. (it also talks to itself since
        // degrees of freedom are sometimes quite egocentric.)
        for (size_t dofIdx = 0; dofIdx < rows_; ++ dofIdx)
        {
            auto nIt    = sparsityPattern[dofIdx].begin();
            auto nEndIt = sparsityPattern[dofIdx].end();
            for (; nIt != nEndIt; ++nIt)
            {
                matrix_->addindex(dofIdx, *nIt);
            }
        }
        matrix_->endindices();
    }

    /*!
     * \brief Return constant reference to matrix implementation.
     */
    inline Matrix& matrix() { return *matrix_; }
    inline const Matrix& matrix() const { return *matrix_; }

    /*!
     * \brief Return number of rows of the matrix.
     */
    inline size_t rows () const { return rows_; }

    /*!
     * \brief Return number of columns of the matrix.
     */
    inline size_t cols () const { return columns_; }

    /*!
     * \brief Set all matrix entries to zero.
     */
    inline void clear()
    {
        (*matrix_) = typename block_type :: field_type(0);
    }

    /*!
     * \brief Set given row to zero except for the diagonal entry which is set to one.
     */
    inline void unitRow( const size_t row )
    {
        block_type idBlock( 0 );
        for (int i = 0; i < idBlock.rows; ++i)
            idBlock[i][i] = 1.0;

        auto& matRow = (*matrix_)[ row ];
        auto colIt = matRow.begin();
        const auto& colEndIt = matRow.end();
        for (; colIt != colEndIt; ++colIt)
        {
            if( colIt.index() == row )
                *colIt = idBlock;
            else
                *colIt = 0.0;
        }
    }

    /*!
     * \brief Fill given block with entries stored in the matrix.
     */
    inline void getBlock( const size_t row, const size_t col, block_type& entry ) const
    {
        entry = (*matrix_)[ row ][ col ];
    }

    /*!
     * \brief Set matrix block to given block.
     */
    inline void setBlock( const size_t row, const size_t col, const block_type& entry )
    {
        (*matrix_)[ row ][ col ] = entry;
    }

    /*!
     * \brief Add block to matrix block.
     */
    inline void addBlock( const size_t row, const size_t col, const block_type& entry )
    {
        (*matrix_)[ row ][ col ] += entry;
    }

    /*!
     * \brief Synchronize matrix and finalize building stage.
     */
    inline void communicate()
    {
        // nothing to do here
        // may call compress when implicit build mode is used
    }

protected:
    size_t rows_;
    size_t columns_;

    std::unique_ptr< Matrix > matrix_;
};

} // namespace Linear
} // namespace Ewoms

#endif
