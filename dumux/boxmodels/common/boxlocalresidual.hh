// $Id$
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief Caculates the residual of models based on the box scheme element-wise.
 */
#ifndef DUMUX_BOX_LOCAL_RESIDUAL_HH
#define DUMUX_BOX_LOCAL_RESIDUAL_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/valgrind.hh>
#include <dune/grid/common/genericreferenceelements.hh>

#include <boost/format.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/matrix.hh>

#include "boxproperties.hh"
#include "boxfvelementgeometry.hh"
#include "boxelementboundarytypes.hh"

namespace Dumux
{
/*!
 * \ingroup BoxModel
 * \brief Element-wise caculation of the residual matrix for models
 *        based on the box scheme .
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class BoxLocalResidual
{
private:
    typedef BoxLocalResidual<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalResidual)) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GridView::Grid::ctype CoordScalar;

    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container ReferenceElements;
    typedef typename RefElemProp::ReferenceElement ReferenceElement;

    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Element::Geometry Geometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSolutionVector)) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;

    typedef Dune::FieldMatrix<Scalar, numEq, numEq>  MatrixBlock;
    typedef Dune::Matrix<MatrixBlock> LocalBlockMatrix;

    // copying the local residual class is not a good idea
    BoxLocalResidual(const BoxLocalResidual &);

public:
    BoxLocalResidual()
    { }

    ~BoxLocalResidual()
    { }

    void init(Problem &prob)
    { problemPtr_ = &prob; }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero.
     */
    void eval(const Element &element)
    {
        FVElementGeometry fvGeom;
        fvGeom.update(gridView_(), element);
        ElementVolumeVariables volVarsPrev, volVarsCur;
        volVarsPrev.update(problem_(),
                           element,
                           fvGeom,
                           true /* oldSol? */);
        volVarsCur.update(problem_(),
                          element,
                          fvGeom,
                          false /* oldSol? */);
        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem_(), element, fvGeom);

        // this is pretty much a HACK because the internal state of
        // the problem is not supposed to be changed during the
        // evaluation of the residual. (Reasons: It is a violation of
        // abstraction, makes everything more prone to errors and is
        // not thread save.) The real solution are context objects!
        problem_().updateCouplingParams(element);

        asImp_().eval(element, fvGeom, volVarsPrev, volVarsCur, bcTypes);
    }

    /*!
     * \brief Compute the storage term for the current solution.
     *
     * This can be used to figure out how much of each conservation
     * quantity is inside the element.
     */
    void evalStorage(const Element &element)
    {
        FVElementGeometry fvGeom;
        fvGeom.update(gridView_(), element);
        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem_(), element, fvGeom);
        ElementVolumeVariables volVars;
        volVars.update(problem_(), element, fvGeom, false);
        
        residual_.resize(fvGeom.numVertices);
        residual_ = 0;

        elemPtr_ = &element;
        fvElemGeomPtr_ = &fvGeom;
        bcTypesPtr_ = &bcTypes;
        prevVolVarsPtr_ = 0;
        curVolVarsPtr_ = &volVars;
        asImp_().evalStorage_();
    }

    /*!
     * \brief Compute the flux term for the current solution.
     */
    void evalFluxes(const Element &element,
                    const ElementVolumeVariables &curVolVars)
    {
        FVElementGeometry fvGeom;
        fvGeom.update(gridView_(), element);
        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem_(), element, fvGeom);

        residual_.resize(fvGeom.numVertices);
        residual_ = 0;

        elemPtr_ = &element;
        fvElemGeomPtr_ = &fvGeom;
        bcTypesPtr_ = &bcTypes;
        prevVolVarsPtr_ = 0;
        curVolVarsPtr_ = &curVolVars;
        asImp_().evalFluxes_();
    }


    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero.
     */
    void eval(const Element &element,
              const FVElementGeometry &fvGeom,
              const ElementVolumeVariables &prevVolVars,
              const ElementVolumeVariables &curVolVars,
              const ElementBoundaryTypes &bcTypes)
    {
        //Valgrind::CheckDefined(fvGeom);
        Valgrind::CheckDefined(prevVolVars);
        Valgrind::CheckDefined(curVolVars);

#if HAVE_VALGRIND
        for (int i=0; i < fvGeom.numVertices; i++) {
            Valgrind::CheckDefined(prevVolVars[i]);
            Valgrind::CheckDefined(curVolVars[i]);
        }
#endif // HAVE_VALGRIND

        elemPtr_ = &element;
        fvElemGeomPtr_ = &fvGeom;
        bcTypesPtr_ = &bcTypes;
        prevVolVarsPtr_ = &prevVolVars;
        curVolVarsPtr_ = &curVolVars;

        // reset residual
        int numVerts = fvElemGeom_().numVertices;
        residual_.resize(numVerts);
        residual_ = 0;

        asImp_().evalFluxes_();
        asImp_().evalVolumeTerms_();

        // evaluate the boundary
        asImp_().evalBoundary_();

#if HAVE_VALGRIND
        for (int i=0; i < fvElemGeom_().numVertices; i++)
            Valgrind::CheckDefined(residual_[i]);
#endif // HAVE_VALGRIND
    }

    /*!
     * \brief Returns the local residual for a given sub-control
     *        volume of the element.
     */
    const ElementSolutionVector &residual() const
    { return residual_; }

    /*!
     * \brief Returns the local residual for a given sub-control
     *        volume of the element.
     */
    const PrimaryVariables residual(int scvIdx) const
    { return residual_[scvIdx]; }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    void evalBoundary_()
    {
        Dune::GeometryType geoType = elem_().geometry().type();
        const ReferenceElement &refElem = ReferenceElements::general(geoType);

        // evaluate boundary conditions for all intersections of
        // the current element
        IntersectionIterator isIt = gridView_().ibegin(elem_());
        const IntersectionIterator &endIt = gridView_().iend(elem_());
        for (; isIt != endIt; ++isIt)
        {
            // handle only faces on the boundary
            if (!isIt->boundary())
                continue;

            // Assemble the boundary for all vertices of the current
            // face
            int faceIdx = isIt->indexInInside();
            int numFaceVerts = refElem.size(faceIdx, 1, dim);
            for (int faceVertIdx = 0;
                 faceVertIdx < numFaceVerts;
                 ++faceVertIdx)
            {
                int elemVertIdx = refElem.subEntity(faceIdx,
                                                    1,
                                                    faceVertIdx,
                                                    dim);

                int boundaryFaceIdx =
                    fvElemGeom_().boundaryFaceIndex(faceIdx, faceVertIdx);

                // add the neuman fluxes of a single boundary segment
                evalBoundarySegment_(isIt,
                                     elemVertIdx,
                                     boundaryFaceIdx);
            }
        }
    }

    // handle boundary conditions for a single sub-control volume face
    void evalBoundarySegment_(const IntersectionIterator &isIt,
                              int scvIdx,
                              int boundaryFaceIdx)
    {
        // temporary vector to store the neumann boundary fluxes
        PrimaryVariables values(0.0);
        const BoundaryTypes &bcTypes = bcTypes_(scvIdx);

        // deal with neumann boundaries
        if (bcTypes.hasNeumann()) {
            problem_().neumann(values,
                               elem_(),
                               fvElemGeom_(),
                               *isIt,
                               scvIdx,
                               boundaryFaceIdx);
            values *= fvElemGeom_().boundaryFace[boundaryFaceIdx].area;
            Valgrind::CheckDefined(values);
            residual_[scvIdx] += values;
        }

        // deal with dirichlet boundaries
        if (bcTypes.hasDirichlet()) {
            problem_().dirichlet(values,
                                 elem_(),
                                 fvElemGeom_(),
                                 *isIt,
                                 scvIdx,
                                 boundaryFaceIdx);

            Valgrind::CheckDefined(values);
            for (int i = 0; i < numEq; ++i) {
                if (!bcTypes.isDirichlet(i))
                    continue;

                int pvIdx = bcTypes.eqToDirichletIndex(i);
                residual_[scvIdx][i] =
                    curPrimaryVars_(scvIdx)[pvIdx] - values[pvIdx];
            }
        }
    }

    void evalFluxes_()
    {
        // calculate the mass flux over the faces and subtract
        // it from the local rates
        for (int k = 0; k < fvElemGeom_().numEdges; k++)
        {
            int i = fvElemGeom_().subContVolFace[k].i;
            int j = fvElemGeom_().subContVolFace[k].j;

            PrimaryVariables flux;
            Valgrind::SetUndefined(flux);
            this->asImp_().computeFlux(flux, k);
            Valgrind::CheckDefined(flux);

            // subtract fluxes from the local mass rates of the
            // respective sub control volume adjacent to the face.
            for (int eq = 0; eq < numEq; ++ eq) {
                residual_[i][eq] -= flux[eq];
                residual_[j][eq] += flux[eq];
            }
        }
    }

    void evalStorage_()
    {
        // calculate the amount of conservation each quantity inside
        // all sub control volumes
        for (int i=0; i < fvElemGeom_().numVertices; i++) {
            residual_[i] = 0.0;
            this->asImp_().computeStorage(residual_[i], i, false /*isOldSol*/);
            residual_[i] *= fvElemGeom_().subContVol[i].volume;
        }
    }

    void evalVolumeTerms_()
    {
        // evaluate the volume terms (storage + source terms)
        for (int i=0; i < fvElemGeom_().numVertices; i++)
        {
            PrimaryVariables massContrib(0), tmp(0);

            // mass balance within the element. this is the
            // $\frac{m}{\partial t}$ term if using implicit
            // euler as time discretization.
            //
            // TODO (?): we might need a more explicit way for
            // doing the time discretization...
            this->asImp_().computeStorage(massContrib, i, false);
            this->asImp_().computeStorage(tmp, i, true);

            massContrib -= tmp;
            massContrib *=
                fvElemGeom_().subContVol[i].volume
                /
                problem_().timeManager().timeStepSize();
            
            for (int j = 0; j < numEq; ++j)
                residual_[i][j] += massContrib[j];

            // subtract the source term from the local rate
            PrimaryVariables source;
            this->asImp_().computeSource(source, i);
            source *= fvElemGeom_().subContVol[i].volume;

            for (int j = 0; j < numEq; ++j) {
                residual_[i][j] -= source[j];

                // make sure that only defined quantities where used
                // to calculate the residual.
                Valgrind::CheckDefined(residual_[i][j]);
            }
        }
    }

    /*!
     * \brief Returns a reference to the problem.
     */
    const Problem &problem_() const
    { return *problemPtr_; };

    /*!
     * \brief Returns a reference to the model.
     */
    const Model &model_() const
    { return problem_().model(); };

    /*!
     * \brief Returns a reference to the vertex mapper.
     */
    const VertexMapper &vertexMapper_() const
    { return problem_().vertexMapper(); };
    const GridView &gridView_() const
    { return problem_().gridView(); }

    const Element &elem_() const
    {
        Valgrind::CheckDefined(elemPtr_);
        return *elemPtr_;
    }

    const FVElementGeometry &fvElemGeom_() const
    {
        Valgrind::CheckDefined(fvElemGeomPtr_);
        return *fvElemGeomPtr_;
    }

    const PrimaryVariables &curPrimaryVars_(int i) const
    {
        return curVolVars_(i).primaryVars();
    }

    const ElementVolumeVariables &curVolVars_() const
    {
        Valgrind::CheckDefined(curVolVarsPtr_);
        return *curVolVarsPtr_;
    }
    const VolumeVariables &curVolVars_(int i) const
    {
        return curVolVars_()[i];
    }

    const ElementVolumeVariables &prevVolVars_() const
    {
        Valgrind::CheckDefined(prevVolVarsPtr_);
        return *prevVolVarsPtr_;
    }
    const VolumeVariables &prevVolVars_(int i) const
    {
        return prevVolVars_()[i];
    }

    const ElementBoundaryTypes &bcTypes_() const
    {
        Valgrind::CheckDefined(bcTypesPtr_);
        return *bcTypesPtr_;
    }
    const BoundaryTypes &bcTypes_(int i) const
    {
        return bcTypes_()[i];
    }

protected:
    ElementSolutionVector residual_;

    // The problem we would like to solve
    Problem *problemPtr_;

    const Element *elemPtr_;
    const FVElementGeometry *fvElemGeomPtr_;

    // current and previous secondary variables for the element
    const ElementVolumeVariables *prevVolVarsPtr_;
    const ElementVolumeVariables *curVolVarsPtr_;

    const ElementBoundaryTypes *bcTypesPtr_;
};
}

#endif
