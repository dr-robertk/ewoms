// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2011 by Andreas Lauser                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This file contains the parts of the local residual to
 *        calculate the heat flux in the fully coupled two-phase
 *        N-component model
 */
#ifndef DUMUX_MPNC_LOCAL_RESIDUAL_ENERGY_HH
#define DUMUX_MPNC_LOCAL_RESIDUAL_ENERGY_HH

#include <dumux/boxmodels/mpnc/mpncproperties.hh>

#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \brief Specialization of the energy module for the isothermal case.
 *
 * This class just does nothing.
 */
template <class TypeTag, bool enableEnergy/*=false*/>
class MPNCLocalResidualEnergy
{

    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { numPhases        = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents    = GET_PROP_VALUE(TypeTag, NumComponents) };
    typedef typename Dune::FieldVector<Scalar, numComponents>  ComponentVector;

public:
    static void computeStorage(EqVector &storage,
                               const VolumeVariables &volVars)
    {
        // do nothing, we're isothermal!
    }


    static void addPhaseStorage(EqVector &storage,
                                const VolumeVariables &volVars,
                                int phaseIdx)
    {
        // do nothing, we're isothermal!
    }

    static void phaseEnthalpyFlux(EqVector &flux,
                                  const EqVector &compMolFlux,
                                  const ElementContext &elemCtx,
                                  int scvfIdx,
                                  int timeIdx,
                                  int phaseIdx)
    {
        // do nothing, we're isothermal!
    }

    static void heatConduction(EqVector &result,
                               const ElementContext &elemCtx,
                               int scvfIdx,
                               int timeIdx)
    {
        // do nothing, we're isothermal!
    }


    static void computeFlux(EqVector & flux,
                            const ElementContext &elemCtx,
                            int scvfIdx,
                            int timeIdx,
                            const ComponentVector *molarFluxes)
    {
        // do nothing, we're isothermal!
    }

    static void computeSource(EqVector &result,
                              const ElementContext &elemCtx,
                              int scvIdx,
                              int timeIdx)
    {
        // do nothing, we're isothermal!
    }
};


template <class TypeTag>
class MPNCLocalResidualEnergy<TypeTag, /*enableEnergy=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim = GridView::dimension };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { energyEqIdx = Indices::energyEqIdx };

    typedef typename Dune::FieldVector<Scalar, numComponents> ComponentVector;

public:
    static void computeStorage(EqVector &storage,
                               const VolumeVariables &volVars)
    {
        storage[energyEqIdx] = 0;

        // energy of the fluids
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            addPhaseStorage(storage, volVars, phaseIdx);
        }

        // handle the heat capacity of the solid
        storage[energyEqIdx] +=
            volVars.fluidState().temperature(/*phaseIdx=*/0)
            * volVars.heatCapacitySolid()
            * (1.0 - volVars.porosity());
    }

    static void addPhaseStorage(EqVector &storage,
                                const VolumeVariables &volVars,
                                int phaseIdx)
    {
        const typename VolumeVariables::FluidState &fs =
            volVars.fluidState();

        // add the internal energy of the phase
        storage[energyEqIdx] +=
            volVars.porosity()
            * fs.density(phaseIdx)
            * fs.internalEnergy(phaseIdx)
            * fs.saturation(phaseIdx);
    }

    static void computeFlux(EqVector &flux,
                            const ElementContext &elemCtx,
                            int scvfIdx,
                            int timeIdx,
                            const ComponentVector molarFlux[numPhases])
    {
        flux[energyEqIdx] = 0.0;

        // fluid phases transport enthalpy individually
        for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx)
            computePhaseEnthalpyFlux(flux,
                                     elemCtx,
                                     scvfIdx,
                                     timeIdx,
                                     phaseIdx,
                                     molarFlux[phaseIdx]);

        //conduction is treated lumped in this model
        computeHeatConduction(flux,
                              elemCtx,
                              scvfIdx,
                              timeIdx);
    }

    static void computePhaseEnthalpyFlux(EqVector &flux,
                                         const ElementContext &elemCtx,
                                         int scvfIdx,
                                         int timeIdx,
                                         const int phaseIdx,
                                         const ComponentVector &molarPhaseFlux)
    {
        // use the phase enthalpy of the upstream vertex to calculate
        // the enthalpy transport
        const auto &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);

        int upIdx = fluxVars.upstreamIdx(phaseIdx);
        int dnIdx = fluxVars.downstreamIdx(phaseIdx);

        const VolumeVariables &up = elemCtx.volVars(upIdx, timeIdx);
        const VolumeVariables &dn = elemCtx.volVars(dnIdx, timeIdx);

        Scalar massFlux = 
            fluxVars.filterVelocityNormal(phaseIdx)
            * (up.fluidState().density(phaseIdx) * fluxVars.upstreamWeight(phaseIdx)
               + dn.fluidState().density(phaseIdx) * fluxVars.downstreamWeight(phaseIdx));

        flux[energyEqIdx] +=
            massFlux
            * (fluxVars.upstreamWeight(phaseIdx)*
               up.fluidState().enthalpy(phaseIdx)
               +
               fluxVars.downstreamWeight(phaseIdx)*
               dn.fluidState().enthalpy(phaseIdx));
    }

    static void computeHeatConduction(EqVector &flux,
                                      const ElementContext &elemCtx,
                                      int scvfIdx,
                                      int timeIdx)
    {
        const FluxVariables &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);

        // diffusive heat flux
        flux[energyEqIdx] +=
            - fluxVars.energyVars().temperatureGradNormal()
            * fluxVars.energyVars().heatConductivity();
    }



    static void computeSource(EqVector &source,
                              const ElementContext &elemCtx,
                              int scvIdx,
                              int timeIdx)
    {
        source[energyEqIdx] = 0.0;
    }
};
} // namespace Dumux

#endif // DUMUX_MPNC_LOCAL_RESIDUAL_ENERGY_HH
