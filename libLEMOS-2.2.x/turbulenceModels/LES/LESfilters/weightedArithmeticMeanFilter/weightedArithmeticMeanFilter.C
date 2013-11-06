/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 LEMOS, University Rostock
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "weightedArithmeticMeanFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "OFstream.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(weightedArithmeticMeanFilter, 0);
addToRunTimeSelectionTable(LESfilter, weightedArithmeticMeanFilter, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
weightedArithmeticMeanFilter::weightedArithmeticMeanFilter
(
    const fvMesh& mesh,
    scalar alpha,
    scalar beta
)
:
    LESfilter(mesh),
    addressing_(centredCPCCellToCellExtStencilObject::New(mesh)),
    alpha_(alpha),
    beta_(beta)
{}


weightedArithmeticMeanFilter::weightedArithmeticMeanFilter(const fvMesh& mesh, const dictionary& dict)
:
    LESfilter(mesh),
    addressing_(centredCPCCellToCellExtStencilObject::New(mesh)),
    alpha_(readScalar(dict.subDict(type() + "Coeffs").lookup("alpha"))),
    beta_(readScalar(dict.subDict(type() + "Coeffs").lookup("beta")))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void weightedArithmeticMeanFilter::read(const dictionary& dict)
{
    dict.subDict(type() + "Coeffs").lookup("alpha") >> alpha_;
    dict.subDict(type() + "Coeffs").lookup("beta") >> beta_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


Foam::tmp<Foam::volVectorField> Foam::weightedArithmeticMeanFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const 
{
    List<List<vector> > stencilData(mesh().nCells());

    addressing_.collectData
    (
        unFilteredField(),
        stencilData
    );

    tmp<volVectorField> filteredField = unFilteredField();

    forAll(filteredField(), cellI)
    {
        List<vector>& cStencil = stencilData[cellI];

	vector& vvf = filteredField()[cellI];
	vvf = pTraits<vector>::zero;
        
	if(cStencil.size() > 0)
	{
	    label cStencilI = 0;

	    vvf = alpha_ * cStencil[cStencilI++];

            while(cStencilI < cStencil.size())
	    {
	        vvf += beta_ * cStencil[cStencilI++];
	    }
	
	    vvf /= (alpha_ + beta_ * (cStencil.size()-1));
	}
    }
  
    filteredField().correctBoundaryConditions();

    unFilteredField.clear();

    return filteredField();
}

Foam::tmp<Foam::volScalarField> Foam::weightedArithmeticMeanFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    List<List<scalar> > stencilData(mesh().nCells());

    addressing_.collectData
    (
        unFilteredField(),
        stencilData
    );

    tmp<volScalarField> filteredField = unFilteredField();

    forAll(filteredField(), cellI)
    {
        List<scalar>& cStencil = stencilData[cellI];

        scalar& vvf = filteredField()[cellI];
        vvf = 0.0;

        if(cStencil.size() > 0)
        {
            label cStencilI = 0;

            vvf = alpha_ * cStencil[cStencilI++];

            while(cStencilI < cStencil.size())
            {
                vvf += beta_ * cStencil[cStencilI++];
            }

	    vvf /= (alpha_ + beta_ * (cStencil.size()-1));
        }
    }

    filteredField().correctBoundaryConditions();

    unFilteredField.clear();

    return unFilteredField;
}

Foam::tmp<Foam::volSymmTensorField> Foam::weightedArithmeticMeanFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    List<List<symmTensor> > stencilData(mesh().nCells());

    addressing_.collectData
    (
        unFilteredField(),
        stencilData
    );

    tmp<volSymmTensorField> filteredField = unFilteredField();

    forAll(filteredField(), cellI)
    {
        List<symmTensor>& cStencil = stencilData[cellI];

        symmTensor& vvf = filteredField()[cellI];
        vvf = pTraits<symmTensor>::zero;;

        if(cStencil.size() > 0)
        {
            label cStencilI = 0;

            vvf = alpha_ * cStencil[cStencilI++];

            while(cStencilI < cStencil.size())
            {
                vvf += beta_ * cStencil[cStencilI++];
            }

	    vvf /= (alpha_ + beta_ * (cStencil.size()-1));
        }
    }

    filteredField().correctBoundaryConditions();

    unFilteredField.clear();

    return filteredField;
}

Foam::tmp<Foam::volTensorField> Foam::weightedArithmeticMeanFilter::operator()
(
    const tmp<volTensorField>& unFilteredField
) const
{
    List<List<tensor> > stencilData(mesh().nCells());

    addressing_.collectData
    (
        unFilteredField(),
        stencilData
    );

    tmp<volTensorField> filteredField = unFilteredField();

    forAll(filteredField(), cellI)
    {
        List<tensor>& cStencil = stencilData[cellI];

        tensor& vvf = filteredField()[cellI];
        vvf = pTraits<tensor>::zero;;

        if(cStencil.size() > 0)
        {
            label cStencilI = 0;

            vvf = alpha_ * cStencil[cStencilI++];

            while(cStencilI < cStencil.size())
            {
                vvf += beta_ * cStencil[cStencilI++];
            }

	    vvf /= (alpha_ + beta_ * (cStencil.size()-1));
        }
    }

    filteredField().correctBoundaryConditions();

    unFilteredField.clear();


    return filteredField;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
