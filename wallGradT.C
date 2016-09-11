/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Application
    wallGradT

Description
    Calculates and writes the gradient of T at the wall. Derived from wallGradU
    by Hasan Celik. This code is tested on OpenFOAM version 3.0.1.
    
    2016, Hamamatsu, Japan

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject Theader
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check T exists
        if (Theader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading T" << endl;
            volScalarField T(Theader, mesh);

            Info<< "    Calculating wallGradT" << endl;

            volScalarField wallGradT
            (
                IOobject
                (
                    "wallGradT",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "wallGradT",
                    T.dimensions()/dimLength,
                    0
                )
            );

            const fvPatchList& patches = mesh.boundary();

            forAll(wallGradT.boundaryField(), patchi)
            {
                const fvPatch& currPatch = patches[patchi];

                if (isA<wallFvPatch>(currPatch))
                {
                    wallGradT.boundaryField()[patchi] =
                        -T.boundaryField()[patchi].snGrad();
                }
            }

            wallGradT.write();
        }
        else
        {
            Info<< "    No T" << endl;
        }
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
