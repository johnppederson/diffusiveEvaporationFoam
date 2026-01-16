/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    diffusiveEvaporationFoam

Description
    Passive scalar transport equation solver for evaporation simulations.

    \heading Solver details
    The equation is given by:

    \f[
        \ddt{C} + \div \left(\vec{U} C\right) - \div \left(D_C \grad C \right)
        = S_{C}
    \f]

    Where:
    \vartable
        C       | Concentration scalar
        D_C     | Diffusion coefficient
        S_C     | Mass source
    \endvartable

    \heading Required fields
    \plaintable
        C       | Concentration scalar
        U       | Velocity [m/s]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Passive scalar transport equation solver for evaporation simulations."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    const bool adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
    const scalar maxFo = runTime.controlDict().lookupOrDefault<scalar>("maxFo", 0);
    const scalar maxDeltaT = runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

    if (adjustTimeStep)
    {
        Info<< "\nDetermining time-step based on non-zero Fourier number: " << maxFo << endl;
        const scalarField& V = mesh.V();
        const scalar minV = gMin(V);
        const int dimensionality = runTime.controlDict().lookupOrDefault("dimensionality", 3);
        dimensionedScalar minL2(dimArea, 0);
        dimensionedScalar avgL2(dimArea, 0);
        if (dimensionality == 1)
        {
            minL2 = dimensionedScalar(dimArea, Foam::pow(minV, 2.0));
            avgL2 = dimensionedScalar(dimArea, Foam::sum(V*Foam::pow(V, 2.0)) / Foam::sum(V));
        }
        else if (dimensionality == 2)
        {
            minL2 = dimensionedScalar(dimArea, minV);
            avgL2 = dimensionedScalar(dimArea, Foam::sum(V*V) / Foam::sum(V));
        }
        else if (dimensionality == 3)
        {
            minL2 = dimensionedScalar(dimArea, Foam::pow(minV, 2.0/3.0));
            avgL2 = dimensionedScalar(dimArea, Foam::sum(V*Foam::pow(V, 2.0/3.0)) / Foam::sum(V));
        }
        else
        {
            FatalErrorInFunction << "dimensionality = " << dimensionality << " invalid, ";
            FatalErrorInFunction << "must be one of (1, 2, 3)!" << exit(FatalError);
        }
        Info<< "Selected dimensionality: " << dimensionality;
        Info<< ", assuming mesh width along 'empty' axes is 1 m" << endl;
        dimensionedScalar dt = maxFo * minL2 / DC;
        runTime.setDeltaT(Foam::min(dt.value(), maxDeltaT));
        Info<< "New deltaT: " << runTime.deltaTValue() << endl;
        Info<< "Max Fo: " << (DC * runTime.deltaTValue() / minL2).value() << endl;
        Info<< "Avg Fo: " << (DC * runTime.deltaTValue() / avgL2).value() << endl;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix CEqn
            (
                fvm::ddt(C)
              + fvm::div(phi, C)
              - fvm::laplacian(DC, C)
             ==
                fvOptions(C)
            );

            CEqn.relax();
            fvOptions.constrain(CEqn);
            CEqn.solve();
            fvOptions.correct(C);
        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
