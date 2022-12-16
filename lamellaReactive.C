/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    scalarTransportFoam

Description
    Solves the steady or transient transport equation for a passive scalar.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "createFvOptions.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

  

    vector x(1,0,0);
    const volScalarField& xCoord = mesh.C().component(vector::X);

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix AEqn
            (
                fvm::ddt(A)
	      - fvm::laplacian(DA,A)
	      + fvm::Sp(k*B,A)
	      - mesh.C().component(vector::X)*gamma*(x&fvc::grad(A))
             ==
                fvOptions(A)
            );

            AEqn.relax();
            fvOptions.constrain(AEqn);
            AEqn.solve();
            fvOptions.correct(A);

	    fvScalarMatrix BEqn
            (
                fvm::ddt(B)
	      - fvm::laplacian(DB,B)
	      + fvm::Sp(k*A,B)
	      - mesh.C().component(vector::X)*gamma*(x&fvc::grad(B))
             ==
                fvOptions(B)
            );

            BEqn.relax();
            fvOptions.constrain(BEqn);
            BEqn.solve();
            fvOptions.correct(B);

	    fvScalarMatrix CEqn
            (
                fvm::ddt(C)
              - k*A*B
            );

            CEqn.relax();
            fvOptions.constrain(CEqn);
            CEqn.solve();
            fvOptions.correct(C);
        }
	Rate = fvc::ddt(C);

        runTime.write();
	
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
