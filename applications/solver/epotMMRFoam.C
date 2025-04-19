/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    epotMMRFoam


Description
    Transient solver for incompressible, laminar, MHD flow of micropolar fluids including micromagnetorotation (using low-Rem approximation for MHD).
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, laminar MHD flow"
        " of micropolar fluids including MMR."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    //Constitutive Magnetization equation M
    M = (((1.0/(mag(H0))* Mo) * H0) - (((1.0/(mag(H0)))*(Mo*t)) * (H0 ^ N))
    
    //Magnetic flux density term B
    dimensionedVector B = mu * H0;
    
    //Lorentz force term initialization
    volVectorField L = sigma * (-fvc::grad(PotE) ^ B) + sigma * ((U ^ B) ^ B);
    
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"
       
        // Momentum predictor
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu+ku, U)
          - 2.0 * ku * fvc::curl(N)
          - (1.0/rho) * L
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }
        
        //Interpolating cross product u x B over mesh faces
	surfaceScalarField psiub = fvc::interpolate(U ^ B) & mesh.Sf();

	//Poisson equation for electric potential
	fvScalarMatrix PotEEqn
	(
	   fvm::laplacian(PotE)
	   ==
	   fvc::div(psiub)
	);
      
        //Reference potential
	PotEEqn.setReference(PotERefCell, PotERefValue);

	//Solving Poisson equation
	PotEEqn.solve();

	//Computation of current density at cell faces
	surfaceScalarField jn = -(fvc::snGrad(PotE) * mesh.magSf()) + psiub;

	//Current density at face center
	surfaceVectorField jnv = jn * mesh.Cf();

	//Interpolation of current density at cell center
	volVectorField jfinal = fvc::surfaceIntegrate(jnv) - (fvc::surfaceIntegrate(jn) * mesh.C());

	//Update current density distribution and boundary condition
	jfinal.correctBoundaryConditions();

	//Lorentz force computation
	L = sigma * (jfinal ^ B);
        
        // Microrotation equation
           fvVectorMatrix NEqn
           (
               j * fvm::ddt(N)
             + j * fvm::div(phi, N)
             - 2.0 * ku * fvc::curl(U)
             + 4.0 * ku * N
             - fvm::laplacian(gu, N)
             - (1.0/rho) * (M ^ H0); 
        
           NEqn.solve();   

        runTime.write();

	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	<< " ClockTime = " << runTime.elapsedClockTime() << " s"
	<< nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
