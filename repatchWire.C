#include "fvCFD.H"
#include "boundaryMesh.H"
//#include "triSurface.H"
#include "repatchPolyTopoChanger.H"

using namespace Foam;
	
bool aligned (const vector& a,const vector& b,const double tol)
{
    //- Test block
    scalar dotProduct = a & b;
    scalar magValue = mag(1 - dotProduct);

	//if(Foam::debug)
	//{
		//Info << "The normal " << a[0] << " " << a[1] << " " << a[2] << " is aligned with " << b[0] << " " << b[1] << " " << b[2] << endl;
		//Info << "a&b = " << dotProduct << endl;
	    //Info << "mag = "<< magValue << ", tol = " << tol << endl;
	//}

	//Info << "a = " << a[0] << " " << a[1] << " " << a[2] << endl;
	//Info << "b = " << b[0] << " " << b[1] << " " << b[2] << endl;

	//Info << "mag(mag(a[0]/(mag(a) + VSMALL)) = " << mag(mag(a[0]/(mag(a) + VSMALL))) << endl;
	//bool xOk = mag(mag(a[0]/(mag(a) + VSMALL)) - mag(b[0]/(mag(b) + VSMALL))) < tol;
	//Info << "mag(mag(a[1]/(mag(a) + VSMALL)) = " << mag(mag(a[1]/(mag(a) + VSMALL))) << endl;
	//bool yOk = mag(mag(a[1]/(mag(a) + VSMALL)) - mag(b[1]/(mag(b) + VSMALL))) < tol;
	//Info << "mag(mag(a[2]/(mag(a) + VSMALL)) = " << mag(mag(a[2]/(mag(a) + VSMALL))) << endl;
	//bool zOk = mag(mag(a[2]/(mag(a) + VSMALL)) - mag(b[2]/(mag(b) + VSMALL))) < tol;

	//bool xOk = mag(mag(a) - mag(b)) < tol;
	//bool yOk(true);
	//bool zOk(true);
    
    return (magValue < tol);
    //return (xOk && yOk && zOk);
}

void createPatch(boundaryMesh& bMesh, string s1, string s2)
{
    word patchName(s1);
    Info << "patchName created" << endl;
    
    bMesh.addPatch(patchName);
    Info << "patch added" << endl;
    
    bMesh.changePatchType(patchName, s2);
    Info << "patchType changed" << endl;
}

int main(int argc, char* argv[])
{
	
	bool debug = false;
	bool debug2 = false;

	argList args(argc, argv);
	Time runTime
	(
	    "controlDict", 
	    args.rootPath(),
	    args.caseName(),
	    "system",
	    "constant",
	    !args.optionFound("noFunctionObjects")
	 );

    //polyMesh meshP
    fvMesh mesh
	(
	    IOobject
		(
			fvMesh::defaultRegion,
		    runTime.timeName(),
			runTime,
			IOobject::MUST_READ//,
		)
	);

	//fvMesh mesh(meshP);

	Info << "Mesh read" << endl;

	//surfaceVectorField N
    surfaceVectorField N// = mesh.Sf()/mesh.magSf();
	(
	    IOobject
		(
			"N", 
		    runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedVector("zero", 
		                   dimless, vector::zero)
	);
	Info << "N field created" << endl;

	// Calculate normals on the new mesh
    N = mesh.Sf()/mesh.magSf();
	Info << "New normals calculated" << endl;


	//- Create and set boundaryMesh 

    boundaryMesh bMesh;
	Info << "bMesh created" << endl;

	bMesh.read(mesh);
	Info << "boundaryMesh read" << endl;

	//
	// Fill patchIDs with old patchIDs
	//
    labelList patchIDs(bMesh.mesh().size() + bMesh.meshFace()[0], -1);             // Alright
    //DynamicList<label> patchIDs;//(bMesh.mesh().size(), -1); // <- segfaults
	//labelList copyPatchIDs(1,-1);
    //dynamicLabelList patchIDs(copyPatchIDs);//(bMesh.mesh().size(), -1);

	label patchID = 0;
	label nFaces = bMesh.meshFace()[0];
    label face0 = bMesh.meshFace()[0];
	//label nFaces = 0;

	forAll(bMesh.patches(), patchI)
	{
	    forAll(bMesh.patches()[patchI], faceI)
		{
			if (patchI == patchID)
			{
	            patchIDs[nFaces + faceI] = patchID; 

	            if (debug)
				{
				    Info << "patchIDs[" << nFaces + faceI << "] = "
					                    << patchID << endl; 
				}
			}
			else
			{
			    FatalErrorIn("Assign old patchIDs")
				<< "patchI = " << patchI << "while"
				<< "patchID = " << patchID << abort(FatalError);
			}
		}
        nFaces = nFaces + bMesh.patches()[patchI].size();
		patchID++;
	}

	//-b. Init patchIDs; autoPatch
	//DynamicList<label, 1> patchIDs(bMesh.mesh());
	
	//patchIDsNew[bMesh.mesh().size()-1] = -1;
    //labelList patchIDsNew(bMesh.mesh().size(), -1);

    //
    // Fill patchIDs with values for every face by floodfilling without
    // crossing feature edge.
    //

    // Current patch number.

	/*
	//-b. Find faces for the new patch: autoPatch way

    // Find first unset face.
    label unsetFaceI = findIndex(patchIDsNew, -1);

    if (unsetFaceI == -1)
    {}

	else
	{
	*/

	//
	// Create newPatchIDs for mesh faces
	//
	
	
	    //-a. add and rename new patch: my way

	Info << "patchIDs initialized" << endl;
	    //- Add patchName for all faces facing positive Z direction
		createPatch(bMesh, "symmetryPlaneZ1", "symmetryPlane");
		createPatch(bMesh, "symmetryPlaneY1", "symmetryPlane");
		createPatch(bMesh, "wireContact1", "patch");

    label newPatchI = bMesh.patches().size()-1; //  last patchID is newPatchI

        // Fill visited with all faces reachable from unsetFaceI.
        boolList visited(bMesh.mesh().size(), false); // <- If the second arg is not provided, elements will be assigned randomly, and anything but 0 is true.

        //-a. my face marking
        //vector zVector(0, 0, 1);
        
        //-a.2 my forAll, idea taken from autoPatch
		//label count = 0;
	
	/*
		double tol = 1e-5;
		vector nZ (0, 0, -1);
		vector nY (0, -1, 0);


        forAll(bMesh.mesh(), faceI)
        {
		    if // Only split wireContact
			(
			    mesh.boundaryMesh()
				[
				    bMesh.whichPatch(faceI)
				].name() == "wireContact"
			)
			{
        	    //- Calculate face normal

	    	    //-a.2 my-autoPatch way
				face f = bMesh.mesh()[faceI];//.reverseFace();
        	    vector normal = f.normal
        	    (
        	        mesh.points()
        	    );
        
        	    normal /= mag(normal) + VSMALL;

				if(debug)
				{
        	        //Info << "normal = " << normal[0] << " " << normal[1] << " " << normal[2] << endl;
        	        //Info << "mag(normal) = " << mag(normal) << endl;
				}


				//N[face0 + faceI] = normal;
        	    //
        	    // Assign patchIDs
        	    //
        
        	    if(debug)
        	    {
                    //Info << "N set " << endl;
                }

		        if(aligned(normal, nZ, tol)) // newPatchI-2
			    {
        	        //sFaces.append(faceI);
        	    	if(debug)
        	    	{
		                Info <<faceI + face0 <<" " << "( "<< normal[0] << " " << normal[1] << " " << normal[2] << " )"  << endl;
		                //Info << "The normal " << normal[0] << " " << normal[1] << " " << normal[2] << " is aligned with " << nZ[0] << " " << nZ[1] << " " << nZ[2] << endl;
						//Info << "mag(normal) = " << mag(normal) << endl;
                        //Info << "(mag(1 - (normal & vector(0,0,1))) < 1e-5) is true" << endl;
                        //Info << "ma" << endl;
        	    	}
        
	    	    	//-a.2 my-autoPatch way
	    	    	//visited[faceI] = true;
			    	patchIDs[faceI + face0] = newPatchI-2;
        
        	    	if(debug)
        	    	{
			    	    //Info << "patchIDs[" << faceI + face0 << "] = " << patchIDs[ faceI + face0  ] << endl;
			    		//Info << "visited[" << faceI << "] = " << visited[faceI] << endl;
        	    	    //Info << "patchIDs[" << visitedFaces + faceI<< "] = " << newPatchID << endl;
        	    	}
			    }

		        else if(aligned(normal, nY, tol)) // newPatchI-1
			    {
        	        //sFaces.append(faceI);
        	    	if(debug)
        	    	{
		                Info <<faceI + face0 <<" " << "( "<< normal[0] << " " << normal[1] << " " << normal[2] << " )"  << endl;
                        //Info << "(mag(1 - (normal & vector(0,0,1))) < 1e-5) is true" << endl;
		                //Info << "The normal " << normal[0] << " " << normal[1] << " " << normal[2] << " is aligned with " << nY[0] << " " << nY[1] << " " << nY[2] << endl;
						//Info << "mag(normal) = " << mag(normal) << endl;
        	    	}
        
	    	    	//-a.2 my-autoPatch way
	    	    	//visited[faceI] = true;
			    	patchIDs[faceI + face0] = newPatchI-1;
        
        	    	if(debug)
        	    	{
			    	    //Info << "patchIDs[" << faceI + face0 << "] = " << patchIDs[faceI + face0]<< endl;
			    		//Info << "visited[" << faceI << "] = " << visited[faceI] << endl;
        	    	    //Info << "patchIDs[" << visitedFaces + faceI<< "] = " << newPatchID << endl;
        	    	}
			    }
			    else // newPatchI
			    {
        	        //sFaces.append(faceI);
        	    	if(debug)
        	    	{
                        //Info << "(mag(1 - (normal & vector(0,0,1))) < 1e-5) is true" << endl;
						//Info << "Normal alignment did not meet" << endl;
        	    	}
        
	    	    	//-a.2 my-autoPatch way
	    	    	//visited[faceI] = true;
			    	patchIDs[faceI + face0] = newPatchI;
        
        	    	if(debug)
        	    	{
			    	    //Info << "patchIDs[" << faceI + face0 << "] = " << patchIDs[faceI + face0] << endl;
			    		//Info << "visited[" << faceI << "] = " << visited[faceI] << endl;
        	    	    //Info << "patchIDs[" << visitedFaces + faceI<< "] = " << newPatchID << endl;
        	    	}
			    }

		        //if(count == bMesh.mesh().size())
			    //{
			    //    Info << "patchIDs.size() = " << patchIDs.size() << endl;
			    //}
		    }
        }
        //Info << "patchIDs.size() outside the loop = " << patchIDs.size() << endl;

	Info << "new patchIDs set" << endl;

	if(debug)
	{
	    label ps = 0 + face0;
	    //nFaces = 0;

	    forAll(bMesh.patches(), patchI)
	    {
	        forAll(bMesh.patches()[patchI], faceI)
	    	{
	    	    Info << "patchIDs [" << ps << "] = " << patchIDs[ps] << endl;
	    		ps++;
	    	}
	    	//nFaces = nFaces + bMesh.patches()[patchI].size();
	    }
	}

	*/

	//
	// Create new patchList
	//
	
   //-c. Copy from autoPatch

    //const PtrList<boundaryPatch>& patches = bMesh.patches();

    // Create new list of patches with old ones first
    //List<polyPatch*> newPatchPtrList(patches.size());
    List<polyPatch*> newPatchPtrList(bMesh.patches().size());


	forAll(bMesh.patches(), patchI)
	{
	    newPatchPtrList[patchI] = new polyPatch
		(
	        bMesh.patches()[patchI].name(),
	        //bMesh.patches()[bMesh.patches().size()-1].size(),
	        bMesh.patches()[patchI].size(), // new patch added with size zero
	        //bMesh.patches()[bMesh.patches().size()-1].start(),
	        bMesh.patches()[patchI].start()+face0,
	        patchI,
	        mesh.boundaryMesh()
		);
        if (debug)
		{
	        //newPatchPtrList[patchI] = new polyPatch
		    //(
	            Info << "patchName = " << bMesh.patches()[patchI].name() << " "; 
	            Info << "patchSize = " << bMesh.patches()[patchI].size() << " ";// new patch added with size zero
	            //bMesh.patches()[bMesh.patches().size()-1].start(),
	            Info << "patchStart = " << bMesh.patches()[patchI].start()+face0 << " ";
	            Info << "patchIndex = " << patchI << endl;
	            //mesh.boundaryMesh()
		    //);
		}
	}

	

	if(debug2)
	{
	    Info << "bMesh.meshFace()" << endl;
	    forAll(bMesh.meshFace(), I)
	    {
	        Info << "bMesh.meshFace()[" << I << "] = " << bMesh.meshFace()[I] << endl;
	    	//Info << bMesh.meshFace() << endl;
	    }
	}
		
	if(debug)
	{
	    forAll(mesh.boundaryMesh(), patchI)
	    {
	        Info << "patch "<< mesh.boundaryMesh()[patchI].name() << " had " << mesh.boundaryMesh()[patchI].size()<< " faces" << endl;
	    }
	}

	//
	// Create a new mesh and use addPatches, or
	// use the patcher object to change patches of the current mesh
	//

	//newMesh.addPatches(patchList);

	repatchPolyTopoChanger patcher(mesh);
	patcher.changePatches(newPatchPtrList);

	//-b. face ordering: autoPatch
    const labelList& meshFace = bMesh.meshFace();

	if(patchIDs.size() <=0)
	{
	    FatalErrorIn("patchIDs.size() is invalid") << abort(FatalError);
	}
	else
	{
	    Info << "patchIDs.size() before repatch is " << patchIDs.size() << endl;
	}

    for(label faceI=face0; faceI < patchIDs.size(); ++faceI)
    //forAll(patchIDs, faceI)
    {

        label meshFaceI = meshFace[faceI-face0];
        //label meshFaceI = meshFace[faceI];

		if(debug2)
		{
			Info << "********************************************************" << endl;
		    Info << "meshFaceI = " << meshFace[faceI-face0] << endl;
		    //Info << "meshFaceI = " << meshFace[faceI] << endl;
		}

		if(debug2)
		{
		    Info << "patchIDs["<< faceI << "] changed from " << patchIDs[faceI] << " to ";
		}


        //if (
		patcher.changePatchID(meshFaceI, patchIDs[faceI]); // <- Not convertible to bool
		//);
		//{
		//    FatalErrorIn("patcher.changePatchID(meshFaceI, patchIDs[faceI])")<<
		//	                                       << ("Cannot set new patchIDs")
		//										   << abort(FatalError);
		//}
	    
		if(debug2)
		{
		    //Info << "face " << faceI << " set to patch " << 
			Info << patchIDs[faceI] << endl;
		}
    }

	patcher.repatch();
	Info << "newPatchIDs assigned" << endl;
	Info << "patchify done" << endl;
	
	if(debug)
	{
	    forAll(mesh.boundaryMesh(), patchI)
	    {
	        Info << "patch "<< mesh.boundaryMesh()[patchI].name() << " has " << mesh.boundaryMesh()[patchI].size()<< " faces"<< endl;
	    }
	}


	runTime++;
	Info << "Time step advanced" << endl;

	mesh.write();
	Info << "newMesh written in " << runTime.timeName() << endl;
}
