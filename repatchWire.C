/**
* @file repatchWire.C
* @brief A test tool to re-patch an existing mesh based on some criteria;      
* to be embeded in bekaert-v2.0.                                        

* @author Maalik (ali@tensorfields.com)
* @date 17 Sep 2024
*/

#include "fvCFD.H"
#include "bMesh.H"
#include "boundaryMesh.H"
//#include "triSurface.H"
#include "repatchPolyTopoChanger.H"

using namespace Foam;
	
/**
 * @brief Find if the two vector are aligned together in the specified tolerance
 * @param a A vector
 * @param b Another vector
 * @param tol Tolerance
 * @return True if a and b are aligned
 */

typedef bMesh bMeshType;

class createPatch//(boundaryMesh& bMesh, string s1, string s2)
{
    word patchName_;
    word patchType_;
    int patchID_;
	vector criterion_;

    bool aligned (const vector& a,const vector& b,const double tol)
    {
        // Test block
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

	public:
        /** Constructors */
        createPatch
        (
        	boundaryMesh& bMesh,
            const word& s1,
        	const word& s2,
        	const vector& v = vector(0, 0, 0)
        ):
        patchName_(s1),
        patchType_(s2),
        patchID_(bMesh.patches().size()),
        criterion_(v)
        {
            bMesh.addPatch(patchName_);
            bMesh.changePatchType(patchName_, patchType_);
        }
		
		/*
        createPatch
        (
        	boundaryMesh& bMesh,
            const word& s1,
        	const word& s2,
        	const vector& v = vector(0, 0, 0)
        ):
        createPatch
        (
        	const boundaryMesh& bMesh,
            const word& s1,
        	const word& s2
        ),
        criterion_(v)
        {}
		*/
        
        int patchID()
        {
            return patchID_;
        }
        
        bool criterion(const vector& v, const double tol)
        {
            return aligned(criterion_, v, tol);
        }
};

/** Define criterion class */
/*
class criterion: public boundaryPatch 
{
	private:
		bool holds_;
	    int newPatchID_;
		boundaryPatch& bp_;

    public:
		inline bool holds()
		{return holds_;}

	    inline int newPatchID()
		{return newPatchID_;}

	    criterion(boundaryPatch& bp, int newPatchID, bool holds=false):
		holds_(holds),
		newPatchID_(newPatchID),
		boundaryPatch_(bp)
		{}
};
*/

/**
 * @brief Driver:
 */
int main(int argc, char* argv[])
{
	
	bool debug = true;   
	bool debug2 = true;  

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
	/**
     * \anchor meshAnch 
	 * Create an `fvMesh` object called `mesh` 
	 *
	 *     |
	 *     V
	 */
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

	/**
	 * @brief create a `surfaceVectorField` called `N`
	 *
	 *       |
	 *       V
	 */
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
		/**
		 * @brief Initialize `N` with **O**
		 *
		 *       |
		 *       V
		 */
		dimensionedVector("zero", 
		                   dimless, vector::zero)
	);
	Info << "N field created" << endl;

	/**
	 * @brief Calculate `N` (face normals) for \ref meshAnch "\c mesh"
	 *
	 *        |
	 *        V
	 */
    N = mesh.Sf()/mesh.magSf();
	Info << "New normals calculated" << endl;


	/** 
	 * @brief Create \c boundaryMesh object called \c bMesh
	 *
	 *       |
	 *       V
	 */

    boundaryMesh bMesh;
	Info << "bMesh created" << endl;

	/** 
	 * @brief Init @c bMesh from @c mesh.
	 *
	 *        |
     *        V
	 */
	bMesh.read(mesh);
	Info << "boundaryMesh read" << endl;

    /**
	 * @brief Init a <b>\c labelList</b> object called \anchor patchIDs_anchor <b>\c patchIDs</b>
	 * with all elements being invalid (-1).
	 *
	 * This list is meant to correspond each face with the correct patch. Init with the size of the sum of two sizes:
	 */

    labelList patchIDs(
	                  /**
	                  * - boundary size of <b>\c mesh</b>, 
                      *
					  * \note <b>\c bMesh.mesh().size()</b> 
	                  * where <b>\c bMesh.mesh()</b> returns a <b>\c bMesh</b> object which includes all boundary data.
					  * 
				      * \note If the size of the boundary is, say, 10, this
					  * value is therefore 10.
	                  */
	                      bMesh.mesh().size() 
	                  /** 
					  * - label of the first boundary face in the original \ref meshAnch "\c mesh". 

					  * \note Due to the counting convention of C++, this
					  * value is actually the number of the \c internalFace s
                      * . So If the number of the \c internalFace s is 10,
					  * this 
					  * value is also 10 (the 11th face in the mesh).

					  * \note In such a scenario, the original \ref meshAnch "\c
					  * mesh" has had 20
					  * faces, and \ref patchIDs_anchor "\c patchIDs" also
					  * has 20 elements.

					  * \note The size of the \ref patchIDs_anchor "\c
					  * patchIDs" equals the number of faces in the \ref
					  * meshAnch "\c mesh"
					  */
	                 + bMesh.meshFace()[0], -1); 
					 /**
					 \varbatim
					 |
					 V
					 \endverbatim
					 */

    //DynamicList<label> patchIDs;//(bMesh.mesh().size(), -1); // <- segfaults
	//labelList copyPatchIDs(1,-1);
    //dynamicLabelList patchIDs(copyPatchIDs);//(bMesh.mesh().size(), -1);

	/**
	 * \brief Init a few <b>\c label</b>s, namely
	 *
	 * \anchor patchID_anchor*
	 * Set patchID to 0
	 *     |
	 *     V
	 */
	label patchID = 0;
	label nInternalFaces = bMesh.meshFace()[0];
	label nVisitedFaces = bMesh.meshFace()[0];
	//label nFaces = 0;

	/** 
	 * Init \c patchIDs with old face-patch correspondence
	 *
	 * For each patch in the mesh (\ref meshAnch "\c mesh")
	 *    |
	 *    V
     */
	forAll(bMesh.patches(), patchI)
	{
		/**
		 * For each face in the patch
		 *     |
		 *     V
		 */
	    forAll(bMesh.patches()[patchI], faceI)
		{ 
			/**
			 * Check if the outer loop counter (patch counter) is patchID (\ref patchID_anchor "\c patchID")
			 *    |
			 *    V
			 */
		    if (patchI == patchID)
			{
				/**
				 * The patchIDs elements corresponding to \c internalFaces remains invalid, and start changing patchIDs
				 * for the face 0 of the patch 0 until completing patch 0, then face0 of patch 1 until completion, etc. 
				 *     |
				 *     V
				 */
	            patchIDs[nVisitedFaces + faceI] = patchID; 

	            if (debug)
				{
				    Info << "patchIDs[" << nInternalFaces + faceI << "] = "
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
	    /** 
		 * update \c visitedFaces
		 * 
		 * Now you have counted all faces in the patch \c patchI
		 * \note \c visitedFaces is initially \c nInternalFaces 
		 * 
		 * \note \c size + \c index = \c index,
		 * where size in this case is the number of faces, and index
		 * is face index starting from 0
		 *     |
		 *     V
		 */
         nVisitedFaces = nVisitedFaces + bMesh.patches()[patchI].size();
		 /** Increment patchID 
		  *     |
		  *     V
		  */
		 patchID++;
	}
	//DynamicList<label, 1> patchIDs(bMesh.mesh());
	
	//
	// Create newPatchIDs for mesh faces
	//

	/**
	 * newPatchI created: The index of the last patch in \c bMesh
	 */
    label newPatchI = bMesh.patches().size()-1; //  
	
	/**
	 * Add new patches with zero faces to \c bMesh
	 */
	Info << "patchIDs initialized" << endl;
	    //- Add patchName for all faces facing positive Z direction
		createPatch symmetryPlaneZ1(bMesh, word("symmetryPlaneZ1"), word("symmetryPlane"), vector(0, 0, -1));
		createPatch symmetryPlaneY1(bMesh, word("symmetryPlaneY1"), word("symmetryPlane"), vector(0, -1, 0));
		createPatch    wireContact1(bMesh, word("wireContact1"),    word("patch"));

        // Fill visited with all faces reachable from unsetFaceI.
		/** Create \c boolList \c visited with the size of \c bMesh*/
        boolList visited(bMesh.mesh().size(), false); // <- If the second arg is not provided, elements will be assigned randomly, and anything but 0 is true.
	
		const double tol = 1e-5;

		/**
		 * Create the list of the boundary faces 
		 */
		const bMeshType& boundaryFaces = bMesh.mesh();
		const pointField& meshPoints = mesh.points();
		const PtrList<boundaryPatch> patches = bMesh.patches();

		/** Loop over boundary faces of *wireContact* to correspond them to an appropriate new patch based on their normal direction */
        forAll(boundaryFaces, faceI)
        {
		    if // Only split wireContact
			(
			    patches
				[
					/** return patch index */
				    bMesh.whichPatch(faceI)
				].name() == "wireContact"
			)
			{
				/** Start criterion;
				 * 
				 * Criterion can be a self-contained, general class..?
				 */

        	    //- Calculate face normal

				face f = boundaryFaces[faceI];//.reverseFace();
        	    vector normal = f.normal(meshPoints);
        
        	    normal /= mag(normal) + VSMALL;

		        if(symmetryPlaneZ1.criterion(normal, tol))
			    {
	    	    	//visited[faceI] = true;
			    	patchIDs[faceI + nInternalFaces] = symmetryPlaneZ1.patchID();
        
        	    	if(debug)
        	    	{
			    	    //Info << "patchIDs[" << faceI + face0 << "] = " << patchIDs[ faceI + face0  ] << endl;
			    		//Info << "visited[" << faceI << "] = " << visited[faceI] << endl;
        	    	    //Info << "patchIDs[" << visitedFaces + faceI<< "] = " << newPatchID << endl;
        	    	}
			    }
				/*

		        else if(aligned(normal, nY, tol)) // newPatchI-1
			    {
        	        //sFaces.append(faceI);
        	    	if(debug)
        	    	{
		                Info <<faceI + nInternalFaces <<" " << "( "<< normal[0] << " " << normal[1] << " " << normal[2] << " )"  << endl;
                        //Info << "(mag(1 - (normal & vector(0,0,1))) < 1e-5) is true" << endl;
		                //Info << "The normal " << normal[0] << " " << normal[1] << " " << normal[2] << " is aligned with " << nY[0] << " " << nY[1] << " " << nY[2] << endl;
						//Info << "mag(normal) = " << mag(normal) << endl;
        	    	}
        
	    	    	//-a.2 my-autoPatch way
	    	    	//visited[faceI] = true;
			    	patchIDs[faceI + nInternalFaces] = patchI_symmY;
        
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
			    	patchIDs[faceI + nInternalFaces] = patchI_wireContact;
        
        	    	if(debug)
        	    	{
			    	    //Info << "patchIDs[" << faceI + face0 << "] = " << patchIDs[faceI + face0] << endl;
			    		//Info << "visited[" << faceI << "] = " << visited[faceI] << endl;
        	    	    //Info << "patchIDs[" << visitedFaces + faceI<< "] = " << newPatchID << endl;
        	    	}
			    }
				*/

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
	/*
	    label ps = 0 + nInternalFaces;
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
	*/
	}

	

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
	        bMesh.patches()[patchI].start()+nInternalFaces,
	        patchI,
	        mesh.boundaryMesh()
		);
        if (debug)
		{
		/*
	        //newPatchPtrList[patchI] = new polyPatch
		    //(
	            Info << "patchName = " << bMesh.patches()[patchI].name() << " "; 
	            Info << "patchSize = " << bMesh.patches()[patchI].size() << " ";// new patch added with size zero
	            //bMesh.patches()[bMesh.patches().size()-1].start(),
	            Info << "patchStart = " << bMesh.patches()[patchI].start()+nInternalFaces << " ";
	            Info << "patchIndex = " << patchI << endl;
	            //mesh.boundaryMesh()
		    //);
		*/
		}
	}

	

	if(debug2)
	{
	/*
	    Info << "bMesh.meshFace()" << endl;
	    forAll(bMesh.meshFace(), I)
	    {
	        Info << "bMesh.meshFace()[" << I << "] = " << bMesh.meshFace()[I] << endl;
	    	//Info << bMesh.meshFace() << endl;
	    }
	*/
	}
		
	if(debug)
	{
	/*
	    forAll(mesh.boundaryMesh(), patchI)
	    {
	        Info << "patch "<< mesh.boundaryMesh()[patchI].name() << " had " << mesh.boundaryMesh()[patchI].size()<< " faces" << endl;
	    }
	*/
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

    for(label faceI=nInternalFaces; faceI < patchIDs.size(); ++faceI)
    //forAll(patchIDs, faceI)
    {

        label meshFaceI = meshFace[faceI-nInternalFaces];
        //label meshFaceI = meshFace[faceI];

		if(debug2)
		{
			Info << "********************************************************" << endl;
		    Info << "meshFaceI = " << meshFace[faceI-nInternalFaces] << endl;
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
