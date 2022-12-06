//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// Cleaver - A MultiMaterial Tetrahedral Mesher
// -- Mesher Class
//
//  Author: Jonathan Bronson (bronson@sci.utah.edu)
//
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
//  Copyright (C) 2011, 2012, Jonathan Bronson
//  Scientific Computing & Imaging Institute
//  University of Utah
//
//  Permission is  hereby  granted, free  of charge, to any person
//  obtaining a copy of this software and associated documentation
//  files  ( the "Software" ),  to  deal in  the  Software without
//  restriction, including  without limitation the rights to  use,
//  copy, modify,  merge, publish, distribute, sublicense,  and/or
//  sell copies of the Software, and to permit persons to whom the
//  Software is  furnished  to do  so,  subject  to  the following
//  conditions:
//
//  The above  copyright notice  and  this permission notice shall
//  be included  in  all copies  or  substantial  portions  of the
//  Software.
//
//  THE SOFTWARE IS  PROVIDED  "AS IS",  WITHOUT  WARRANTY  OF ANY
//  KIND,  EXPRESS OR IMPLIED, INCLUDING  BUT NOT  LIMITED  TO THE
//  WARRANTIES   OF  MERCHANTABILITY,  FITNESS  FOR  A  PARTICULAR
//  PURPOSE AND NONINFRINGEMENT. IN NO EVENT  SHALL THE AUTHORS OR
//  COPYRIGHT HOLDERS  BE  LIABLE FOR  ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
//  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
//  USE OR OTHER DEALINGS IN THE SOFTWARE.
//-------------------------------------------------------------------
//-------------------------------------------------------------------

#ifndef CLEAVER_H
#define CLEAVER_H

#include <cwgUtilities.h>
#include "../TSPVideo.h"
#include "../VideoClass.h"

#include <vector>
#include "ScalarField.h"
#include "Volume.h"
#include "TetMesh.h"

#include <string>

namespace Cleaver
{            
    class CleaverMesherImp;
	class CleaverMesherImpLite;

	class CleaverMesherLite
	{
	public:
		CleaverMesherLite(const LabelVolume* tspvideo, double cleaver_tet_ratio);
		~CleaverMesherLite();

		void createTetMesh(bool verbose);
		TetMesh* getTetMesh() const;
		void cleanup();

	private:
		CleaverMesherImpLite* m_pimpl;
	};

	class CleaverMesher
	{
	public:
	    CleaverMesher(const AbstractVolume *volume);
	    ~CleaverMesher();
	
	    void createTetMesh(bool verbose);
	    TetMesh* getTetMesh() const;
	
	    void setVolume(const AbstractVolume *volume);
	    const AbstractVolume* getVolume() const;
	
	    void cleanup();
	
	private:
	    CleaverMesherImp *m_pimpl;
	};

	// THIS IS THE MOST IMPORTANT FUNCTION THAT SHOULD BE USED
	TetMesh* createMeshFromVolume(const AbstractVolume *volume, bool verbose = false);
	TetMesh* createMeshFromTSPVideo(LabelVolume* tsp_video, double cleaver_tet_ratio);
	int	find_tet_storage_index(const Cleaver::Vertex3D* vi, const Cleaver::Tet* tet, bool& found);

	extern const std::string VersionNumber;
	extern const std::string VersionDate;
	extern const std::string Version;
	
	extern const float DefaultAlphaShort;
	extern const float DefaultAlphaLong;
}

#endif // CLEAVER_H
