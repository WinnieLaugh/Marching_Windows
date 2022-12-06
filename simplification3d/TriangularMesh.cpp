#include "stdafx.h"
#include "TriangularMesh.h"
#include "TetrahedralMesh.h"

using namespace std;
using namespace Cleaver;

extern float Cleaver::GLOBAL_INTERFACE_COLORS[12][3];

bool cwg::TriangularMesh::save_ply( const std::string& fn, cwg::TetrahedralMesh* tetmesh ) const
{
	// THIS METHOD CAN ONLY SAVE INTERNAL BOUNDARIES
	ofstream ofile(fn);
	if (!ofile)
		return false;

	vector<bool> chosen_verts(tetmesh->ResultV.size(), false);
	vector<cv::Vec3b> TrisColors(this->Tris.size());

	int intersurf_count = 0;
	for (int i=0; i<this->Tris.size(); i++)
	{
		
			int t0 = this->Tris[i]->Tets[0];
			int t1 = this->Tris[i]->Tets[1];
			int lb0 = t0 != -1 ? tetmesh->ResultT[t0].get_label() : -1;
			int lb1 = t1 != -1 ? tetmesh->ResultT[t1].get_label() : -1;

			double dr = GLOBAL_INTERFACE_COLORS[(lb0+12) % 12][0] + GLOBAL_INTERFACE_COLORS[(lb1+12) % 12][0];
			double dg = GLOBAL_INTERFACE_COLORS[(lb0+12) % 12][1] + GLOBAL_INTERFACE_COLORS[(lb1+12) % 12][1];
			double db = GLOBAL_INTERFACE_COLORS[(lb0+12) % 12][2] + GLOBAL_INTERFACE_COLORS[(lb1+12) % 12][2];
			TrisColors[i] = cv::Vec3b(db*127.5, dg*127.5, db*127.5);

			chosen_verts[ this->Tris[i]->Vert(0) ] = true;
			chosen_verts[ this->Tris[i]->Vert(1) ] = true;
			chosen_verts[ this->Tris[i]->Vert(2) ] = true;
		
		
//			if (this->Tris[i]->interSurf)
			{
				++intersurf_count;
			}
	}

	int nchosen_verts = 0;
	for (int i=0; i<chosen_verts.size(); i++)
		nchosen_verts += (chosen_verts[i] ? 1 : 0);

	//write_ply_file_head(ofile, nchosen_verts, this->Tris.size());
	write_ply_file_head(ofile, nchosen_verts, intersurf_count);
	int new_index = 0;
	map<int, int> VIndexMap;
	for (int i=0; i<chosen_verts.size(); i++)
	{
		if (chosen_verts[i])
		{
			ofile << tetmesh->ResultV[i].pos().x << " " 
				<< tetmesh->ResultV[i].pos().y << " " 
				<< tetmesh->ResultV[i].pos().z << " ";
			ofile << (unsigned int)(tetmesh->ResultV[i].color().x*255.0) << " "
				<< (unsigned int)(tetmesh->ResultV[i].color().y*255.0) << " " 
				<< (unsigned int)(tetmesh->ResultV[i].color().z*255.0) << endl;

			VIndexMap[i] = new_index++;
		}
	}

	for (int i=0; i<this->Tris.size(); i++)
	{
//		if (this->Tris[i]->interSurf)
		{
			int new_tids[3] = {0};
			for (int j=0; j<3; j++)
				new_tids[j] = VIndexMap.at( this->Tris[i]->Vert(j) );

			ofile << "3 "
				<< new_tids[0] << " "
				<< new_tids[1] << " "
				<< new_tids[2] << " "
				<< (unsigned int)TrisColors[i][2] << " "
				<< (unsigned int)TrisColors[i][1] << " "
				<< (unsigned int)TrisColors[i][0] << endl;
		}
		
	}

	ofile.close();
	return true;
}

cwg::TriangularMesh::~TriangularMesh()
{
	delete_vector_elems(Tris);
}

