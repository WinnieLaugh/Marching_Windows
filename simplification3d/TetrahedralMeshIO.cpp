#include "stdafx.h"
#include "TetrahedralMesh.h"
#include "Utils/MemorySurveillance.h"
#include <algorithm>

using namespace cv;
using namespace std;
using namespace Cleaver;

extern float Cleaver::GLOBAL_INTERFACE_COLORS[12][3];


std::vector<int> find_neighbour(std::vector<std::vector<int>>& indexes, int query_index) {
	std::vector<int> result;
	for (std::size_t i = 0; i < 4; i++) {
		std::vector<int> vs;
		vs.push_back(indexes[query_index][(i + 1) % 4]);
		vs.push_back(indexes[query_index][(i + 2) % 4]);
		vs.push_back(indexes[query_index][(i + 3) % 4]);

		for (std::size_t j = 0; j < indexes.size(); j++) {
			if (j != query_index) {
				auto findIter = std::find(indexes[j].begin(), indexes[j].end(), vs[0]);

				if (findIter != indexes[j].end()) {
					findIter = std::find(indexes[j].begin(), indexes[j].end(), vs[1]);

					if (findIter != indexes[j].end()) {
						findIter = std::find(indexes[j].begin(), indexes[j].end(), vs[2]);

						if (findIter != indexes[j].end()) {
							result.push_back(j);
							break;
						} // matches 3 vertices
					}// matches 2 vertices
				}// matches 1 vertex
			}// query 8 tetraheda
		}// the four indices in the querying tetrahedon
	}

	if (result[0] == result[1]){
		std::cout << query_index << std::endl;
		std::cout << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
		std::cout << indexes[query_index][0] << " " << indexes[query_index][1] << " " << indexes[query_index][2] << " " << indexes[query_index][3] << std::endl;
	}

	for(int i=0; i<4; i++){
		if(result[i] < 0)
			std::cout << result[i] << std::endl;
	}
	return result;
}



bool check_valid(Cleaver::vec3 &p, Cleaver::vec3 &q, Cleaver::vec3 &r, Cleaver::vec3 &s) 
  {
      double px=p.x, py=p.y, pz=p.z, qx=q.x, qy=q.y, qz=q.z, rx=r.x, ry=r.y, rz=r.z, sx=s.x, sy=s.y, sz=s.z;

      double pqx = qx - px;
      double pqy = qy - py;
      double pqz = qz - pz;
      double prx = rx - px;
      double pry = ry - py;
      double prz = rz - pz;
      double psx = sx - px;
      double psy = sy - py;
      double psz = sz - pz;

      // CGAL::abs uses fabs on platforms where it is faster than 
      // Then semi-static filter.

      double maxx = CGAL::abs(pqx);
      double maxy = CGAL::abs(pqy);
      double maxz = CGAL::abs(pqz);

      double aprx = CGAL::abs(prx);
      double apsx = CGAL::abs(psx);

      double apry = CGAL::abs(pry);
      double apsy = CGAL::abs(psy);

      double aprz = CGAL::abs(prz);
      double apsz = CGAL::abs(psz);

      if (maxx < aprx) maxx = aprx;
      if (maxx < apsx) maxx = apsx;
      if (maxy < apry) maxy = apry;
      if (maxy < apsy) maxy = apsy;
      if (maxz < aprz) maxz = aprz;
      if (maxz < apsz) maxz = apsz;

      double det = CGAL::determinant(pqx, pqy, pqz,
                                     prx, pry, prz,
                                     psx, psy, psz);

      double eps = 5.1107127829973299e-15 * maxx * maxy * maxz;

      if (det > eps)  return true;
      if (det < -eps) return false;

  }

bool cwg::TetrahedralMesh::save_txt_tetm( const std::string& filename ) const
{
	throw exception("save_txt_tetm: not be tested!");
	//create_path_from_full_path(filename);
	ofstream ofile(filename);
	if (!ofile)
		return false;

	ofile << m_storage_verts.size() << " " << m_storage_tets.size() << endl;
	for (int i=0; i<m_storage_verts.size(); i++)
	{
		
		ofile << m_storage_verts[i] << endl;
	}
	for (int i=0; i<m_storage_tets.size(); i++)
	{
		
		ofile << m_storage_tets[i] << endl;
	}

	ofile.close();
	return true;
}

bool cwg::TetrahedralMesh::load_txt_tetm( const std::string& filename )
{
	throw exception("load_txt_tetm: not be tested!");
	release();
	ifstream ifile;
	ifile.open(filename);
	if (ifile.is_open() == false)
		return false;

	int nverts, ntets;
	ifile >> nverts >> ntets;
	m_storage_verts.resize(nverts);//m_storage_verts = new TetVertex [nverts];
	m_storage_tets.resize(ntets);//m_storage_tets = new Tetrahedron [ntets];

	int id, bordercode;
	double x, y, z, fval, w;
	for (int i=0; i<nverts; i++)
	{
		TetVertex* v = &m_storage_verts[i];//m_storage_verts + i;
		ifile >> *v;
		m_storage_verts[v->ID] = *v;
	}
	int vids[4];
	int label;
	for (int i=0; i<ntets; i++)
	{
		Tetrahedron* t = &m_storage_tets[i];//m_storage_tets + i;
		ifile >> *t;
		m_storage_tets[t->ID] = *t;
		t->correlate_verts(this);
	}
	ifile.close();
	return true;
}


bool cwg::TetrahedralMesh::load_mesh( const std::string& filename )
{
	//throw exception("load_mesh: not be tested!");
	release();
	ifstream ifile;
	ifile.open(filename);
	if (ifile.is_open() == false)
		return false;

	int nverts, ntets;
	vec3 inputV;
	int label;
	
	string txt;
	while (ifile >> txt)
	{
		if(txt == "Vertices")
		{
			ifile >> nverts;
			m_storage_verts.resize(nverts, NULL);
			m_storage_verts.resize(nverts);

			int id, bordercode;
			double x, y, z, fval, w;
			for (int i=0; i<nverts; i++)
			{
				TetVertex* v = &m_storage_verts[i];//m_storage_verts + i;
				ifile >> inputV.x;
				ifile >> inputV.y;
				ifile >> inputV.z;
				ifile >> label;

				v->ID = i;
				v->setPos(inputV);
				
				m_storage_verts[v->ID] = *v;
			}
		}

		else if(txt == "Tetrahedra")
		{
			ifile >> ntets;
			
			m_storage_tets.resize(ntets);

			int vids[4];
			int label;
			for (int i=0; i<ntets; i++)
			{
				Tetrahedron* t = &m_storage_tets[i];//m_storage_tets + i;
				ifile >> t->Verts[0];
				ifile >> t->Verts[1];
				ifile >> t->Verts[2];
				ifile >> t->Verts[3];

				ifile >> label;

				t->Verts[0]--;
				t->Verts[1]--;
				t->Verts[2]--;
				t->Verts[3]--;

				t->ID = i;
				t->set_label(label);
				
				m_storage_tets[t->ID] = *t;
				t->correlate_verts(this);
			}
		}
	}


	
	
	ifile.close();
	return true;
}

bool cwg::TetrahedralMesh::load_mesh_surf( const std::string& filename )
{
	//throw exception("load_mesh: not be tested!");
	release();
	ifstream ifile;
	ifile.open(filename);
	if (ifile.is_open() == false)
		return false;

	int nverts, ntri;
	vec3 inputV;
	int label;

	int v0,v1,v2;

	string txt;
	while (ifile >> txt)
	{
		if(txt == "Vertices")
		{
			ifile >> nverts;
			m_storage_verts.resize(nverts, NULL);
			m_storage_verts.resize(nverts);

			int id, bordercode;
			double x, y, z, fval, w;
			for (int i=0; i<nverts; i++)
			{
				TetVertex* v = &m_storage_verts[i];//m_storage_verts + i;
				ifile >> inputV.x;
				ifile >> inputV.y;
				ifile >> inputV.z;
				ifile >> label;

				v->ID = i;
				v->setPos(inputV);

				m_storage_verts[v->ID] = *v;
			}
		}

		else if(txt == "Triangles")
		{
			ifile >> ntri;
			Tris.resize(ntri, NULL);
			m_storage_tris.resize(ntri);

			int vids[4];
			int label;
			for (int i=0; i<ntri; i++)
			{
				TetTriangle* t = &m_storage_tris[i];//m_storage_tets + i;
				ifile >> v0;
				ifile >> v1;
				ifile >> v2;
				
				t->SetVertices(v0-1,v1-1,v2-1);

				ifile >> label;

				t->ID = i;
				
				Tris[t->ID]=t;
			}
		}
	}




	ifile.close();
	return true;
}

bool cwg::TetrahedralMesh::save_txt_surf_ply( const std::string& filename ) const
{
	ofstream ofile(filename);
	if (!ofile)
		return false;

	// 1. find referenced verts:
	vector<int> SurfTris;
	vector<cv::Vec3b> SurfTris_Color;
	SurfTris.reserve(Tris.size());
	SurfTris_Color.reserve(Tris.size());

	vector<bool> chosen_verts(ResultV.size(), false);
	for (int i=0; i<Tris.size(); i++)
	{
		int t0 = Tris[i]->Tets[0];
		int t1 = Tris[i]->Tets[1];

		bool is_a_surface = false;
		bool is_domain_surface = false;
		cv::Vec3b color;

		if (t0 < 0 && t1 < 0)
			throw exception("bool cwg::TetrahedralMesh::save_txt_surf_ply( const std::string& filename ) const:\
							t0 < 0 && t1 < 0");
		/*else if (t0 >= 0 && t1 < 0)
		{
			is_a_surface = true;
			int lb0 = ResultT[t0].get_label();
			color = cv::Vec3b(	GLOBAL_INTERFACE_COLORS[lb0 % 12][0]*255,
				GLOBAL_INTERFACE_COLORS[lb0 % 12][1]*255,
				GLOBAL_INTERFACE_COLORS[lb0 % 12][2]*255	);
		}
		else if (t1 >= 0 && t0 < 0)
		{
			is_a_surface = true;
			int lb1 = ResultT[t1].get_label();
			color = cv::Vec3b(	GLOBAL_INTERFACE_COLORS[lb1 % 12][0]*255,
				GLOBAL_INTERFACE_COLORS[lb1 % 12][1]*255,
				GLOBAL_INTERFACE_COLORS[lb1 % 12][2]*255	);
		}*/
		else
		{
			int lb0 = ResultT[t0].get_label();
			int lb1 = ResultT[t1].get_label();
			if (lb0 != lb1)
			{
				is_a_surface = true;
				double dr = GLOBAL_INTERFACE_COLORS[lb0 % 12][0] + GLOBAL_INTERFACE_COLORS[lb1 % 12][0];
				double dg = GLOBAL_INTERFACE_COLORS[lb0 % 12][1] + GLOBAL_INTERFACE_COLORS[lb1 % 12][1];
				double db = GLOBAL_INTERFACE_COLORS[lb0 % 12][2] + GLOBAL_INTERFACE_COLORS[lb1 % 12][2];
				color = cv::Vec3b(dr*127.5, dg*127.5, db*127.5);
			}
		}

		if (is_a_surface)
		{
			SurfTris.push_back(i);
			chosen_verts[ Tris[i]->Vert(0) ] = true;
			chosen_verts[ Tris[i]->Vert(1) ] = true;
			chosen_verts[ Tris[i]->Vert(2) ] = true;
			SurfTris_Color.push_back(color);
		}
	}

	int nchosen_verts = 0;
	for (int i=0; i<chosen_verts.size(); i++)
		nchosen_verts += (chosen_verts[i] ? 1 : 0); 

	map<int,int> VIndexMap;
	int new_index = 0;
	write_ply_file_head(ofile, nchosen_verts, SurfTris.size());
	for (int i=0; i<chosen_verts.size(); i++)
	{
		if (chosen_verts[i])
		{
			ofile << ResultV[i].pos().x << " " << ResultV[i].pos().y << " " << ResultV[i].pos().z << " ";
			VIndexMap[i] = new_index++;
			ofile << "192 214 203" << endl;
		}
	}
	for (int i=0; i<SurfTris.size(); i++)
	{
		int new_tids[3] = {0};
		for (int j=0; j<3; j++)
			new_tids[j] = VIndexMap.at( Tris[ SurfTris[i] ]->Vert(j) );

		ofile << "3 " 
			<< new_tids[0] << " "
			<< new_tids[1] << " "
			<< new_tids[2] << " "
			<< (unsigned int)SurfTris_Color[ i ][2] << " "
			<< (unsigned int)SurfTris_Color[ i ][1] << " "
			<< (unsigned int)SurfTris_Color[ i ][0] << endl;
	}
}


bool cwg::TetrahedralMesh::save_txt_vol_ply_weight( const std::string& filename ) const
{
	

	ofstream ofile(filename);
	if (!ofile)
		return false;

	bool vis_vert_bordercode(false), vis_vert_dsvol(true);


	double max_dsvol = -1.0;
	double min_dsvol = numeric_limits<double>::max();
	Mat tet_vis_table = Mat::zeros(ResultT.size(), 1, CV_64FC1);
/*	for (int i=0; i<ResultT.size(); i++)
	{
		double tmp_dsvol = ResultT[i].compute_delta_supervolume(this);
		tet_vis_table.at<double>(i, 0) = tmp_dsvol;
		min_dsvol = min(min_dsvol, tmp_dsvol);
		max_dsvol = max(max_dsvol, tmp_dsvol);
	}
	cout << "min_dsvol: " << min_dsvol << "\t" << "max_dsvol: " << max_dsvol << endl;
	double range_dsvol = max_dsvol - min_dsvol;
*/


	Mat intensities(256, 1, CV_8UC1);
	for (int i=0; i<256; i++)
		intensities.at<unsigned char>(i) = i;
	Mat colormap_intensities;
	cv::applyColorMap(intensities, colormap_intensities, cv::COLORMAP_JET);


	write_ply_file_head(ofile, ResultV.size(), Tris.size());

	//====================================================
	//==========max and min weight========================
	double maxw = 0, minw = 10000;
	for (int i=0; i<ResultV.size(); i++)
	{
		double w = (ResultV[i].weight());
		if (w > maxw)
		{
			maxw = w;
		}
		if (w < minw)
		{
			minw = w;
		}
	}

	cout<<maxw<<" "<<minw<<endl;

	double range = maxw-minw;

	int nn = 12;

	double step = range/nn;

	double diffuse[4];
	

	for (int i=0; i<ResultV.size(); i++)
	{
		diffuse[0] = 0.00f;
		diffuse[1] = 0.00f;
		diffuse[2] = 0.00f;
		diffuse[3] = 1.0f;


		ofile << ResultV[i].pos().x << " " << ResultV[i].pos().y << " " << ResultV[i].pos().z << " ";

		double k = (ResultV[i].weight()) - minw;
		//double k = ResultV[i].weight() - minw;
		

		if (k<step)
		{
			diffuse[1] = nn*k/range;
			diffuse[2] = 1.0;
		}
		else if (k<step*2)
		{
			diffuse[1] = 1.0;
			diffuse[2] = 1.0+nn*(step-k)/range;
		}
		else if (k<step*3)
		{
			diffuse[1] = 1.0;
			diffuse[0] = nn*(k-step*2)/range;
		}
		else
		{
			diffuse[1] = 1.0+nn*(step*3-k)/range;
			diffuse[0] = 1.0;
		}
/*

		if (k<step)
		{
			diffuse[1] = nn*k/range;
			diffuse[0] = 1.0;
		}
		else if (k<step*2)
		{
			diffuse[1] = 1.0;
			diffuse[0] = 1.0+nn*(step-k)/range;
		}
		else if (k<step*3)
		{
			diffuse[1] = 1.0;
			diffuse[2] = nn*(k-step*2)/range;
		}
		else
		{
			diffuse[1] = 1.0+nn*(step*3-k)/range;
			diffuse[2] = 1.0;
		}
*/
		ofile 
			<< (int)(diffuse[0]*255.0) << " "
			<< (int)(diffuse[1]*255.0) << " " 
			<< (int)(diffuse[2]*255.0) << endl;


		

/*		if(Verts[i]->weighted)
		{
			ofile 
				<< 255 << " "
				<< 0 << " " 
				<< 0 << endl;
		}
		else
		{
			ofile 
				<< 0 << " "
				<< 0 << " " 
				<< 255 << endl;
		}
	*/	


	}

	for (int i=0; i<Tris.size(); i++)
	{


		int t0 = Tris[i]->Tets[0];
		int t1 = Tris[i]->Tets[1];

		if (t0 < 0 && t1 < 0)
			throw exception("bool cwg::TetrahedronMesh::save_txt_ply( const std::string& filename ) const:\
							t0 < 0 && t1 < 0");
		else if (t0 >= 0 && t1 < 0)
		{
			int lb0 = ResultT[t0].get_label();
			//if(lb0 != 0)
			{
				ofile << "3 " << Tris[i]->Vert(0) << " " << Tris[i]->Vert(1) << " " << Tris[i]->Vert(2) << " ";
				ofile 
					<< (int)(GLOBAL_INTERFACE_COLORS[lb0 % 12][0]*255.0) << " "
					<< (int)(GLOBAL_INTERFACE_COLORS[lb0 % 12][1]*255.0) << " " 
					<< (int)(GLOBAL_INTERFACE_COLORS[lb0 % 12][2]*255.0) << endl;
			}


		}
		else if (t1 >= 0 && t0 < 0)
		{
			int lb1 = ResultT[t1].get_label();
			//if(lb1 != 0)
			{
				ofile << "3 " << Tris[i]->Vert(0) << " " << Tris[i]->Vert(1) << " " << Tris[i]->Vert(2) << " ";
				ofile 
					<< (int)(GLOBAL_INTERFACE_COLORS[lb1 % 12][0]*255.0) << " "
					<< (int)(GLOBAL_INTERFACE_COLORS[lb1 % 12][1]*255.0) << " " 
					<< (int)(GLOBAL_INTERFACE_COLORS[lb1 % 12][2]*255.0) << endl;
			}

		}
		else
		{

			int lb0 = ResultT[t0].get_label();
			int lb1 = ResultT[t1].get_label();

			//if(lb1!=0 || lb0!=0)
			{

				ofile << "3 " << Tris[i]->Vert(0) << " " << Tris[i]->Vert(1) << " " << Tris[i]->Vert(2) << " ";

				double dr = GLOBAL_INTERFACE_COLORS[lb0 % 12][0] + GLOBAL_INTERFACE_COLORS[lb1 % 12][0];
				double dg = GLOBAL_INTERFACE_COLORS[lb0 % 12][1] + GLOBAL_INTERFACE_COLORS[lb1 % 12][1];
				double db = GLOBAL_INTERFACE_COLORS[lb0 % 12][2] + GLOBAL_INTERFACE_COLORS[lb1 % 12][2];

				ofile << (int)(dr*127.5f) << " " << (int)(dg*127.5f) << " " << (int)(db*127.5f) << endl;
			}

		}

	}
	ofile.close();
	return true;
}


bool cwg::TetrahedralMesh::save_txt_vol_ply_without0_weights( const std::string& filename ) const
{
	ofstream ofile(filename);
	if (!ofile)
		return false;

	bool vis_vert_bordercode(false), vis_vert_dsvol(true);


/*	double max_dsvol = -1.0;
	double min_dsvol = numeric_limits<double>::max();
	Mat tet_vis_table = Mat::zeros(ResultT.size(), 1, CV_64FC1);
	for (int i=0; i<ResultT.size(); i++)
	{
		double tmp_dsvol = ResultT[i].compute_delta_supervolume(this);
		tet_vis_table.at<double>(i, 0) = tmp_dsvol;
		min_dsvol = min(min_dsvol, tmp_dsvol);
		max_dsvol = max(max_dsvol, tmp_dsvol);
	}
	cout << "min_dsvol: " << min_dsvol << "\t" << "max_dsvol: " << max_dsvol << endl;
	double range_dsvol = max_dsvol - min_dsvol;
*/


	Mat intensities(256, 1, CV_8UC1);
	for (int i=0; i<256; i++)
		intensities.at<unsigned char>(i) = i;
	Mat colormap_intensities;
	cv::applyColorMap(intensities, colormap_intensities, cv::COLORMAP_JET);

	
	int ts = 0;
	for (int i=0; i<Tris.size(); i++)
	{
		int t0 = Tris[i]->Tets[0];
		int t1 = Tris[i]->Tets[1];
		if(t0>=0 && t1<0)
		{
			int lb0 = ResultT[t0].get_label();
			if(lb0 != 0)
			{
				ts++;
			}
		}
		else if(t0<0 && t1>=0)
		{
			int lb1 = ResultT[t1].get_label();
			if(lb1 != 0)
			{
				ts++;
			}
		}
		else
		{
			int lb0 = ResultT[t0].get_label();
			int lb1 = ResultT[t1].get_label();
			if(lb0 != 0 || lb1 !=0)
			{
				ts++;
			}
		}
	}

	write_ply_file_head(ofile, ResultV.size(), ts);

	//====================================================
	//==========max and min weight========================
	double maxw = 0, minw = 10000;
	for (int i=0; i<ResultV.size(); i++)
	{
		double w = (ResultV[i].weight());
		if (w > maxw)
		{
			maxw = w;
		}
		if (w < minw)
		{
			minw = w;
		}
	}

	cout<<maxw<<" "<<minw<<endl;

	double range = maxw-minw;

	int nn=18;
	double step = range/nn;

	double diffuse[4];

	vis_vert_bordercode = true;
	for (int i=0; i<ResultV.size(); i++)
	{
		
		diffuse[0] = 0.00f;
		diffuse[1] = 0.00f;
		diffuse[2] = 0.00f;
		diffuse[3] = 1.0f;


		ofile << ResultV[i].pos().x << " " << ResultV[i].pos().y << " " << ResultV[i].pos().z << " ";

		float k = (ResultV[i].weight()) - minw;



		if (k<step)
		{
			diffuse[1] = nn*k/range;
			diffuse[2] = 1.0;
		}
		else if (k<step*2)
		{
			diffuse[1] = 1.0;
			diffuse[2] = 1.0+nn*(step-k)/range;
		}
		else if (k<step*3)
		{
			diffuse[1] = 1.0;
			diffuse[0] = nn*(k-step*2)/range;
		}
		else
		{
			diffuse[1] = 1.0+nn*(step*3-k)/range;
			diffuse[0] = 1.0;
		}

		ofile 
			<< (int)(diffuse[0]*255.0) << " "
			<< (int)(diffuse[1]*255.0) << " " 
			<< (int)(diffuse[2]*255.0) << endl;
		
		
	}

	for (int i=0; i<Tris.size(); i++)
	{
		
		
		int t0 = Tris[i]->Tets[0];
		int t1 = Tris[i]->Tets[1];

		if (t0 < 0 && t1 < 0)
			throw exception("bool cwg::TetrahedronMesh::save_txt_ply( const std::string& filename ) const:\
							t0 < 0 && t1 < 0");
		else if (t0 >= 0 && t1 < 0)
		{
			int lb0 = ResultT[t0].get_label();
			if(lb0 != 0)
			{
				ofile << "3 " << Tris[i]->Vert(0) << " " << Tris[i]->Vert(1) << " " << Tris[i]->Vert(2) << " ";
				ofile 
					<< (int)(GLOBAL_INTERFACE_COLORS[lb0 % 12][0]*255.0) << " "
					<< (int)(GLOBAL_INTERFACE_COLORS[lb0 % 12][1]*255.0) << " " 
					<< (int)(GLOBAL_INTERFACE_COLORS[lb0 % 12][2]*255.0) << endl;
			}

			
		}
		else if (t1 >= 0 && t0 < 0)
		{
			int lb1 = ResultT[t1].get_label();
			if(lb1 != 0)
			{
				ofile << "3 " << Tris[i]->Vert(0) << " " << Tris[i]->Vert(1) << " " << Tris[i]->Vert(2) << " ";
				ofile 
					<< (int)(GLOBAL_INTERFACE_COLORS[lb1 % 12][0]*255.0) << " "
					<< (int)(GLOBAL_INTERFACE_COLORS[lb1 % 12][1]*255.0) << " " 
					<< (int)(GLOBAL_INTERFACE_COLORS[lb1 % 12][2]*255.0) << endl;
			}
			
		}
		else
		{
			
			int lb0 = ResultT[t0].get_label();
			int lb1 = ResultT[t1].get_label();

			if(lb1!=0 || lb0!=0)
			{

				ofile << "3 " << Tris[i]->Vert(0) << " " << Tris[i]->Vert(1) << " " << Tris[i]->Vert(2) << " ";

				double dr = GLOBAL_INTERFACE_COLORS[lb0 % 12][0] + GLOBAL_INTERFACE_COLORS[lb1 % 12][0];
				double dg = GLOBAL_INTERFACE_COLORS[lb0 % 12][1] + GLOBAL_INTERFACE_COLORS[lb1 % 12][1];
				double db = GLOBAL_INTERFACE_COLORS[lb0 % 12][2] + GLOBAL_INTERFACE_COLORS[lb1 % 12][2];

				ofile << (int)(dr*127.5f) << " " << (int)(dg*127.5f) << " " << (int)(db*127.5f) << endl;
			}

		}

	}
	ofile.close();
	return true;
}


bool cwg::TetrahedralMesh::save_txt_vol_stellar_mesh( const std::string& filename ) const
{
	ofstream ofile(filename);
	if (!ofile)
		return false;

	bool vis_vert_bordercode(false), vis_vert_dsvol(true);


	double max_dsvol = -1.0;
	double min_dsvol = numeric_limits<double>::max();
	Mat tet_vis_table = Mat::zeros(ResultT.size(), 1, CV_64FC1);
	for (int i=0; i<ResultT.size(); i++)
	{
		double tmp_dsvol = ResultT[i].compute_delta_supervolume(this);
		tet_vis_table.at<double>(i, 0) = tmp_dsvol;
		min_dsvol = min(min_dsvol, tmp_dsvol);
		max_dsvol = max(max_dsvol, tmp_dsvol);
	}
	cout << "min_dsvol: " << min_dsvol << "\t" << "max_dsvol: " << max_dsvol << endl;
	double range_dsvol = max_dsvol - min_dsvol;



	Mat intensities(256, 1, CV_8UC1);
	for (int i=0; i<256; i++)
		intensities.at<unsigned char>(i) = i;
	Mat colormap_intensities;
	cv::applyColorMap(intensities, colormap_intensities, cv::COLORMAP_JET);




	//write_ply_file_head(ofile, Verts.size(), Tris.size());

	ofile << ResultV.size()<<endl;
	vis_vert_bordercode = true;
	for (int i=0; i<ResultV.size(); i++)
	{
		
		ofile << ResultV[i].pos().x << " " << ResultV[i].pos().y << " " << ResultV[i].pos().z << " " << ResultV[i].get_label(0)<<endl;
		


	}
	ofile << ResultT.size()<<endl;
	for (int i=0; i<ResultT.size(); i++)
	{
		ofile << ResultT[i].Verts[0] << " " << ResultT[i].Verts[1] << " " << ResultT[i].Verts[2] << " " <<ResultT[i].Verts[3]<<" "<< ResultT[i].get_label()<<endl;
	}

	ofile << Tris.size()<<endl;
	for (int i=0; i<Tris.size(); i++)
	{
		ofile << Tris[i]->Vert(0) << " " << Tris[i]->Vert(1) << " " << Tris[i]->Vert(2) << " " << ResultT[Tris[i]->Tets[0]].get_label()<<endl;
	}
	ofile.close();
	return true;
}


bool cwg::TetrahedralMesh::save_txt_vol_ply( const std::string& filename ) const
{
	ofstream ofile(filename);
	if (!ofile)
		return false;

	bool vis_vert_bordercode(false), vis_vert_dsvol(true);


	double max_dsvol = -1.0;
	double min_dsvol = numeric_limits<double>::max();
	Mat tet_vis_table = Mat::zeros(ResultT.size(), 1, CV_64FC1);
	for (int i=0; i<ResultT.size(); i++)
	{
		double tmp_dsvol = ResultT[i].compute_delta_supervolume_sliding_result(this);
		tet_vis_table.at<double>(i, 0) = tmp_dsvol;
		min_dsvol = min(min_dsvol, tmp_dsvol);
		max_dsvol = max(max_dsvol, tmp_dsvol);
	}
	cout << "min_dsvol: " << min_dsvol << "\t" << "max_dsvol: " << max_dsvol << endl;
	double range_dsvol = max_dsvol - min_dsvol;



	Mat intensities(256, 1, CV_8UC1);
	for (int i=0; i<256; i++)
		intensities.at<unsigned char>(i) = i;
	Mat colormap_intensities;
	cv::applyColorMap(intensities, colormap_intensities, cv::COLORMAP_JET);


	

	write_ply_file_head(ofile, ResultV.size(), Tris.size());

	vis_vert_bordercode = true;
	for (int i=0; i<ResultV.size(); i++)
	{

		ofile << ResultV[i].pos().x << " " << ResultV[i].pos().y << " " << ResultV[i].pos().z << " ";
		if (vis_vert_bordercode)
		{
			unsigned char tmp = ResultV[i].border_type();
			switch (tmp)
			{
			case 3:
				ofile << "255 0 0" << endl; break;
			case 2:
				ofile << "0 255 0" << endl; break;
			case 1:
				ofile << "255 255 0" << endl; break;
			case 0:
				ofile << "0 255 255" << endl; break;
			}
		}
		else if (vis_vert_dsvol)
		{
			double avg_dsvol = 0.0;
			for (vector<int>::const_iterator it = ResultV[i].Tets.begin(); it != ResultV[i].Tets.end(); it++)
				avg_dsvol += tet_vis_table.at<double>(*it);
			avg_dsvol /= ResultV[i].Tets.size();
			unsigned char colorindex = (avg_dsvol - min_dsvol) / range_dsvol * 255.0;

			ofile << (unsigned int)colormap_intensities.at<cv::Vec3b>(colorindex)[2] << " "
				<< (unsigned int)colormap_intensities.at<cv::Vec3b>(colorindex)[1] << " "
				<< (unsigned int)colormap_intensities.at<cv::Vec3b>(colorindex)[0] << endl;
		}


	}

	for (int i=0; i<Tris.size(); i++)
	{


		int t0 = Tris[i]->Tets[0];
		int t1 = Tris[i]->Tets[1];

		if (t0 < 0 && t1 < 0)
			throw exception("bool cwg::TetrahedronMesh::save_txt_ply( const std::string& filename ) const:\
							t0 < 0 && t1 < 0");
		else if (t0 >= 0 && t1 < 0)
		{
			int lb0 = ResultT[t0].get_label();
			//if(lb0 != 0)
			{
				ofile << "3 " << Tris[i]->Vert(0) << " " << Tris[i]->Vert(1) << " " << Tris[i]->Vert(2) << " ";
				ofile 
					<< (int)(GLOBAL_INTERFACE_COLORS[lb0 % 12][0]*255.0) << " "
					<< (int)(GLOBAL_INTERFACE_COLORS[lb0 % 12][1]*255.0) << " " 
					<< (int)(GLOBAL_INTERFACE_COLORS[lb0 % 12][2]*255.0) << endl;
			}


		}
		else if (t1 >= 0 && t0 < 0)
		{
			int lb1 = ResultT[t1].get_label();
			//if(lb1 != 0)
			{
				ofile << "3 " << Tris[i]->Vert(0) << " " << Tris[i]->Vert(1) << " " << Tris[i]->Vert(2) << " ";
				ofile 
					<< (int)(GLOBAL_INTERFACE_COLORS[lb1 % 12][0]*255.0) << " "
					<< (int)(GLOBAL_INTERFACE_COLORS[lb1 % 12][1]*255.0) << " " 
					<< (int)(GLOBAL_INTERFACE_COLORS[lb1 % 12][2]*255.0) << endl;
			}

		}
		else
		{

			int lb0 = ResultT[t0].get_label();
			int lb1 = ResultT[t1].get_label();

			//if(lb1!=0 || lb0!=0)
			{

				ofile << "3 " << Tris[i]->Vert(0) << " " << Tris[i]->Vert(1) << " " << Tris[i]->Vert(2) << " ";

				double dr = GLOBAL_INTERFACE_COLORS[lb0 % 12][0] + GLOBAL_INTERFACE_COLORS[lb1 % 12][0];
				double dg = GLOBAL_INTERFACE_COLORS[lb0 % 12][1] + GLOBAL_INTERFACE_COLORS[lb1 % 12][1];
				double db = GLOBAL_INTERFACE_COLORS[lb0 % 12][2] + GLOBAL_INTERFACE_COLORS[lb1 % 12][2];

				ofile << (int)(dr*127.5f) << " " << (int)(dg*127.5f) << " " << (int)(db*127.5f) << endl;
			}

		}

	}
	ofile.close();
	return true;
}


bool cwg::TetrahedralMesh::save_txt_vol_ply_quality(const std::string& filename, float min_angle, float max_angle) const
{
	ofstream ofile(filename);
	if (!ofile)
		return false;

	bool vis_vert_bordercode(false), vis_vert_dsvol(true);

	double max_dsvol = -1.0;
	double min_dsvol = numeric_limits<double>::max();
	Mat tet_vis_table = Mat::zeros(ResultT.size(), 1, CV_64FC1);
	for (int i=0; i<ResultT.size(); i++)
	{
		double tmp_dsvol = ResultT[i].compute_delta_supervolume(this);
		tet_vis_table.at<double>(i, 0) = tmp_dsvol;
		min_dsvol = min(min_dsvol, tmp_dsvol);
		max_dsvol = max(max_dsvol, tmp_dsvol);
	}
	cout << "min_dsvol: " << min_dsvol << "\t" << "max_dsvol: " << max_dsvol << endl;
	double range_dsvol = max_dsvol - min_dsvol;
	Mat intensities(256, 1, CV_8UC1);
	for (int i=0; i<256; i++)
		intensities.at<unsigned char>(i) = i;
	Mat colormap_intensities;
	cv::applyColorMap(intensities, colormap_intensities, cv::COLORMAP_JET);
	write_ply_file_head(ofile, ResultV.size(), Tris.size());

	vis_vert_bordercode = true;
	for (int i=0; i<ResultV.size(); i++)
	{
		ofile << ResultV[i].pos().x << " " << ResultV[i].pos().y << " " << ResultV[i].pos().z << " ";
		if (vis_vert_bordercode)
		{
			unsigned char tmp = ResultV[i].border_type();
			switch (tmp)
			{
			case 3:
				ofile << "255 0 0" << endl; break;
			case 2:
				ofile << "0 255 0" << endl; break;
			case 1:
				ofile << "255 255 0" << endl; break;
			case 0:
				ofile << "0 255 255" << endl; break;
			}
		}
		else if (vis_vert_dsvol)
		{
			double avg_dsvol = 0.0;
			for (vector<int>::const_iterator it = ResultV[i].Tets.begin(); it != ResultV[i].Tets.end(); it++)
				avg_dsvol += tet_vis_table.at<double>(*it);
			avg_dsvol /= ResultV[i].Tets.size();
			unsigned char colorindex = (avg_dsvol - min_dsvol) / range_dsvol * 255.0;

			ofile << (unsigned int)colormap_intensities.at<cv::Vec3b>(colorindex)[2] << " "
				<< (unsigned int)colormap_intensities.at<cv::Vec3b>(colorindex)[1] << " "
				<< (unsigned int)colormap_intensities.at<cv::Vec3b>(colorindex)[0] << endl;
		}
	}

	float range = max_angle-min_angle;

	float step = range/4.0;

	float diffuse[4];


	float diffuse2[4];




	for (int i=0; i<Tris.size(); i++)
	{
		diffuse[0] = 0.00f;
		diffuse[1] = 0.00f;
		diffuse[2] = 0.00f;
		diffuse[3] = 1.0f;

		
		diffuse2[0] = 0.00f;
		diffuse2[1] = 0.00f;
		diffuse2[2] = 0.00f;
		diffuse2[3] = 1.0f;

		ofile << "3 " << Tris[i]->Vert(0) << " " << Tris[i]->Vert(1) << " " << Tris[i]->Vert(2) << " ";
		int t0 = Tris[i]->Tets[0];
		int t1 = Tris[i]->Tets[1];

		if (t0 < 0 && t1 < 0)
			throw exception("bool cwg::TetrahedronMesh::save_txt_ply( const std::string& filename ) const:\
							t0 < 0 && t1 < 0");
		else if (t0 >= 0 && t1 < 0)
		{
			int lb0 = ResultT[t0].get_label();



			/*	ofile 
			<< (int)(GLOBAL_INTERFACE_COLORS[lb0 % 12][0]*255.0) << " "
			<< (int)(GLOBAL_INTERFACE_COLORS[lb0 % 12][1]*255.0) << " " 
			<< (int)(GLOBAL_INTERFACE_COLORS[lb0 % 12][2]*255.0) << endl;
			*/

			float k = ResultT[t0].get_min_angle() - min_angle;
			if (k<step)
			{
				diffuse[1] = 4*k/range;
				diffuse[2] = 1.0;
			}
			else if (k<step*2)
			{
				diffuse[1] = 1.0;
				diffuse[2] = 1.0+4*(step-k)/range;
			}
			else if (k<step*3)
			{
				diffuse[1] = 1.0;
				diffuse[0] = 4*(k-step*2)/range;
			}
			else
			{
				diffuse[1] = 1+4*(step*3-k)/range;
				diffuse[0] = 1.0;
			}

			ofile 
				<< (int)(diffuse[0]*255.0) << " "
				<< (int)(diffuse[1]*255.0) << " " 
				<< (int)(diffuse[2]*255.0) << endl;


		}
		else if (t1 >= 0 && t0 < 0)
		{
			int lb1 = m_storage_tets[t1].get_label();

			float k = m_storage_tets[t1].get_min_angle() - min_angle;
			if (k<step)
			{
				diffuse[1] = 4*k/range;
				diffuse[2] = 1.0;
			}
			else if (k<step*2)
			{
				diffuse[1] = 1.0;
				diffuse[2] = 1.0+4*(step-k)/range;
			}
			else if (k<step*3)
			{
				diffuse[1] = 1.0;
				diffuse[0] = 4*(k-step*2)/range;
			}
			else
			{
				diffuse[1] = 1+4*(step*3-k)/range;
				diffuse[0] = 1.0;
			}

			ofile 
				<< (int)(diffuse[0]*255.0) << " "
				<< (int)(diffuse[1]*255.0) << " " 
				<< (int)(diffuse[2]*255.0) << endl;

		}
		else
		{
			int lb0 = ResultT[t0].get_label();
			int lb1 = ResultT[t1].get_label();

			float k =  ResultT[t0].get_min_angle() - min_angle;
			float k2 = ResultT[t1].get_min_angle() - min_angle;

			if (k<step)
			{
				diffuse[1] = 4*k/range;
				diffuse[2] = 1.0;
			}
			else if (k<step*2)
			{
				diffuse[1] = 1.0;
				diffuse[2] = 1.0+4*(step-k)/range;
			}
			else if (k<step*3)
			{
				diffuse[1] = 1.0;
				diffuse[0] = 4*(k-step*2)/range;
			}
			else
			{
				diffuse[1] = 1+4*(step*3-k)/range;
				diffuse[0] = 1.0;
			}

			if (k2<step)
			{
				diffuse2[1] = 4*k2/range;
				diffuse2[2] = 1.0;
			}
			else if (k2<step*2)
			{
				diffuse2[1] = 1.0;
				diffuse2[2] = 1.0+4*(step-k2)/range;
			}
			else if (k2<step*3)
			{
				diffuse2[1] = 1.0;
				diffuse2[0] = 4*(k2-step*2)/range;
			}
			else
			{
				diffuse2[1] = 1+4*(step*3-k2)/range;
				diffuse2[0] = 1.0;
			}

			/*	ofile 
			<< (int)(diffuse[0]*255.0) << " "
			<< (int)(diffuse[1]*255.0) << " " 
			<< (int)(diffuse[2]*255.0) << endl;

			*/

			double dr = diffuse[0] + diffuse2[0];
			double dg = diffuse[1] + diffuse2[1];
			double db = diffuse[2] + diffuse2[2];

			ofile << (int)(dr*127.5f) << " " << (int)(dg*127.5f) << " " << (int)(db*127.5f) << endl;

		}



	}
	ofile.close();
	return true;
}


bool cwg::TetrahedralMesh::save_obj_quality(const std::string& filename) const
{
	

	vector<int> badTet;
	for (int i=0; i<ResultT.size(); i++)
	{
		if (ResultT[i].get_min_angle() < 0.26)	// < 15 c
		{
			badTet.push_back(i);
		}
	}

	ofstream ofile(filename);
	if (!ofile)
		return false;

	//v
	for (int i=0; i<ResultV.size(); i++)
	{
		ofile << "v "<< ResultV[i].pos().x << " " << ResultV[i].pos().y << " " << ResultV[i].pos().z << endl;
	}

	//tet
	for (int i=0; i<badTet.size(); i++)
	{
		ofile << "f "<< ResultT[badTet[i]].Verts[0] + 1 <<" "<< ResultT[badTet[i]].Verts[1] + 1 <<" "<< ResultT[badTet[i]].Verts[2] + 1 <<endl;
		ofile << "f "<< ResultT[badTet[i]].Verts[0] + 1 <<" "<< ResultT[badTet[i]].Verts[2] + 1 <<" "<< ResultT[badTet[i]].Verts[3] + 1 <<endl;
		ofile << "f "<< ResultT[badTet[i]].Verts[1] + 1 <<" "<< ResultT[badTet[i]].Verts[2] + 1 <<" "<< ResultT[badTet[i]].Verts[3] + 1 <<endl;
		ofile << "f "<< ResultT[badTet[i]].Verts[0] + 1 <<" "<< ResultT[badTet[i]].Verts[1] + 1 <<" "<< ResultT[badTet[i]].Verts[3] + 1 <<endl;
	}

	ofile.close();
}


bool cwg::TetrahedralMesh::save_bin_tetm( const std::string& filename ) const
{
	//create_path_from_full_path(filename);
	ofstream ofile(filename, ios::binary);
	if (!ofile)
		return false;

	int vsz = ResultV.size();
	int tsz = ResultT.size();
	ofile.write((char*)&vsz, sizeof(int));
	ofile.write((char*)&tsz, sizeof(int));
	for (int i=0; i<ResultV.size(); i++)
	{
		ResultV[i].write(ofile);
	}
	for (int i=0; i<ResultT.size(); i++)
	{
		const Tetrahedron* t = &ResultT[i];
		t->write(ofile);
	}

	ofile.close();
	return true;
}


bool cwg::TetrahedralMesh::save_result_v_t( const std::string& filename_verts, const std::string& filename_tets)
{
	//create_path_from_full_path(filename);
	ofstream ofile(filename_verts);
	if (!ofile)
		return false;
	
	int vsz = ResultV.size();
	int tsz = ResultT.size();

	ofile << "Verts. " << vsz << " \n";
	std::cout << "Verts. " << vsz << " \n"; 
	for (int i=0; i<ResultV.size(); i++)
	{
		ofile << ResultV[i].ID << "\n";
		ofile << ResultV[i].pos().x << "\t" << ResultV[i].pos().y << "\t" << ResultV[i].pos().z << "\n";
		ofile.flush();
	}
	ofile.flush();
	ofile.close();

	ofstream ofile_tet(filename_tets);
	if (!ofile_tet)
		return false;
	std::cout << "Saved all the verts!" << std::endl;
	int check_num;
	std::cin >> check_num;
	ofile_tet << "Tets. " << tsz << " \n";
	std::cout << "Tets. " << tsz << " \n";
	for (int i=0; i<ResultT.size(); i++)
	{
		if(ResultT[i].get_label() != 0){
			ofile_tet << ResultT[i].ID << "\n";
			ofile_tet << ResultT[i].get_label() << "\n";

			for(int j=0; j<4; j++){
				ofile_tet << ResultT[i].Verts[j] << "\t";
			}

			ofile_tet << "\n";
			ofile_tet.flush();
		}
	}

	std::cout << "Saved all the tets!" << std::endl;
	
	
	return true;
}


bool cwg::TetrahedralMesh::save_result_v_t_label( const std::string& filename_verts, const std::string& filename_tets, int label)
{
	//create_path_from_full_path(filename);
	ofstream ofile(filename_verts);
	if (!ofile)
		return false;
	
	std::map<int, int> save_vert_idx;

	int vsz = ResultV.size();
	int tsz = ResultT.size();

	ofile << "Verts. " << vsz << " \n";
	//std::cout << "Verts. " << vsz << " \n"; 
	int vert_idx = 0;
	for (int i=0; i<ResultV.size(); i++)
	{
		ofile << ResultV[i].ID << "\n";
		ofile << ResultV[i].pos().x << "\t" << ResultV[i].pos().y << "\t" << ResultV[i].pos().z << "\n";
		save_vert_idx.insert(make_pair(i, vert_idx));
		//std::cout << i << " ";
		vert_idx++;
		ofile.flush();
	}
	ofile.flush();
	ofile.close();

	ofstream ofile_tet(filename_tets);
	if (!ofile_tet)
		return false;
	//std::cout << "Saved all the verts!" << std::endl;
	//int check_num;
	//std::cin >> check_num;
	ofile_tet << "Tets. " << tsz << " \n";
	//std::cout << "Tets. " << tsz << " \n";
	for (int i=0; i<ResultT.size(); i++)
	{
		if(ResultT[i].get_label() == label){
			ofile_tet << ResultT[i].ID << "\n";
			ofile_tet << ResultT[i].get_label() << "\n";

			for(int j=0; j<4; j++){
				ofile_tet << ResultT[i].Verts[j] << "\t";
			}

			ofile_tet << "\n";
			ofile_tet.flush();
		}
	}

	//std::cout << "Saved all the tets!" << std::endl;
	
	
	return true;
}

bool cwg::TetrahedralMesh::load_bin_tetm( const std::string& filename )
{
	release();
	ifstream ifile;
	ifile.open(filename, ios::binary);
	if (ifile.is_open() == false)
		return false;

	int nverts, ntets;
	ifile.read( (char*)&nverts, sizeof(int) );
	ifile.read( (char*)&ntets, sizeof(int) );
	//m_storage_verts = new TetVertex [nverts];
	//m_storage_tets = new Tetrahedron [ntets];
	std::cout << "num of verts: " << nverts << std::endl;
	m_storage_verts.resize(nverts);
	m_storage_tets.resize(ntets);
	

	for (int i=0; i<nverts; i++)
	{
		TetVertex* v = &m_storage_verts[i];
		v->read(ifile);
		m_storage_verts[i] = *v;
	}

	for (int i=0; i<ntets; i++)
	{
		Tetrahedron* t = &m_storage_tets[i];
		t->read(ifile);
		t->correlate_verts(this);
		m_storage_tets[i] = *t;
	}
	ifile.close();
	return true;
}

bool cwg::TetrahedralMesh::load_bin_tetm_to_result( const std::string& filename )
{
	release();
	ifstream ifile;
	ifile.open(filename, ios::binary);
	if (ifile.is_open() == false)
		return false;

	int nverts, ntets;
	ifile.read( (char*)&nverts, sizeof(int) );
	ifile.read( (char*)&ntets, sizeof(int) );
	//m_storage_verts = new TetVertex [nverts];
	//m_storage_tets = new Tetrahedron [ntets];
	m_storage_verts.resize(nverts);
	m_storage_tets.resize(ntets);

	ResultV.resize(nverts);
	ResultT.resize(ntets);

	for (int i=0; i<nverts; i++)
	{
		TetVertex* v = &m_storage_verts[i];
		v->read(ifile);
		m_storage_verts[i] = *v;
		ResultV[i] = *v;
	}

	for (int i=0; i<ntets; i++)
	{
		Tetrahedron* t = &m_storage_tets[i];
		t->read(ifile);
		t->correlate_verts(this);
		m_storage_tets[i] = *t;
		ResultT[i] = *t;
	}

	ifile.close();
	return true;
}

bool cwg::TetrahedralMesh::save_ascii_cgal_from_result(const std::string& filename) {
		//create_path_from_full_path(filename);
	ofstream ofile(filename);
	if (!ofile){
		std::cout << "could not open file: " << filename << std::endl;
		return false;
	}
		
	int vsz = ResultV.size();
	int tsz = ResultT.size();
	
	ofile << "CGAL c3t3" << std::endl;
	ofile << "3" << std::endl;
	ofile << vsz << std::endl;

	std::cout << "Vert size: " << vsz << std::endl;
	std::cout << "Tet size: " << tsz << std::endl;

	for (int i=0; i<ResultV.size(); i++)
	{
		ofile << ResultV[i].pos().x << " " << ResultV[i].pos().y << " " << ResultV[i].pos().z << " -1" << " 0" << std::endl;
	}

	std::vector<int> border_tris;
	for (int i=0; i<Tris.size(); i++){
		int t0 = Tris[i]->Tets[0];
		int t1 = Tris[i]->Tets[1];

		if ( (t0 < 0 ) && (t1 < 0)){
			throw exception("bool cwg::TetrahedralMesh::save_txt_surf_ply( const std::string& filename ) const:\
							t0 < 0 && t1 < 0");
		}else if((t0 < 0 ) || (t1 < 0)){
			border_tris.push_back(i);
		}
	}

	std::vector<std::vector<int>> indexes_vector(tsz+border_tris.size());
	std::vector<int> tet_labels(tsz+border_tris.size());

	//push back finite cells
	for (int i=0; i<ResultT.size(); i++){
		std::vector<int> tet_indexes(4);
		for(int j=0; j<4; j++){
			tet_indexes[j] = ResultT[i].Verts[j];
		}
		indexes_vector[i].push_back(tet_indexes[0]+1);
		indexes_vector[i].push_back(tet_indexes[1]+1);

		if (check_valid(ResultV[tet_indexes[0]].pos(), ResultV[tet_indexes[1]].pos(), ResultV[tet_indexes[2]].pos(), ResultV[tet_indexes[3]].pos())){
			indexes_vector[i].push_back(tet_indexes[2]+1);
			indexes_vector[i].push_back(tet_indexes[3]+1);
		}else{
			indexes_vector[i].push_back(tet_indexes[3]+1);
			indexes_vector[i].push_back(tet_indexes[2]+1);
		}
		tet_labels[i] = ResultT[i].get_label()+1;
	}

	//push back infinite cells
	for(int i=0; i<border_tris.size(); i++){
		int t0 = Tris[border_tris[i]]->Tets[0];
		int t1 = Tris[border_tris[i]]->Tets[1];

		int t = t0 < 0 ? t1:t0;
		//get the neighbor of the infinite cell
		std::vector<int> valid_tet_verts(4);
		for (int j=0; j<4; j++){
			valid_tet_verts[j] = indexes_vector[t][j]-1;
		}

		//get the position of the other point in the neighbor of infinite cell
		int tet_vert_sum = 0;
		for (int j=0; j<3; j++){
			tet_vert_sum += std::find(valid_tet_verts.begin(), valid_tet_verts.end(), Tris[border_tris[i]]->Vert(j)) - valid_tet_verts.begin(); 
		}

		int infinite_vert_index = 6 - tet_vert_sum;

		//set infinite cell
		std::vector<int> tet_indexes(4);
		for (int j=0; j<4; j++){
			if(j != infinite_vert_index){
				tet_indexes[j] = valid_tet_verts[j]+1;
			}else{
				tet_indexes[j] = 0;
			}
		}

		indexes_vector[ResultT.size()+i].push_back(tet_indexes[0]);
		indexes_vector[ResultT.size()+i].push_back(tet_indexes[1]);
		indexes_vector[ResultT.size()+i].push_back(tet_indexes[3]);
		indexes_vector[ResultT.size()+i].push_back(tet_indexes[2]);

		tet_labels[ResultT.size()+i] = 0;

	}
	
	//write out tets
	ofile << (tsz+border_tris.size()) << std::endl;
	for (int i=0; i<(tsz+border_tris.size()); i++){
		ofile << indexes_vector[i][0] << " " << indexes_vector[i][1] << " " << indexes_vector[i][2] << " " << indexes_vector[i][3] << " " << std::endl;
	}

	//write out neighbours
	for(int i=0; i<(tsz+border_tris.size()); i++){
		std::vector<int> result = find_neighbour(indexes_vector, i);
		ofile << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
	}

	//write out labels
	for(int i=0; i<((tsz+border_tris.size())); i++){
		ofile << tet_labels[i] << " 0 0 0 0" << std::endl; 

	}

	ofile.close();
	return true;
}

bool cwg::TetrahedralMesh::save_ascii_cgal(const std::string& filename) {
		//create_path_from_full_path(filename);
	ofstream ofile(filename);
	if (!ofile){
		std::cout << "could not open file: " << filename << std::endl;
		return false;
	}
		

	int vsz = m_storage_verts.size();
	int tsz = m_storage_tets.size();
	
	
	ofile << "CGAL c3t3" << std::endl;
	ofile << "3" << std::endl;
	ofile << vsz << std::endl;

	std::cout << "Vert size: " << vsz << std::endl;
	std::cout << "Tet size: " << tsz << std::endl;

	for (int i=0; i<m_storage_verts.size(); i++)
	{
		ofile << m_storage_verts[i].pos().x << " " << m_storage_verts[i].pos().y << " " << m_storage_verts[i].pos().z << " -1" << " 0" << std::endl;
	}

	std::vector<int> border_tris;
	for (int i=0; i<Tris.size(); i++){
		int t0 = Tris[i]->Tets[0];
		int t1 = Tris[i]->Tets[1];

		if ( (t0 < 0 ) && (t1 < 0)){
			throw exception("bool cwg::TetrahedralMesh::save_txt_surf_ply( const std::string& filename ) const:\
							t0 < 0 && t1 < 0");
		}else if((t0 < 0 ) || (t1 < 0)){
			border_tris.push_back(i);
		}
	}

	std::vector<std::vector<int>> indexes_vector(tsz+border_tris.size());
	std::vector<int> tet_labels(tsz+border_tris.size());

	//push back finite cells
	for (int i=0; i<m_storage_tets.size(); i++){
		std::vector<int> tet_indexes(4);
		for(int j=0; j<4; j++){
			tet_indexes[j] = m_storage_tets[i].Verts[j];
		}
		indexes_vector[i].push_back(tet_indexes[0]+1);
		indexes_vector[i].push_back(tet_indexes[1]+1);

		if (check_valid(m_storage_verts[tet_indexes[0]].pos(), m_storage_verts[tet_indexes[1]].pos(), m_storage_verts[tet_indexes[2]].pos(), m_storage_verts[tet_indexes[3]].pos())){
			indexes_vector[i].push_back(tet_indexes[2]+1);
			indexes_vector[i].push_back(tet_indexes[3]+1);
		}else{
			indexes_vector[i].push_back(tet_indexes[3]+1);
			indexes_vector[i].push_back(tet_indexes[2]+1);
		}
		tet_labels[i] = m_storage_tets[i].get_label()+1;
	}

	//push back infinite cells
	for(int i=0; i<border_tris.size(); i++){
		int t0 = Tris[border_tris[i]]->Tets[0];
		int t1 = Tris[border_tris[i]]->Tets[1];

		int t = t0 < 0 ? t1:t0;
		//get the neighbor of the infinite cell
		std::vector<int> valid_tet_verts(4);
		for (int j=0; j<4; j++){
			valid_tet_verts[j] = indexes_vector[t][j]-1;
		}

		//get the position of the other point in the neighbor of infinite cell
		int tet_vert_sum = 0;
		for (int j=0; j<3; j++){
			tet_vert_sum += std::find(valid_tet_verts.begin(), valid_tet_verts.end(), Tris[border_tris[i]]->Vert(j)) - valid_tet_verts.begin(); 
		}

		int infinite_vert_index = 6 - tet_vert_sum;

		//set infinite cell
		std::vector<int> tet_indexes(4);
		for (int j=0; j<4; j++){
			if(j != infinite_vert_index){
				tet_indexes[j] = valid_tet_verts[j]+1;
			}else{
				tet_indexes[j] = 0;
			}
		}

		indexes_vector[m_storage_tets.size()+i].push_back(tet_indexes[0]);
		indexes_vector[m_storage_tets.size()+i].push_back(tet_indexes[1]);
		indexes_vector[m_storage_tets.size()+i].push_back(tet_indexes[3]);
		indexes_vector[m_storage_tets.size()+i].push_back(tet_indexes[2]);

		tet_labels[m_storage_tets.size()+i] = 0;

	}


	
	//write out tets
	ofile << (tsz+border_tris.size()) << std::endl;
	for (int i=0; i<(tsz+border_tris.size()); i++){
		ofile << indexes_vector[i][0] << " " << indexes_vector[i][1] << " " << indexes_vector[i][2] << " " << indexes_vector[i][3] << " " << std::endl;
	}

	//write out neighbours
	for(int i=0; i<(tsz+border_tris.size()); i++){
		std::vector<int> result = find_neighbour(indexes_vector, i);
		ofile << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
	}

	//write out labels
	for(int i=0; i<((tsz+border_tris.size())); i++){
		ofile << tet_labels[i] << " 0 0 0 0" << std::endl; 

	}

	ofile.close();
	return true;
}

bool cwg::TetrahedralMesh::save_mesh_from_result(const std::string& filename) {

	//create_path_from_full_path(filename);
	ofstream ofile(filename);
	if (!ofile){
		std::cout << "could not open file: " << filename << std::endl;
		return false;
	}
		
	int vsz = ResultV.size();
	int tsz = ResultT.size();
	
	ofile << "MeshVersionFormatted 1" << std::endl;
	ofile << "Dimension 3" << std::endl;
	

	std::cout << "Vert size: " << vsz << std::endl;
	std::cout << "Tet size: " << tsz << std::endl;

	ofile << std::endl;
	ofile << "Vertices" << std::endl;
	ofile << vsz << std::endl;

	for (int i=0; i< ResultV.size(); i++)
	{
		int label_this = 0;
		if (!ResultV[i].m_labels.empty()){
			label_this = ResultV[i].m_labels[0]+1;
		}
		ofile << ResultV[i].pos().x << " " << ResultV[i].pos().y << " " << ResultV[i].pos().z << "  " << label_this << std::endl;
	}

	std::vector<std::vector<int>> indexes_vector(tsz);
	std::vector<int> tet_labels(tsz);

	//push back finite cells
	for (int i=0; i<ResultT.size(); i++){
		std::vector<int> tet_indexes(4);
		for(int j=0; j<4; j++){
			tet_indexes[j] = ResultT[i].Verts[j];
		}
		indexes_vector[i].push_back(tet_indexes[0]+1);
		indexes_vector[i].push_back(tet_indexes[1]+1);

		if (check_valid(ResultV[tet_indexes[0]].pos(), ResultV[tet_indexes[1]].pos(), ResultV[tet_indexes[2]].pos(), ResultV[tet_indexes[3]].pos())){
			indexes_vector[i].push_back(tet_indexes[2]+1);
			indexes_vector[i].push_back(tet_indexes[3]+1);
		}else{
			indexes_vector[i].push_back(tet_indexes[3]+1);
			indexes_vector[i].push_back(tet_indexes[2]+1);
		}
		tet_labels[i] = ResultT[i].get_label()+1;
	}
	
	//write out tets
	ofile << std::endl;
	ofile << "Tetrahedra" << std::endl;
	ofile << tsz << std::endl;
	for (int i=0; i<(tsz); i++){
		ofile << indexes_vector[i][0] << " " << indexes_vector[i][1] << " " 
			<< indexes_vector[i][2] << " " << indexes_vector[i][3] << " " << tet_labels[i] << std::endl;
	}

	
	/*for (int i=0; i<ResultT.size(); i++){
		ofile << ResultT[i].Verts[0] + 1 << " " << ResultT[i].Verts[1] + 1 << " " 
			  << ResultT[i].Verts[2] + 1 << " " << ResultT[i].Verts[3] + 1 << " " << ResultT[i].get_label() << std::endl;
	}*/

	/*ofile << std::endl;
	ofile << "Triangles" << std::endl;
	ofile << tsz << std::endl;

	for (int i=0; i<Tris.size(); i++) {
		
	}*/
	ofile << "End" << std::endl;
	ofile.close();
	return true;	
}


bool cwg::TetrahedralMesh::save_mesh(const std::string& filename) {
	ofstream ofile(filename);
	if (!ofile){
		std::cout << "could not open file: " << filename << std::endl;
		return false;
	}
	
	int valid_verts=0, valid_tets=0;

	int vert_storage_sz = m_storage_verts.size();
	int tet_storage_sz = m_storage_tets.size();

	int valid_vert_num=0, valid_tet_num=0;
	for(int vert_index=0; vert_index < vert_storage_sz; vert_index++){
		if(m_storage_verts[vert_index].ID != -1){
			valid_vert_num++;
		}
	}

	for(int tet_index=0; tet_index < tet_storage_sz; tet_index++){
		if(m_storage_tets[tet_index].ID != -1){
			valid_tet_num++;
		}
	}
	
	ofile << "MeshVersionFormatted 1" << std::endl;
	ofile << "Dimension 3" << std::endl;
	
	std::cout << "Vert size: " << valid_vert_num << std::endl;
	std::cout << "Tet size: " << valid_tet_num << std::endl;

	ofile << std::endl;
	ofile << "Vertices" << std::endl;
	ofile << valid_vert_num << std::endl;

	std::map<int, int> vert_dict;

	int valid_vert_index=0;
	for (int i=0; i< m_storage_verts.size(); i++)
	{
		if(m_storage_verts[i].ID != -1){
			ofile << m_storage_verts[i].pos().x << " " << m_storage_verts[i].pos().y << " " << m_storage_verts[i].pos().z << "  " << m_storage_verts[i].m_labels[0]+1 << std::endl;
			vert_dict[i] = valid_vert_index;
			valid_vert_index ++;
		}
		
	}

	std::vector<std::vector<int>> indexes_vector(valid_tet_num);
	std::vector<int> tet_labels(valid_tet_num);

	//push back finite cells
	int valid_tet_index = 0;
	for (int i=0; i<m_storage_tets.size(); i++){
		if (m_storage_tets[i].ID != -1){
			std::vector<int> tet_indexes(4);
			for(int j=0; j<4; j++){
				tet_indexes[j] = m_storage_tets[i].Verts[j];
			}
			indexes_vector[valid_tet_index].push_back(vert_dict[tet_indexes[0]]+1);
			indexes_vector[valid_tet_index].push_back(vert_dict[tet_indexes[1]]+1);

			if (check_valid(m_storage_verts[tet_indexes[0]].pos(), m_storage_verts[tet_indexes[1]].pos(), 
							m_storage_verts[tet_indexes[2]].pos(), m_storage_verts[tet_indexes[3]].pos())){
				indexes_vector[valid_tet_index].push_back(vert_dict[tet_indexes[2]]+1);
				indexes_vector[valid_tet_index].push_back(vert_dict[tet_indexes[3]]+1);
			}else{
				indexes_vector[valid_tet_index].push_back(vert_dict[tet_indexes[3]]+1);
				indexes_vector[valid_tet_index].push_back(vert_dict[tet_indexes[2]]+1);
			}
			tet_labels[valid_tet_index] = m_storage_tets[i].get_label()+1;

			valid_tet_index++;
		}
		
	}
	
	//write out tets
	//ofile << std::endl;
	ofile << "Tetrahedra" << std::endl;
	ofile << valid_tet_index << std::endl;
	for (int i=0; i<valid_tet_index; i++){
		ofile << indexes_vector[i][0] << " " << indexes_vector[i][1] << " " 
			<< indexes_vector[i][2] << " " << indexes_vector[i][3] << " " << tet_labels[i] << std::endl;
	}

	
	/*for (int i=0; i<ResultT.size(); i++){
		ofile << ResultT[i].Verts[0] + 1 << " " << ResultT[i].Verts[1] + 1 << " " 
			  << ResultT[i].Verts[2] + 1 << " " << ResultT[i].Verts[3] + 1 << " " << ResultT[i].get_label() << std::endl;
	}*/

	/*ofile << std::endl;
	ofile << "Triangles" << std::endl;
	ofile << tsz << std::endl;

	for (int i=0; i<Tris.size(); i++) {
		
	}*/
	ofile << "End" << std::endl;
	ofile.close();
	return true;	
}

void cwg::write_ply_file_head(std::ofstream& ofile, int vertex_count, int face_count)
{
	ofile << "ply" << endl;
	ofile << "format ascii 1.0" << endl;
	ofile << "element vertex " << vertex_count << endl;
	ofile << "property float64 x " << endl;
	ofile << "property float64 y " << endl;
	ofile << "property float64 z " << endl;
	ofile << "property uint8 red " << endl;
	ofile << "property uint8 green " << endl;
	ofile << "property uint8 blue " << endl;
	ofile << "element face " << face_count << endl;
	ofile << "property list uint8 int32 vertex_index" << endl;
	ofile << "property uint8 red " << endl;
	ofile << "property uint8 green " << endl;
	ofile << "property uint8 blue " << endl;
	ofile << "end_header" << endl;
}
