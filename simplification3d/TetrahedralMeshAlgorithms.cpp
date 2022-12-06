#include "stdafx.h"
#include "TetrahedralMesh.h"
#include "Utils/MemorySurveillance.h"

#include "H5Cpp.h"
using namespace H5;

using namespace cv;
using namespace std;
using namespace Cleaver;


SIMPLIFICATION3D_API bool glb_boundary_fixed = false;
SIMPLIFICATION3D_API float glb_theta = 0.8;
//Tet adjust

const vector<string> explode(const string& s, const char& c)
{
	string buff("");
	vector<string> v;

	for (auto n : s)
	{
		if (n != c) buff += n; else
			if (n == c && buff != "") { v.push_back(buff); buff = ""; }
	}
	if (buff != "") v.push_back(buff);

	return v;
}


cwg::TetrahedralMesh::TetrahedralMesh()
{
	m_cleaver_tet_mesh = NULL;
	AdjMatrix = NULL;
	m_labelvolume = NULL;
	m_tar_verts = -1;
	m_final_percentage = 1.0;
	m_X = m_Y = m_Z = -1;

	InternalMesh = new TriangularMesh();

	process_all = false;
	heap_time = 0;
	keep_in_previous = false;
	result_v_id = result_t_id = 0;
	render_mode = false;

	for(int i=0; i< 3; i++){
		walls_begin[i] = 0;
		walls_end[i] = 0;
	}

}


cwg::TetrahedralMesh::~TetrahedralMesh()
{
	SAFE_DELETE(m_cleaver_tet_mesh);
	release();
}


void cwg::TetrahedralMesh::init_memory(int slidex, int slidey, int slidez)
{
	// slides storage
	m_storage_verts.resize((slidex + 1)*(slidey + 1)*(slidez + 1));
	m_storage_tets.resize( slidex *slidey * slidez * tetn);		//=======5 -- number of tet===============

	m_storage_edges.resize(m_storage_tets.size() * 10);	//roughly total number of edges

	m_storage_tets_tmp.resize( m_storage_tets.size());

	vector<Tetrahedron>::iterator& it = m_storage_tets_tmp.begin();
	int size = m_storage_tets.size();

	for(int tid_tmp = 0; tid_tmp < size; tid_tmp++){

		it = m_storage_tets_tmp.begin() + tid_tmp;
		it->ID = -1;
		for(int i=0; i<4; i++){
			it->Verts[i] = -1;
		}
	}

	AdjMatrix =  new Sparse_Matrix(m_storage_verts.size(), m_storage_verts.size(), SYM_BOTH, false, CRS);
	AdjMatrix->begin_fill_entry();

	to_debug_bool = false;

}


void cwg::TetrahedralMesh::clear_memory()
{
	m_storage_verts.clear();
	m_storage_verts.shrink_to_fit();
	m_storage_tets.clear();
	m_storage_tets.shrink_to_fit();
	m_storage_edges.clear();
	m_storage_edges.shrink_to_fit();
	m_storage_tets_tmp.clear();
	m_storage_tets_tmp.shrink_to_fit();
	delete AdjMatrix;

	SAFE_DELETE(m_cleaver_tet_mesh);
	EdgeSet.clear();
	EdgeSet.empty();
}

//=====================wenhua============================ 

//*****************************this one is tested by wenhua, it seems right**************
void cwg::TetrahedralMesh::tetrahedronlize(int i, int j, int k, TetCube& tmp_cube, vector<Tetrahedron>::iterator& it, int& tID){
	if (tetn == 6)
	{
		tmp_cube.ComputeTetrahedra(it, tID);
	}
	else{
	//==============
		if((k&1) == 0 && (i&1) == 0 && (j&1) ==0)
		{
			tmp_cube.ComputeTetrahedra2(it, tID);
		}
		else if((k&1) == 0 && (i&1) == 1 && (j&1) ==0)
		{
			tmp_cube.ComputeTetrahedra3(it, tID);
		}
		else if((k&1) == 0 && (i&1) == 0 && (j&1) ==1)
		{
			tmp_cube.ComputeTetrahedra3(it, tID);
		}
		else if((k&1) == 0 && (i&1) == 1 && (j&1) ==1)
		{
			tmp_cube.ComputeTetrahedra2(it, tID);
		}
		else if((k&1) == 1 && (i&1) == 0 && (j&1) ==0)
		{
			tmp_cube.ComputeTetrahedra3(it, tID);
		}
		else if((k&1) == 1 && (i&1) == 1 && (j&1) ==0)
		{
			tmp_cube.ComputeTetrahedra2(it, tID);
		}
		else if((k&1) == 1 && (i&1) == 0 && (j&1) ==1)
		{
			tmp_cube.ComputeTetrahedra2(it, tID);
		}
		else if((k&1) == 1 && (i&1) == 1 && (j&1) ==1)
		{
			tmp_cube.ComputeTetrahedra3(it, tID);
		}
	}
}


//reset the save_vidx and the save_tidx
void cwg::TetrahedralMesh::walls(Visualize& vis, int startx, int starty, int startz, int endx, int endy, int endz){

	if(startx == 0)
		walls_begin[0] = 1;
	else
		walls_begin[0] = 0;

	if(starty == 0)
		walls_begin[1] = 1;
	else
		walls_begin[1] = 0;

	if(startz == 0)
		walls_begin[2] = 1;
	else
		walls_begin[2] = 0;

	if(endx == m_X)
		walls_end[0] = 1;
	else
		walls_end[0] = 0;

	if(endy == m_Y)
		walls_end[1] = 1;
	else
		walls_end[1] = 0;

	if(endz == m_Z)
		walls_end[2] = 1;
	else
		walls_end[2] = 0;

	std::cout << "starts: " << startx << "\t" << starty << "\t" << startz << std::endl;
	std::cout << "ends: " << endx << "\t" << endy << "\t" << endz << std::endl;

}


//set the slideboundaries for all the windows
void cwg::TetrahedralMesh::sliding_boundaries(Visualize& vis, int startx, int starty, int startz, int endx, int endy, int endz){


	int idx = 0;

	bool end_bool = (endx == m_X) || (endy == m_Y) || (endz == m_Z);

	//the boundary on the face of startz
		//k=startz and k=startz+1
	int slide_x_length = endx - startx;
	int slide_y_length = endy - starty;
	int slide_z_length = endz - startz;
		
	//set for z
	if(!walls_end[2]){
		for(int j=0; j<=slide_x_length; j++){
			for(int i=0; i<=slide_y_length; i++){
				idx = (slide_z_length-1) * (slide_X+1) * (slide_Y+1) + j*(slide_Y+1) + i;
				m_storage_verts[idx].boundary = true;

				idx = (slide_z_length) *(slide_X+1)*(slide_Y+1) + j*(slide_Y+1) + i;
				m_storage_verts[idx].boundary = true;
			}
		}
	}

	//set for x
	if( !walls_end[0]){
		for(int k=0; k<=slide_z_length; k++){
			for(int i=0; i<=slide_y_length; i++){
				idx = k*(slide_Y+1)*(slide_X+1) + (slide_x_length-1)*(slide_Y+1) + i;
				m_storage_verts[idx].slideBoundary = true;

				idx = k*(slide_X+1)*(slide_Y+1) + (slide_x_length)*(slide_Y+1) + i;
				m_storage_verts[idx].slideBoundary = true;
			}
		}
	}

	//y
	if(!walls_end[1]){
		for(int k=0; k<=slide_z_length; k++){
			for(int j=0; j<=slide_x_length; j++){

				idx = k*(slide_Y+1)*(slide_X+1) + j*(slide_Y+1) + (slide_y_length-1);
				m_storage_verts[idx].boundary = true;

				idx = k*(slide_Y+1)*(slide_X+1) + j*(slide_Y+1) + (slide_y_length);
				m_storage_verts[idx].boundary = true;
			}
		}
	}
}


void cwg::TetrahedralMesh::reset_win_verts(int win_inside_index, int startx, int starty, int startz){
	int win_startx = (slide_X/2) * ( (win_inside_index % 4) / 2);
	int win_starty = (slide_Y/2) * ( win_inside_index % 2);
	int win_startz = (slide_Z/2) * ( win_inside_index / 4);
	int win_endx   = win_startx + slide_X/2;
	int win_endy   = win_starty + slide_Y/2;
	int win_endz   = win_startz + slide_Z/2;

	int vert_save_startx, vert_save_starty, vert_save_startz;

	if( ((win_inside_index % 4) / 2) == 0 )
		vert_save_startx = 0;
	else
		vert_save_startx = 1;

	if( ( win_inside_index % 2) == 0)
		vert_save_starty = 0;
	else
		vert_save_starty = 1;

	if( ( win_inside_index / 4) == 0)
		vert_save_startz = 0;
	else
		vert_save_startz = 1;

	for(int k=win_startz+vert_save_startz; k<=win_endz; k++){
		for(int i=win_startx+vert_save_startx; i<=win_endx; i++){
			for(int j=win_starty+vert_save_starty; j<=win_endy; j++){

				int idx = (k)*(slide_X+1)*(slide_Y+1) + (i)*(slide_Y+1) + (j);

				m_storage_verts[idx].pos() = Cleaver::vec3(startx+i,starty+j,startz+k);
				m_storage_verts[idx].update_funcval();
				m_storage_verts[idx].invalidate();//vert.tets.clear(); //vert.m_labels.clear();
				m_storage_verts[idx].slideBoundary = false;
				m_storage_verts[idx].boundary = false;
				m_storage_verts[idx].keep = 0;
				m_storage_verts[idx].ID = -1;
				m_storage_verts[idx].OrigBoundaryTris.clear();
				m_storage_verts[idx].OuterBoundaryTris.clear();
				m_storage_verts[idx].weighted = false;
				m_storage_verts[idx].weightLevel = -1;

			}
		}
	}

	for(int k=win_startz; k<win_endz; k++){
		for(int i=win_startx; i<win_endx; i++){
			for(int j=win_starty; j<win_endy; j++){
				int idx = ((k)*(slide_X)*(slide_Y) + (i)*(slide_Y) + (j))*tetn;

				for(int index0=0; index0<tetn; index0++){
					m_storage_tets[idx+index0].ID = -1;
					for(int vert_index=0; vert_index<4; vert_index++){
						m_storage_tets[idx+index0].Verts[vert_index] = -1;
						m_storage_tets[idx+index0].tmp_Verts[vert_index][0] = -1;
						m_storage_tets[idx+index0].tmp_Verts[vert_index][1] = -1;
						m_storage_tets[idx+index0].tmp_Verts[vert_index][2] = -1;
					}
				}
			}
		}
	}

	//reset m_storage_tets_tmp
	for(int k=win_startz; k<win_endz; k++){
		for(int i=win_startx; i<win_endx; i++){
			for(int j=win_starty; j<win_endy; j++){
				int idx = (k*slide_X*slide_Y + i *slide_Y + j)*tetn;

				for(int index0=0; index0<tetn; index0++){
					m_storage_tets_tmp[idx+index0].ID = -1;
					for(int vert_index=0; vert_index<4; vert_index++){
						m_storage_tets_tmp[idx+index0].Verts[vert_index] = -1;
						m_storage_tets_tmp[idx+index0].tmp_Verts[vert_index][0] = -1;
						m_storage_tets_tmp[idx+index0].tmp_Verts[vert_index][1] = -1;
						m_storage_tets_tmp[idx+index0].tmp_Verts[vert_index][2] = -1;
					}
				}
			}
		}
	}
}


void cwg::TetrahedralMesh::reset_all(){
	int size = m_storage_tets.size();
	
	//reset the m_storage_tet_tmp for future use
	for(int i=0; i<size; i++){
		m_storage_tets_tmp[i].ID = -1;
		for(int index0=0; index0<4; index0++){
			m_storage_tets_tmp[i].Verts[index0] = -1;
			m_storage_tets_tmp[i].tmp_Verts[index0][0] = -1;
			m_storage_tets_tmp[i].tmp_Verts[index0][1] = -1;
			m_storage_tets_tmp[i].tmp_Verts[index0][2] = -1;
		}
	}

	for(int i=0; i<size; i++){
		m_storage_tets[i].ID = -1;
		for(int index0=0; index0<4; index0++){
			m_storage_tets[i].Verts[index0] = -1;
			m_storage_tets[i].tmp_Verts[index0][0] = -1;
			m_storage_tets[i].tmp_Verts[index0][1] = -1;
			m_storage_tets[i].tmp_Verts[index0][2] = -1;
		}
	}

	size = m_storage_verts.size();
	for(int i=0; i<size; i++){
			m_storage_verts[i].pos() = Cleaver::vec3(0,0,0);
			m_storage_verts[i].update_funcval();
			m_storage_verts[i].invalidate();//vert.tets.clear(); //vert.m_labels.clear();
			m_storage_verts[i].slideBoundary = false;
			m_storage_verts[i].boundary = false;
			m_storage_verts[i].keep = 0;
			m_storage_verts[i].ID = -1;
			m_storage_verts[i].OrigBoundaryTris.clear();
			m_storage_verts[i].OuterBoundaryTris.clear();
			m_storage_verts[i].weighted = false;
			m_storage_verts[i].weightLevel = -1;
	}
}


void cwg::TetrahedralMesh::tet_Adjust(){
	int size = m_storage_tets.size();

	std::vector<Tetrahedron> m_storage_tets_tmp_adjust;
	m_storage_tets_tmp_adjust.resize(size);

	vector<Tetrahedron>::iterator& it_adjust = m_storage_tets_tmp_adjust.begin();
	for(int tid_tmp = 0; tid_tmp < size; tid_tmp++){

		it_adjust = m_storage_tets_tmp_adjust.begin() + tid_tmp;
		it_adjust->ID = -1;
		for(int i=0; i<4; i++){
			it_adjust->Verts[i] = -1;
		}
	}

	std::vector<int> replace_array;
	replace_array.resize(size);
	for(int i=0; i<size; i++){
		replace_array[i] = -1;
	}

	for(int k=0; k<2; k++){
		for(int i=0; i<2; i++){
			for(int j=0; j<2; j++){
				int cube_kind_index = k * 2 * 2 + i * 2 + j;
				start_index[cube_kind_index][0] = (i) * slide_X/2; //x
				start_index[cube_kind_index][1] = (j) * slide_Y/2; //y
				start_index[cube_kind_index][2] = (k) * slide_Z/2; //z
				start_index[cube_kind_index][3] = 0;               //index

			}
		}
	}

	for(int ti=0; ti<size; ti++){
		if(m_storage_tets[ti].ID != -1){
			int vert_xyz[4][3];

			for(int vert_index=0; vert_index<4; vert_index++){
				vert_xyz[vert_index][0] = ( (m_storage_tets[ti].Verts[vert_index] % ((slide_X+1)*(slide_Y+1))) / (slide_Y+1) );
				vert_xyz[vert_index][1] = ( (m_storage_tets[ti].Verts[vert_index] % ((slide_X+1)*(slide_Y+1))) % (slide_Y+1) );
				vert_xyz[vert_index][2] = (  m_storage_tets[ti].Verts[vert_index] / ((slide_X+1)*(slide_Y+1)) );
			}

			replace_array[ti] = find_adjust_tet_id(vert_xyz);
		}
	}


	for(int i=0; i<size; i++){
		if(replace_array[i] != -1){
			int new_i = replace_array[i];
			m_storage_tets_tmp_adjust[new_i].ID = new_i;
			m_storage_tets_tmp_adjust[new_i].set_label(m_storage_tets[i].get_label());

			for(int index0=0; index0<4; index0++){
				m_storage_tets_tmp_adjust[new_i].Verts[index0] = m_storage_tets[i].Verts[index0];
			}

		}
	}

	for(int ti=0; ti<size; ti++){
		if(m_storage_tets_tmp_adjust[ti].ID != -1){
			m_storage_tets[ti].ID = ti;
			m_storage_tets[ti].set_label(m_storage_tets_tmp_adjust[ti].get_label());

			for(int index0=0; index0<4; index0++){
				m_storage_tets[ti].Verts[index0] = m_storage_tets_tmp_adjust[ti].Verts[index0];
			}
		}else{
			m_storage_tets[ti].ID = -1;
			for(int index0=0; index0<4; index0++){
				m_storage_tets[ti].Verts[index0] = -1;
			}
		}
	}
			
	m_storage_tets_tmp_adjust.clear();
	m_storage_tets_tmp_adjust.shrink_to_fit();

}

//To adjust the keep number of the verts
//If the verts are in the left part, which means they will be stored in the text
//In this case, if they are not in the same cube, they will be seperated.
//So for the tets that are to be stored and not in the same cube, it will keep++;
void cwg::TetrahedralMesh::keep_Adjust(bool x_end){

	int size = m_storage_tets.size();

	for(int i=0; i<size; i++){		
		if(m_storage_tets[i].ID != -1){
			int cube_index[4];

			//adjusted it to make the middle line right
			for(int index0=0; index0<4; index0++){
			
				int z_left_right =( m_storage_tets[i].Verts[index0] / ((slide_X+1)*(slide_Y+1)) );
				int y_left_right =( (m_storage_tets[i].Verts[index0] % ((slide_X+1)*(slide_Y+1))) % (slide_Y+1) );
				int x_left_right =( (m_storage_tets[i].Verts[index0] % ((slide_X+1)*(slide_Y+1))) / (slide_Y+1) );
				
				int cube_indexes[3];

				cube_indexes[0] = (x_left_right - 1)/(slide_X/2);
				cube_indexes[0] = cube_indexes[0] < 0 ? 0 : cube_indexes[0];

				cube_indexes[1] = (y_left_right - 1)/(slide_Y/2);
				cube_indexes[1] = cube_indexes[1] < 0 ? 0 : cube_indexes[1];

				cube_indexes[2] = (z_left_right - 1)/(slide_Z/2);
				cube_indexes[2] = cube_indexes[2] < 0 ? 0 : cube_indexes[2];

				cube_index[index0] = cube_indexes[2] * 2 * 2 + cube_indexes[0]*2 + cube_indexes[1];
			}
	
			

			if( !x_end){
				bool all_left = true;
				
				for(int vert_i=0; vert_i <4; vert_i++){
					int z_left_right =( m_storage_tets[i].Verts[vert_i] / ((slide_X+1)*(slide_Y+1)) );
					int y_left_right =( (m_storage_tets[i].Verts[vert_i] % ((slide_X+1)*(slide_Y+1))) % (slide_Y+1) );
					int x_left_right =( (m_storage_tets[i].Verts[vert_i] % ((slide_X+1)*(slide_Y+1))) / (slide_Y+1) );

					//all the verts inside the right half the of processing cube
					all_left = all_left && ( z_left_right >= 0  && z_left_right <= slide_Z); 
					all_left = all_left && ( y_left_right >= 0  && y_left_right <= slide_Y); 
					all_left = all_left && ( x_left_right > (slide_X/2) && x_left_right <= slide_X);
				}

				if(!all_left){

					bool all_in_one_bool = (cube_index[0] == cube_index[1]) && (cube_index[1] == cube_index[2]) && (cube_index[2] == cube_index[3]);

					if(!all_in_one_bool){
						for(int index0=0; index0<4; index0++)
							m_storage_verts[m_storage_tets[i].Verts[index0]].keep++;
					}
				}
			}else{
					bool all_in_one_bool = (cube_index[0] == cube_index[1]) && (cube_index[1] == cube_index[2]) && (cube_index[2] == cube_index[3]);

					if(!all_in_one_bool){
						for(int index0=0; index0<4; index0++)
							m_storage_verts[m_storage_tets[i].Verts[index0]].keep++;

					}
			}


		}
	}

}

//the function that would replace loadVolumeDate_sliding
void cwg::TetrahedralMesh::loadVolumeData_sliding_wenhua(Visualize& vis, const std::string& foldername , const std::vector<std::string>& filelist, int startx, int starty, int startz, int endx, int endy, int endz){
	//for test purpose

	walls(vis, startx, starty, startz, endx, endy, endz);

	int win_x_start = startx / (slide_X/2);
	int win_y_start = starty / (slide_Y/2);
	int win_z_start = startz / (slide_Z/2);

	
	if(walls_begin[2]){

		if(walls_begin[1]){
			if(walls_begin[0]){	//situation 0  0, 0, 0
				//load from file 0,1,2,3,4,5,6,7

				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 0);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 1);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 2);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 3);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 4);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 5);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 6);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 7);
			}
			else{                  //situation 1  0, 0, 1
				//load from right (0,1,4,5) from (2,3,6,7)
				loadData_fromRight_wenhua(startx, starty, startz);

				//load_from file 2,3,6,7
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 2);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 3);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 6);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 7);				
			}

		}else{

			if(walls_begin[0]){    //situation 2  0, 1, 0

				reset_all();
				//load from file 1,3,5,7
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 1);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 3);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 5);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 7);

				//load from text 0,2,4,6
				loadData_fromText_wenhua(win_x_start  , win_y_start  , win_z_start  , 0);
				loadData_fromText_wenhua(win_x_start+1, win_y_start  , win_z_start  , 2);
				loadData_fromText_wenhua(win_x_start  , win_y_start  , win_z_start+1, 4);
				loadData_fromText_wenhua(win_x_start+1, win_y_start  , win_z_start+1, 6);

			}else{                 //situation 3 0, 1, 1
				//load from right (0,1,4,5) from (2,3,6,7)
				loadData_fromRight_wenhua(startx, starty, startz);

				//load from file 3,7
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 3);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 7);

				//load from text 2,6
				loadData_fromText_wenhua(win_x_start+1, win_y_start  , win_z_start  , 2);
				loadData_fromText_wenhua(win_x_start+1, win_y_start  , win_z_start+1, 6);
				
			}
		}

	}else{

		if(walls_begin[1]){
			if(walls_begin[0]){   //situation 4  1, 0, 0
				//load from file 4,5,6,7
				reset_all();
				//load from text 0,1,2,3
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 4);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 5);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 6);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 7);

				loadData_fromText_wenhua(win_x_start  , win_y_start  , win_z_start  , 0);
				loadData_fromText_wenhua(win_x_start  , win_y_start+1, win_z_start  , 1);
				loadData_fromText_wenhua(win_x_start+1, win_y_start  , win_z_start  , 2);
				loadData_fromText_wenhua(win_x_start+1, win_y_start+1, win_z_start  , 3);

			}else{                //situation 5  1, 0, 1
				//load from right (0,1,4,5) from (2,3,6,7)
				loadData_fromRight_wenhua(startx, starty, startz);

				//load from file 2,3
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 6);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 7);

				//load from text 6,7
				loadData_fromText_wenhua(win_x_start+1, win_y_start  , win_z_start, 2);
				loadData_fromText_wenhua(win_x_start+1, win_y_start+1, win_z_start, 3);

			}
		}else{
			if(walls_begin[0]){  //situation 6  1, 1, 0
				//load from file 5,7
				reset_all();
				loadData_fromFile_wenhua(foldername, filelist, startx, starty ,startz, 5);
				loadData_fromFile_wenhua(foldername, filelist, startx, starty, startz, 7);

				//load from text 0,1,2,3,4,6
				loadData_fromText_wenhua(win_x_start  , win_y_start  , win_z_start  , 0);
				loadData_fromText_wenhua(win_x_start  , win_y_start+1, win_z_start  , 1);
				loadData_fromText_wenhua(win_x_start+1, win_y_start  , win_z_start  , 2);
				loadData_fromText_wenhua(win_x_start+1, win_y_start+1, win_z_start  , 3);
				loadData_fromText_wenhua(win_x_start  , win_y_start  , win_z_start+1, 4);
				loadData_fromText_wenhua(win_x_start+1, win_y_start  , win_z_start+1, 6);

			}else{               //situation 7  1, 1, 1
				//load from right (0,1,4,5) from (2,3,6,7)
				loadData_fromRight_wenhua(startx, starty, startz);

				//load from file 7
				loadData_fromFile_wenhua(foldername, filelist, startx, starty ,startz, 7);
				//load from text 2,3,6
				loadData_fromText_wenhua(win_x_start+1, win_y_start  , win_z_start  , 2);
				loadData_fromText_wenhua(win_x_start+1, win_y_start+1, win_z_start  , 3);
				loadData_fromText_wenhua(win_x_start+1, win_y_start  , win_z_start+1, 6);

			}
		}
	}

	

	loadData_betweenCubes_wenhua();

	sliding_boundaries(vis, startx, starty, startz, endx, endy, endz);
	
	int size = m_storage_tets.size();
	for(int i=0; i<size; i++){
		if(m_storage_tets[i].ID != -1)
			for (int a=0; a<4; a++)			
				m_storage_verts[ m_storage_tets[i].Verts[a] ].correlate_tet(i);
	}

	compute_verts_border_code();
}

//Vert file format
//==================================================================================
//ID.x ID.y ID.z
//keep_number
//x y z
//label0	label1

//Tet file format
//=================================================================================
//tet_idx tet_idy tet_idz index0
//0/1 (whether the verts are inside the cube, 0 means no, 1 means yes)
//label
//1	vert0.idx vert0.idy vert0.idz (the first 0/1 number shows whether this vert is inside the cube)
//1	vert1.idx vert1.idy vert1.idz
//1	vert2.idx vert2.idy vert2.idz
//1	vert3.idx vert3.idy vert3.idz
//=

//This cube means the big cube, which is the whole m_storage_vert and the m_storage_tet

void cwg::TetrahedralMesh::saveResultData_wenhua(int data_buffer_startx, int data_buffer_starty, int data_buffer_startz, int buffer_win_x, int buffer_win_y, int buffer_win_z){	
	
		
	int win_x = data_buffer_startx/(slide_X/2) + ( buffer_win_x);
	int win_y = data_buffer_starty/(slide_Y/2) + ( buffer_win_y);
	int win_z = data_buffer_startz/(slide_Z/2) + ( buffer_win_z);

	int step_here[3] = {slide_X/2, slide_Y/2, slide_Z/2};
	if(walls_end[0] && buffer_win_x==1)
		step_here[0] =  steps[0];

	if(walls_end[1] && buffer_win_y==1)
		step_here[1] = steps[1];

	if(walls_end[2] && buffer_win_z==1)
		step_here[2] = steps[2];	
	
	int buffer_win_startx = (slide_X/2) * ( buffer_win_x);
	int buffer_win_starty = (slide_Y/2) * ( buffer_win_y);
	int buffer_win_startz = (slide_Z/2) * ( buffer_win_z);

	int buffer_win_endx   = buffer_win_startx + step_here[0];
	int buffer_win_endy   = buffer_win_starty + step_here[1];
	int buffer_win_endz   = buffer_win_startz + step_here[2];

	int vert_save_startx, vert_save_starty, vert_save_startz;
	
	if(buffer_win_x == 0)
		vert_save_startx = 0;
	else
		vert_save_startx = 1;

	if(buffer_win_y == 0)
		vert_save_starty = 0;
	else
		vert_save_starty = 1;

	if(buffer_win_z == 0)
		vert_save_startz = 0;
	else
		vert_save_startz = 1;

	//this will be shared by both vert file and the tet file	
	std::string tmpfilename("");
	tmpfilename += std::to_string(win_x);
	tmpfilename += "_";
	tmpfilename += std::to_string(win_y);
	tmpfilename += "_";
	tmpfilename += std::to_string(win_z);
	tmpfilename += ".txt";

	//Vert file
	std::ofstream vertFile;
	std::string vertFileName = get_foldername() + "Data\\vert_";
	vertFileName += tmpfilename;

	vertFile.open(vertFileName, ios::binary);
	int vert_count = 0;
	int tet_size = m_storage_tets.size();
	for(int k=vert_save_startz; k<=slide_Z/2; k++){
		for(int i=vert_save_startx; i<=slide_X/2; i++){
			for(int j=vert_save_starty; j<=slide_Y/2; j++){
				int vert_idx = (k+buffer_win_startz)*(slide_X+1)*(slide_Y+1) + (i+buffer_win_startx)*(slide_Y+1) + (j+buffer_win_starty);

				if(m_storage_verts[vert_idx].ID != -1){
					vertFile.write((char*)&i, sizeof(int));
					vertFile.write((char*)&j, sizeof(int));
					vertFile.write((char*)&k, sizeof(int));

					int keep_num = m_storage_verts[vert_idx].keep;
					vertFile.write((char*)&keep_num, sizeof(int));

					double vert_pos_x = m_storage_verts[vert_idx].pos().x;
					double vert_pos_y = m_storage_verts[vert_idx].pos().y;
					double vert_pos_z = m_storage_verts[vert_idx].pos().z;
					
					vertFile.write((char*)&vert_pos_x, sizeof(double));
					vertFile.write((char*)&vert_pos_y, sizeof(double));
					vertFile.write((char*)&vert_pos_z, sizeof(double));

					int label_size = m_storage_verts[vert_idx].m_labels.size();
					vertFile.write((char*)&label_size, sizeof(int));

					for(int index0=0; index0<label_size; index0++){
						int m_label_index0 = m_storage_verts[vert_idx].m_labels[index0];
						vertFile.write((char*) &m_label_index0, sizeof(int));
					}

					vert_count ++;
				}	
			}
		}
	}
	vertFile.write((char*) &vert_count, sizeof(int));
	vertFile.close();

	//Tet tmp file
	std::ofstream tetFile;
	std::string tetFilename  = get_foldername() + "Data\\tet_";
	
	tetFilename += tmpfilename;
	tetFile.open(tetFilename, ios::binary);
	int tet_count = 0;

	for(int k=0; k<slide_Z/2; k++){
		for(int i=0; i<slide_X/2; i++){
			for(int j=0; j<slide_Y/2; j++){
				int tet_idx0 = ((k+buffer_win_startz)*(slide_X)*(slide_Y) + (i+buffer_win_startx)*(slide_Y) + (j+buffer_win_starty))*tetn;
	
				for(int index=0; index<tetn; index++){
					int tet_idx = tet_idx0 + index;

					if(m_storage_tets_tmp[tet_idx].ID != -1){
						tetFile.write((char*) &i, sizeof(int));
						tetFile.write((char*) &j, sizeof(int));
						tetFile.write((char*) &k, sizeof(int));
						tetFile.write((char*) &index, sizeof(int));
						
						int tet_tmp_label = m_storage_tets_tmp[tet_idx].get_label();
						int tmp_num = 0;
						tetFile.write((char*) &tmp_num, sizeof(int));
						tetFile.write((char*) &tet_tmp_label, sizeof(int));

						bool inside_cubes_test[4];
						for(int index0=0; index0<4; index0++){

							int vert_idx = m_storage_tets_tmp[tet_idx].tmp_Verts[index0][0];
							int vert_idy = m_storage_tets_tmp[tet_idx].tmp_Verts[index0][1];
							int vert_idz = m_storage_tets_tmp[tet_idx].tmp_Verts[index0][2];
							
							bool inside_cube = true;
							inside_cube = (vert_idx >= (buffer_win_startx+vert_save_startx) && (vert_idx <= buffer_win_endx));
							inside_cube = inside_cube && (vert_idy >= (buffer_win_starty+vert_save_starty)) && (vert_idy <= buffer_win_endy);
							inside_cube = inside_cube && (vert_idz >= (buffer_win_startz+vert_save_startz)) && (vert_idz <= buffer_win_endz);
							
							if(inside_cube){
								//int vert_idx = (k+buffer_win_startz)*(slide_X+1)*(slide_Y+1) + (i+buffer_win_startx)*(slide_Y+1) + (j+buffer_win_starty);
								int tmp_vert_idx = vert_idz * (slide_X+1)*(slide_Y+1) + vert_idx*(slide_Y+1) + vert_idy;
								inside_cube = (m_storage_verts[tmp_vert_idx].ID != -1);
							}

							inside_cubes_test[index0] = inside_cube;
							
							tetFile.write((char*)&inside_cube, sizeof(bool));
							
							int vert_ids[3];
							vert_ids[0] = vert_idx-buffer_win_startx;
							vert_ids[1] = vert_idy-buffer_win_starty;
							vert_ids[2] = vert_idz-buffer_win_startz;
							tetFile.write((char*) vert_ids, sizeof(int)*3);
						}

						/*bool not_inside_cube_cube = !inside_cubes_test[0] && !inside_cubes_test[1] && !inside_cubes_test[2] && !inside_cubes_test[3];
						if(not_inside_cube_cube){
							std::cout << "Are you kidding? You have done tet adjust!" << std::endl;
							std::cout << "buffer win: " << buffer_win_x << "\t" << buffer_win_y << "\t" << buffer_win_z << std::endl;

							for(int index0=0; index0<4; index0++){
								std::cout << "vert_index: ";
								int cube_index = (m_storage_tets_tmp[tet_idx].tmp_Verts[index0][0]);
								std::cout << cube_index << "\t";

								cube_index = (m_storage_tets_tmp[tet_idx].tmp_Verts[index0][1]);
								std::cout << cube_index << "\t";

								cube_index = (m_storage_tets_tmp[tet_idx].tmp_Verts[index0][2]);
								std::cout << cube_index << std::endl;
							}
						}*/
						tet_count ++;
					}	
				}

			}
		}
	}

	//m_storage_tets
	for(int k=0; k<slide_Z/2; k++){
		for(int i=0; i<slide_X/2; i++){
			for(int j=0; j<slide_Y/2; j++){
				int tet_idx0 = ((k+buffer_win_startz)*(slide_X)*(slide_Y) + (i+buffer_win_startx)*(slide_Y) + (j+buffer_win_starty))*tetn;
	
				for(int index=0; index<tetn; index++){
					int tet_idx = tet_idx0 + index;

					if(m_storage_tets[tet_idx].ID != -1){
						int indexes[4];
						indexes[0] = i;
						indexes[1] = j;
						indexes[2] = k;
						indexes[3] = index;
						tetFile.write((char*) indexes, sizeof(int)*4);

						bool inside_this_cube[4] = {true, true, true, true};
						//test===============================
						int cube_index[4];

						for(int index0=0; index0<4; index0++){
							int vert_idx = (m_storage_tets[tet_idx].Verts[index0] % ((slide_X+1)*(slide_Y+1))) / (slide_Y+1);
							int vert_idy = (m_storage_tets[tet_idx].Verts[index0] % ((slide_X+1)*(slide_Y+1))) % (slide_Y+1);
							int vert_idz =  m_storage_tets[tet_idx].Verts[index0] / ((slide_X+1)*(slide_Y+1));

							inside_this_cube[index0] = (vert_idx >= (buffer_win_startx+vert_save_startx) && (vert_idx <= buffer_win_endx));
							inside_this_cube[index0] = inside_this_cube[index0] && (vert_idy >= (buffer_win_starty+vert_save_starty)) && (vert_idy <= buffer_win_endy);
							inside_this_cube[index0] = inside_this_cube[index0] && (vert_idz >= (buffer_win_startz+vert_save_startz)) && (vert_idz <= buffer_win_endz);
						}

						//bool all_inside_cube 0/1
						int all_inside;
						if(inside_this_cube[0] && inside_this_cube[1] && inside_this_cube[2] && inside_this_cube[3]){
							all_inside = 1;
							tetFile.write((char*)&all_inside, sizeof(int));
						}else{
							all_inside = 0;
							tetFile.write((char*)&all_inside, sizeof(int));
						}

						//label
						int tet_label = m_storage_tets[tet_idx].get_label();
						tetFile.write((char*)&tet_label, sizeof(int));

						//verts
						for(int index0=0; index0<4; index0++){
							int vert_idx = (m_storage_tets[tet_idx].Verts[index0] % ((slide_X+1)*(slide_Y+1))) / (slide_Y+1);
							int vert_idy = (m_storage_tets[tet_idx].Verts[index0] % ((slide_X+1)*(slide_Y+1))) % (slide_Y+1);
							int vert_idz =  m_storage_tets[tet_idx].Verts[index0] / ((slide_X+1)*(slide_Y+1));
							
							bool vert_inside_this_cube = inside_this_cube[index0];
							tetFile.write((char*)&vert_inside_this_cube, sizeof(bool));

							int tet_ids[3];
							tet_ids[0] = vert_idx-buffer_win_startx;
							tet_ids[1] = vert_idy-buffer_win_starty;
							tet_ids[2] = vert_idz-buffer_win_startz;
							tetFile.write((char*) tet_ids, sizeof(int)*3);
							
						}

						tet_count ++;
					}
					
				}

			}
		}
	}

	tetFile.write((char*) &tet_count, sizeof(int));
	tetFile.flush();
	tetFile.close();

}


//win_x, win_y, win_z means index of the small cube
void cwg::TetrahedralMesh::loadData_fromText_wenhua(int win_x, int win_y, int win_z, int win_inside_index){
	std::ifstream vertFile;
	
	int data_win_start_x = win_x * (slide_X/2);
	int data_win_start_y = win_y * (slide_Y/2);
	int data_win_start_z = win_z * (slide_Z/2);
	
	int buffer_win_x_start = ((win_inside_index%4)/2) *(slide_X/2);
	int buffer_win_y_start = (win_inside_index%2) * (slide_Y/2);
	int buffer_win_z_start = (win_inside_index/4) * (slide_Z/2);

	int data_cube_start_x = (win_x - ((win_inside_index%4)/2)) * (slide_X/2);
	int data_cube_start_y = (win_y - (win_inside_index%2))     * (slide_Y/2);
	int data_cube_start_z = (win_z - (win_inside_index/4))     * (slide_Z/2);

	reset_win_verts(win_inside_index, data_cube_start_x, data_cube_start_y, data_cube_start_z);

	int data_buffer_startx = data_win_start_x - buffer_win_x_start;
	int data_buffer_starty = data_win_start_y - buffer_win_y_start;
	int data_buffer_startz = data_win_start_z - buffer_win_z_start;

//	reset_win_verts(win_inside_index, false);
	//this will be shared by both vert file and the tet file	
	std::string tmpfilename("");
	tmpfilename += std::to_string(win_x);
	tmpfilename += "_";
	tmpfilename += std::to_string(win_y);
	tmpfilename += "_";
	tmpfilename += std::to_string(win_z);
	tmpfilename += ".txt";

	std::string vertFileName = get_foldername() + "Data\\vert_";
	vertFileName += tmpfilename;

	vertFile.open(vertFileName, ios::binary);
	std::string data;
	if(vertFile.is_open()){
		int i, j, k, keep_number;
		double vert_poses[3];
		int label_size, label_int;
		vertFile.seekg(-(sizeof(int)), vertFile.end);
		int vert_size;
		vertFile.read((char*)&vert_size, sizeof(int));
		vertFile.seekg(0, vertFile.beg);

		for(int vert_count=0; vert_count<vert_size; vert_count++){
			//idx idy idz
			vertFile.read((char*)&i, sizeof(int));
			vertFile.read((char*)&j, sizeof(int));
			vertFile.read((char*)&k, sizeof(int));
			int idx = (buffer_win_z_start+k)*(slide_X+1)*(slide_Y+1) + (buffer_win_x_start+i)*(slide_Y+1) + (buffer_win_y_start+j);

			//keep_number
			vertFile.read((char*)&keep_number, sizeof(int));
			m_storage_verts[idx].keep = keep_number;

			//pos.x pos.y pos.z
			vertFile.read((char*)vert_poses, sizeof(double)*3);
			m_storage_verts[idx].pos() = Cleaver::vec3(vert_poses[0], vert_poses[1], vert_poses[2]);
			m_storage_verts[idx].ID = idx;

			//labels
			vertFile.read((char*)&label_size, sizeof(int));

			for(int label_index=0; label_index<label_size; label_index++){			
				vertFile.read((char*)&label_int, sizeof(int));
				m_storage_verts[idx].assign_label(label_int);
			}
		}
	}else{
		std::cout << "vert file can not open! so it can not be read!" << std::endl;
	}

	vertFile.close();

	std::ifstream tetFile;

	std::string tetFileName  = get_foldername() + "Data\\tet_";
	tetFileName += tmpfilename;
	tetFile.open(tetFileName, ios::binary);

	if(tetFile.is_open()){
		tetFile.seekg(-(sizeof(int)), tetFile.end);
		int tet_size;
		tetFile.read((char*)&tet_size, sizeof(int));
		tetFile.seekg(0, tetFile.beg);

		for(int tet_count=0; tet_count<tet_size; tet_count++){
			int x, y, z, index0;
			tetFile.read((char*)&x, sizeof(int));
			tetFile.read((char*)&y, sizeof(int));
			tetFile.read((char*)&z, sizeof(int));
			tetFile.read((char*)&index0, sizeof(int));

			int idx = ((buffer_win_z_start+z)*(slide_X)*(slide_Y) + (buffer_win_x_start+x)*(slide_Y) + (buffer_win_y_start+y)) * tetn + index0;
			//bool all_inside_cube
			int all_inside_cube_int;
			tetFile.read((char*)&all_inside_cube_int, sizeof(int));

			if(all_inside_cube_int == 1){ //all the verts of this tet is inside the cube
				//label

				int label_int;
				tetFile.read((char*)&label_int, sizeof(int));
				m_storage_tets[idx].set_label(label_int);

				//id
				if(m_storage_tets[idx].ID != -1)
					throw exception("Really? Wrong here again? Please! load from text to tet!");

				m_storage_tets[idx].ID = idx;

				//verts
				int vert_idx_int, vert_idy_int, vert_idz_int;
				bool inside_cube_int;

				//verts
				for(int i=0; i<4; i++){
					tetFile.read((char*)&inside_cube_int, sizeof(bool));
					tetFile.read((char*)&vert_idx_int, sizeof(int));
					tetFile.read((char*)&vert_idy_int, sizeof(int));
					tetFile.read((char*)&vert_idz_int, sizeof(int));
					m_storage_tets[idx].Verts[i] = (vert_idz_int+buffer_win_z_start)*(slide_X+1)*(slide_Y+1) + (vert_idx_int+buffer_win_x_start)*(slide_Y+1) + (vert_idy_int+buffer_win_y_start);

				}
			}
			else{ //some of the verts of the tet is not inside the cube				
				if(m_storage_tets_tmp[idx].ID != -1)
					idx = find_available_tmp_tet_id(buffer_win_x_start+x, buffer_win_y_start+y, buffer_win_z_start+z);

				//label 
				int label_int;
				tetFile.read((char*)&label_int, sizeof(int));
				m_storage_tets_tmp[idx].set_label(label_int);


				m_storage_tets_tmp[idx].ID = idx;
				int vert_idx_int, vert_idy_int, vert_idz_int;
				bool keep_int;
				
				//verts
				for(int i=0; i<4; i++){
					tetFile.read((char*)&keep_int, sizeof(bool));
					tetFile.read((char*)&vert_idx_int, sizeof(int));
					tetFile.read((char*)&vert_idy_int, sizeof(int));
					tetFile.read((char*)&vert_idz_int, sizeof(int));
					m_storage_tets_tmp[idx].tmp_Verts[i][0] = vert_idx_int + buffer_win_x_start;
					m_storage_tets_tmp[idx].tmp_Verts[i][1]  = vert_idy_int + buffer_win_y_start;
					m_storage_tets_tmp[idx].tmp_Verts[i][2]  = vert_idz_int + buffer_win_z_start;
				}
			}
			
		}
	}else{
		std::cout << "tet file can not open! so it can not be read!" << std::endl;
	}

	tetFile.close();
	
	

}


//This cube means the big cube
void cwg::TetrahedralMesh::loadData_fromFile_wenhua(const std::string& foldername , const std::vector<std::string>& filelist, int Data_cube_start_x, int Data_cube_start_y, int Data_cube_start_z, int win_inside_index){
	

	int Buffer_win_start_x = (slide_X/2) * ( (win_inside_index % 4) / 2);
	int Buffer_win_start_y = (slide_Y/2) * ( win_inside_index % 2      );
	int Buffer_win_start_z = (slide_Z/2) * ( win_inside_index / 4      );

	int data_win_start_x = Data_cube_start_x + (slide_X/2) * ( (win_inside_index % 4) / 2);
	int data_win_start_y = Data_cube_start_y + (slide_Y/2) * ( win_inside_index % 2);
	int data_win_start_z = Data_cube_start_z + (slide_Z/2) * ( win_inside_index / 4);

	int data_win_end_x = data_win_start_x + steps[0];
	int data_win_end_y = data_win_start_y + steps[1];
	int data_win_end_z = data_win_start_z + steps[2];

	//Reset the verts and make a clear verts environment for future use
	reset_win_verts(win_inside_index, Data_cube_start_x, Data_cube_start_y, Data_cube_start_z);

	int label;
	vector<Tetrahedron>::iterator it = m_storage_tets.begin();

	for (int k = data_win_start_z; k < data_win_end_z; k++)
	{
		string& filename = foldername+filelist[k ];
		ifstream ifile;
		ifile.open( filename.c_str() , ios::binary );
		ifile.seekg((ifile.beg + 3 + data_win_start_x*(m_Y) + data_win_start_y)*sizeof(int));	//skip first 3 values (size)

		for (int j = data_win_start_x; j < data_win_end_x; j++)
		{
			for (int i = data_win_start_y; i < data_win_end_y; i++)
			{
				ifile.read( (char*)&label , sizeof(int) );
				TetCube tmp_cube(label);
				
				int idx[8] = {
								(slide_Y+1)*(slide_X+1)*(k - data_win_start_z + Buffer_win_start_z ) + (slide_Y+1)*(j - data_win_start_x + Buffer_win_start_x ) +(i  - data_win_start_y + Buffer_win_start_y),
								(slide_Y+1)*(slide_X+1)*(k - data_win_start_z + Buffer_win_start_z ) + (slide_Y+1)*(j - data_win_start_x + Buffer_win_start_x ) +(i+1- data_win_start_y + Buffer_win_start_y),
								(slide_Y+1)*(slide_X+1)*(k - data_win_start_z + Buffer_win_start_z ) + (slide_Y+1)*(j+1-data_win_start_x + Buffer_win_start_x) + (i  - data_win_start_y + Buffer_win_start_y),
								(slide_Y+1)*(slide_X+1)*(k - data_win_start_z + Buffer_win_start_z ) + (slide_Y+1)*(j+1-data_win_start_x + Buffer_win_start_x) + (i+1- data_win_start_y + Buffer_win_start_y),

								(slide_Y+1)*(slide_X+1)*(k+1- data_win_start_z + Buffer_win_start_z) + (slide_Y+1)*(j - data_win_start_x + Buffer_win_start_x ) +(i  - data_win_start_y + Buffer_win_start_y),
								(slide_Y+1)*(slide_X+1)*(k+1- data_win_start_z + Buffer_win_start_z) + (slide_Y+1)*(j - data_win_start_x + Buffer_win_start_x ) +(i+1- data_win_start_y + Buffer_win_start_y),
								(slide_Y+1)*(slide_X+1)*(k+1- data_win_start_z + Buffer_win_start_z) + (slide_Y+1)*(j+1-data_win_start_x + Buffer_win_start_x) + (i  - data_win_start_y + Buffer_win_start_y),
								(slide_Y+1)*(slide_X+1)*(k+1- data_win_start_z + Buffer_win_start_z) + (slide_Y+1)*(j+1-data_win_start_x + Buffer_win_start_x) + (i+1- data_win_start_y + Buffer_win_start_y),							
							};

				for(int index=0; index<8; index++){
					m_storage_verts[idx[index]].assign_label(label);
					m_storage_verts[idx[index]].ID = idx[index];
				}

				//set v to cube
				for(int index=0; index<8; index++){
					tmp_cube.Verts[index] = idx[index];
				}

				//tetrahedronlize
				int tidx = (slide_Y)*(slide_X)*(k - data_win_start_z + Buffer_win_start_z) + (slide_Y)*(j- data_win_start_x + Buffer_win_start_x) +  (i- data_win_start_y + Buffer_win_start_y);
				
				int tID = tidx * tetn;
				it = m_storage_tets.begin() + tID;

				tetrahedronlize(i, j, k, tmp_cube, it, tID);


			}

			ifile.seekg( (m_Y - (data_win_end_y - data_win_start_y) )*sizeof(int),ios_base::cur);	// next line

		}
		ifile.close();
	}

}


void cwg::TetrahedralMesh::loadData_fromRight_wenhua(int startx, int starty, int startz){
	reset_win_verts(0, startx, starty, startz); 
	reset_win_verts(1, startx, starty, startz); 
	reset_win_verts(4, startx, starty, startz); 
	reset_win_verts(5, startx, starty, startz);
	//Reset and make a clear environment

	int from_x_start = slide_X/2;

	//vert
	for(int k=0;  k<=slide_Z; k++){
		for(int i=from_x_start+1; i<=slide_X; i++){
			for(int j=0; j<=slide_Y; j++){

				int idx_from = k*(slide_X+1)*(slide_Y+1) + i*(slide_Y+1) + j;
				int idx_to   = k*(slide_X+1)*(slide_Y+1) + (i-(slide_X/2))*(slide_Y+1) + j;

				if(m_storage_verts[idx_from].ID == -1){
					m_storage_verts[idx_to].ID = -1;
				}else{
					m_storage_verts[idx_to].ID = idx_to;
					m_storage_verts[idx_to].pos().x = m_storage_verts[idx_from].pos().x;
					m_storage_verts[idx_to].pos().y = m_storage_verts[idx_from].pos().y;
					m_storage_verts[idx_to].pos().z = m_storage_verts[idx_from].pos().z;
					m_storage_verts[idx_to].m_labels = m_storage_verts[idx_from].m_labels;
					m_storage_verts[idx_to].keep = m_storage_verts[idx_from].keep;
				}
			}
		}
	}


	//tmp_tet
	for(int k=0;  k<slide_Z; k++){
		for(int i=from_x_start; i<slide_X; i++){
			for(int j=0; j<slide_Y; j++){
				int idx_from = (k*slide_X*slide_Y + i*slide_Y+ j) * tetn;
				int idx_to   = (k*slide_X*slide_Y + (i-(slide_X/2))*slide_Y + j) * tetn;

				for(int index=0; index < tetn; index++){
					if(m_storage_tets_tmp[idx_from+index].ID != -1){
						//set Id and label
						if(m_storage_tets_tmp[idx_to+index].ID == -1){
							m_storage_tets_tmp[idx_to+index].ID = idx_to+index;
							m_storage_tets_tmp[idx_to+index].set_label( m_storage_tets_tmp[idx_from+index].get_label());

							for(int index0=0; index0<4; index0++){			
								m_storage_tets_tmp[idx_to+index].tmp_Verts[index0][0] = m_storage_tets_tmp[idx_from+index].tmp_Verts[index0][0] - slide_X/2;
								m_storage_tets_tmp[idx_to+index].tmp_Verts[index0][1] = m_storage_tets_tmp[idx_from+index].tmp_Verts[index0][1];
								m_storage_tets_tmp[idx_to+index].tmp_Verts[index0][2] = m_storage_tets_tmp[idx_from+index].tmp_Verts[index0][2];
							}
						}else{
							int tet_id = find_available_tmp_tet_id(i-slide_X/2, j, k);

							m_storage_tets_tmp[tet_id].ID = tet_id;
							m_storage_tets_tmp[tet_id].set_label( m_storage_tets_tmp[idx_from+index].get_label());

							for(int index0=0; index0<4; index0++){			
								m_storage_tets_tmp[tet_id].tmp_Verts[index0][0] = m_storage_tets_tmp[idx_from+index].tmp_Verts[index0][0] - slide_X/2;
								m_storage_tets_tmp[tet_id].tmp_Verts[index0][1] = m_storage_tets_tmp[idx_from+index].tmp_Verts[index0][1];
								m_storage_tets_tmp[tet_id].tmp_Verts[index0][2] = m_storage_tets_tmp[idx_from+index].tmp_Verts[index0][2];
							}
				
						}

					}
				}


			}
		}
	}


	//tet
	for(int k=0;  k<slide_Z; k++){
		for(int i=from_x_start; i<slide_X; i++){
			for(int j=0; j<slide_Y; j++){
				int idx_from = (k*slide_X*slide_Y + i*slide_Y+ j) * tetn;
				int idx_to   = (k*slide_X*slide_Y + (i-(slide_X/2))*slide_Y + j) * tetn;

				for(int index=0; index < tetn; index++){

					if(m_storage_tets[idx_from+index].ID != -1){
						bool right_halfs[4] = {true, true, true, true};
						
						int diff = (slide_X/2) * (slide_Y+1);

						bool right_half = true;
						//deal with the tets which contains verts outside the cube, i.e. the origin vert whose x is less than slide_X/2
						for(int index0=0; index0<4; index0++){
							right_halfs[index0] =  ((m_storage_tets[idx_from+index].Verts[index0]%((slide_X+1)*(slide_Y+1)))/(slide_Y+1)) > from_x_start; //x
							right_halfs[index0] =  right_halfs[index0] && ((m_storage_tets[idx_from+index].Verts[index0] % ((slide_X+1)*(slide_Y+1))) % (slide_Y+1)) >= 0;
							right_halfs[index0] =  right_halfs[index0] && (m_storage_tets[idx_from+index].Verts[index0] / ((slide_X+1)*(slide_Y+1))) >= 0;
							right_half = right_half && ( right_halfs[index0]);
						}

						//If all the verts are in the small cube
						if(right_half){
							m_storage_tets[idx_to+index].set_label(m_storage_tets[idx_from+index].get_label());
							m_storage_tets[idx_to+index].ID = idx_to+index;

							for(int index0=0; index0<4; index0++){
								m_storage_tets[idx_to+index].Verts[index0] = m_storage_tets[idx_from+index].Verts[index0] - diff;	
							}
						}else{//if something should be stored in the tmp

							//m_storage_tets[idx_to+index].ID = -1;

							//find the index for the tmp tet
							int m_storage_tets_tmp_index = idx_to+index;

							if(m_storage_tets_tmp[m_storage_tets_tmp_index].ID != -1)
								m_storage_tets_tmp_index = find_available_tmp_tet_id(i-slide_X/2, j, k);

							m_storage_tets_tmp[m_storage_tets_tmp_index].set_label(m_storage_tets[idx_from+index].get_label());
							m_storage_tets_tmp[m_storage_tets_tmp_index].ID = m_storage_tets_tmp_index;

							for(int index0=0; index0<4; index0++){
								int vert_idx = (m_storage_tets[idx_from+index].Verts[index0] % ((slide_X+1)*(slide_Y+1))) / (slide_Y+1);
								int vert_idy = (m_storage_tets[idx_from+index].Verts[index0] % ((slide_X+1)*(slide_Y+1))) % (slide_Y+1);
								int vert_idz =  m_storage_tets[idx_from+index].Verts[index0] / ((slide_X+1)*(slide_Y+1)); 
								
								m_storage_tets_tmp[m_storage_tets_tmp_index].tmp_Verts[index0][0] = vert_idx - slide_X/2;
								m_storage_tets_tmp[m_storage_tets_tmp_index].tmp_Verts[index0][1] = vert_idy;
								m_storage_tets_tmp[m_storage_tets_tmp_index].tmp_Verts[index0][2] = vert_idz;
							}

						}

					}


				}
			}
		}
	}

} 


int cwg::TetrahedralMesh::find_available_tet_id(int orig_tet_x, int orig_tet_y, int orig_tet_z){
	int tet_id;

	int cube_indexes[3];

	cube_indexes[0] = (orig_tet_x)/(slide_X/2);
	cube_indexes[1] = (orig_tet_y)/(slide_Y/2);
	cube_indexes[2] = (orig_tet_z)/(slide_Z/2);

	int tet_start_index[3] = {cube_indexes[0]*(slide_X/2), cube_indexes[1]*(slide_Y/2), cube_indexes[2]*(slide_Z/2)};
	int tet_end_index[3]   = {(cube_indexes[0]+1)*(slide_X/2), (cube_indexes[1]+1)*(slide_Y/2), (cube_indexes[2]+1)*(slide_Z/2)};

	bool found_one = false;
	for(int tet_adjust_k=tet_start_index[2]; tet_adjust_k<tet_end_index[2]; tet_adjust_k++){
		for(int tet_adjust_i=tet_start_index[0]; tet_adjust_i<tet_end_index[0]; tet_adjust_i++){
			for(int tet_adjust_j=tet_start_index[1]; tet_adjust_j<tet_end_index[1]; tet_adjust_j++){

				int tet_id0 = (tet_adjust_k*slide_X*slide_Y + tet_adjust_i*slide_Y + tet_adjust_j)*tetn;

				for(int tet_adjust_index=0; tet_adjust_index<tetn; tet_adjust_index++){
					tet_id = tet_id0 + tet_adjust_index;

					if(m_storage_tets[tet_id].ID == -1)
							found_one = true;

					if(found_one) break;
				}
				if(found_one) break;
			}
			if(found_one) break;
		}
		if(found_one) break;
	 }

	if(!found_one){
		tet_id = -1;
		throw exception("all the tets are occupied and can not accomodate the one from m_storage_tet_tmp");
	}
		
	return tet_id;
}


int cwg::TetrahedralMesh::find_available_tmp_tet_id(int orig_tmp_tet_x, int orig_tmp_tet_y, int orig_tmp_tet_z){
	int tet_id;

	int cube_indexes[3];

	cube_indexes[0] = (orig_tmp_tet_x)/(slide_X/2);
	cube_indexes[1] = (orig_tmp_tet_y)/(slide_Y/2);
	cube_indexes[2] = (orig_tmp_tet_z)/(slide_Z/2);

	int tet_start_index[3] = {cube_indexes[0]*(slide_X/2), cube_indexes[1]*(slide_Y/2), cube_indexes[2]*(slide_Z/2)};
	int tet_end_index[3]   = {(cube_indexes[0]+1)*(slide_X/2), (cube_indexes[1]+1)*(slide_Y/2), (cube_indexes[2]+1)*(slide_Z/2)};

	bool found_one = false;
	for(int tet_adjust_k=tet_start_index[2]; tet_adjust_k<tet_end_index[2]; tet_adjust_k++){
		for(int tet_adjust_i=tet_start_index[0]; tet_adjust_i<tet_end_index[0]; tet_adjust_i++){
			for(int tet_adjust_j=tet_start_index[1]; tet_adjust_j<tet_end_index[1]; tet_adjust_j++){

				int tet_id0 = (tet_adjust_k*slide_X*slide_Y + tet_adjust_i*slide_Y + tet_adjust_j)*tetn;

				for(int tet_adjust_index=0; tet_adjust_index<tetn; tet_adjust_index++){
					tet_id = tet_id0 + tet_adjust_index;

					if(m_storage_tets_tmp[tet_id].ID == -1)
							found_one = true;

					if(found_one) break;
				}
				if(found_one) break;
			}
			if(found_one) break;
		}
		if(found_one) break;
	 }

	if(!found_one){
		tet_id = -1;
		throw exception("all the tets are occupied and can not accomodate the one from m_storage_tet_tmp");
	}
		
	return tet_id;
}


int cwg::TetrahedralMesh::find_adjust_tet_id(int vert[4][3]){
	int tet_id;

	bool vert_tested[4] = {false, false, false, false};
	int vert_x, vert_y, vert_z;

	int smallest_x = 1e10;
	int smallest_y = 1e10;
	int smallest_x_index = -1;
	int smallest_y_index = -1;

	for(int i=0; i<4; i++){
		if(vert[i][0] < smallest_x){
			smallest_x = vert[i][0];
			smallest_x_index = i;
		}
	}


	while(!vert_tested[0] || !vert_tested[1] || !vert_tested[2] || !vert_tested[3]){
			
		int vert_index = rand() % 4;		
			if(vert_tested[vert_index])
				continue;

		int vert_x = vert[vert_index][0];
		int vert_y = vert[vert_index][1];
		int vert_z = vert[vert_index][2];


		int cube_indexes[3];
	
		cube_indexes[0] = (vert_x - 1) / (slide_X/2);
		cube_indexes[0] = cube_indexes[0] < 0 ? 0 : cube_indexes[0]; 

		cube_indexes[1] = (vert_y - 1) / (slide_Y/2);
		cube_indexes[1] = cube_indexes[1] < 0 ? 0 : cube_indexes[1]; 

		cube_indexes[2] = (vert_z - 1) / (slide_Z/2);
		cube_indexes[2] = cube_indexes[2] < 0 ? 0 : cube_indexes[2]; 

		int cube_index = cube_indexes[2] * 2 * 2 + cube_indexes[0]*2 + cube_indexes[1];

		if( start_index[cube_index][2] == (cube_index/(2 * 2)+1)*(slide_Z/2)){
			vert_tested[vert_index] = true;
			continue;
		}

		tet_id = (start_index[cube_index][2]*slide_X*slide_Y + start_index[cube_index][0]*slide_Y + start_index[cube_index][1])*tetn + start_index[cube_index][3];

		if(start_index[cube_index][3] != (tetn-1)){
			start_index[cube_index][3] ++;
		}else{
			start_index[cube_index][3] = 0;
		
			if( ((start_index[cube_index][1]+1) % (slide_Y/2)) != 0 ){
				start_index[cube_index][1] ++;
			}else{
				start_index[cube_index][1] = (start_index[cube_index][1] / (slide_Y/2)) * (slide_Y/2);

				if( (start_index[cube_index][0]+1) % (slide_X/2) != 0){
					start_index[cube_index][0]++;
				}else{
					start_index[cube_index][0] = (start_index[cube_index][0] / (slide_X/2)) * (slide_X/2);
					start_index[cube_index][2]++;
							
				}
			}
		}
		break;
	}
	
	if(vert_tested[0] && vert_tested[1] && vert_tested[2] && vert_tested[3])
		throw exception("Really?? Wrong Tet Adjust! Please use other parameters!");

	return tet_id;
}


void cwg::TetrahedralMesh::loadData_betweenCubes_wenhua(){
	//to test

	int buffer_length[3] = {slide_X, slide_Y, slide_Z};

	int size = m_storage_tets_tmp.size();

	//deal with the verts located between cubes
	//if the verts are all inside the big cube, make the tet valid
	for(int i=0; i<size; i++){
		if(m_storage_tets_tmp[i].ID != -1){

			bool inside_cube[4] = {true, true, true, true};

			bool all_inside = true;
			for(int index0=0; index0<4; index0++){
			
				int vert_ids[3] = { m_storage_tets_tmp[i].tmp_Verts[index0][0], 
										m_storage_tets_tmp[i].tmp_Verts[index0][1], 				   
										m_storage_tets_tmp[i].tmp_Verts[index0][2]};
				
				//std::cout << "vert_id: " << vert_ids[0] <<  " " << vert_ids[1] << " " << vert_ids[2] << std::endl;
				for(int dimen_i=0; dimen_i<3; dimen_i++){
					inside_cube[index0] = inside_cube[index0] && (vert_ids[dimen_i]>=0) && (vert_ids[dimen_i] <= buffer_length[dimen_i]);
				}

				if(inside_cube[index0]){
					int vert_id = vert_ids[2]*(buffer_length[0]+1)*(buffer_length[1]+1) + vert_ids[0]*(buffer_length[1]+1) + vert_ids[1];
					all_inside = all_inside && (m_storage_verts[vert_id].ID != -1);	
				}else{
					all_inside = false;
				}
			}
 
			if(all_inside){
				if(m_storage_tets[i].ID == -1){
					m_storage_tets[i].ID = i;
					m_storage_tets[i].set_label(m_storage_tets_tmp[i].get_label());
					m_storage_tets_tmp[i].ID = -1;

					for(int index0=0; index0<4; index0++){
						int vert_idx = m_storage_tets_tmp[i].tmp_Verts[index0][0];
						int vert_idy = m_storage_tets_tmp[i].tmp_Verts[index0][1];
						int vert_idz = m_storage_tets_tmp[i].tmp_Verts[index0][2];
					
						int vert_id = vert_idz*(buffer_length[0]+1)*(buffer_length[1]+1) + vert_idx*(buffer_length[1]+1) + vert_idy;
						m_storage_tets[i].Verts[index0] = vert_id;
						m_storage_verts[vert_id].keep --;

					}
				}
				else{
					int tet_x = ((i/tetn) % (slide_X * slide_Y)) / slide_Y;
					int tet_y = ((i/tetn) % (slide_X * slide_Y)) % slide_Y;
					int tet_z = ((i/tetn) / (slide_X * slide_Y));

					int tet_id = find_available_tet_id(tet_x, tet_y, tet_z);
					
					m_storage_tets[tet_id].ID = tet_id;
					m_storage_tets[tet_id].set_label(m_storage_tets_tmp[i].get_label());
					m_storage_tets_tmp[i].ID = -1;


					for(int index0=0; index0<4; index0++){
						int vert_idx = m_storage_tets_tmp[i].tmp_Verts[index0][0];
						int vert_idy = m_storage_tets_tmp[i].tmp_Verts[index0][1];
						int vert_idz = m_storage_tets_tmp[i].tmp_Verts[index0][2];
					
						int vert_id = vert_idz*(buffer_length[0]+1)*(buffer_length[1]+1) + vert_idx*(buffer_length[1]+1) + vert_idy;
						m_storage_tets[tet_id].Verts[index0] = vert_id;
						m_storage_verts[vert_id].keep --;
					}
							
				} // if(m_storage_tets[i].ID != -1)
			}// if(all_inside)
			
		}
	}

}


void cwg::TetrahedralMesh::saveData(int data_cube_start_x, int data_cube_start_y, int data_cube_start_z){

	bool x_end = ((data_cube_start_x + slide_X/2 +steps[0]) == m_X);

	//bool x_end = ((data_cube_start_x + slide_X/2 +steps[0]) == slide_X);
	if(!x_end){
		for(int k=0; k<=1; k++){
			for(int j=0; j<=1; j++){
				for(int i=0; i<=1;i++){		
					
					bool dont_save = (i==1);
					if(dont_save)
						continue;

					saveResultData_wenhua(data_cube_start_x, data_cube_start_y,data_cube_start_z, i, j, k);
				}
			}
		}
	}else{

		for(int k=0; k<=1; k++){
			for(int j=0; j<=1; j++){
				for(int i=0; i<=1;i++){

					saveResultData_wenhua(data_cube_start_x, data_cube_start_y,data_cube_start_z, i, j, k);
				}
			}
		}
	}


	//std::cout << "Saved the Data" << std::endl;
}


//Vert file format
//==================================================================================
//x y z
//label0
//label1
//

//Tet file format
//=================================================================================
//tet_idx
//vert0.idx
//vert1.idx
//vert2.idx
//vert3.idx
//=

void cwg::TetrahedralMesh::write_obj_for_debug(int tet_id){
	string keep_file_name = std::to_string(ResultT[tet_id].get_label()) + "_" + std::to_string(tet_id) + ".obj";
	std::cout << keep_file_name << std::endl;
	ofstream output_debug(keep_file_name);

	for(int vert_idx=0; vert_idx<4; vert_idx++){
		output_debug << "v " << ResultV[ResultT[tet_id].Verts[vert_idx]].pos().x << " " << ResultV[ResultT[tet_id].Verts[vert_idx]].pos().y << " "
				     << ResultV[ResultT[tet_id].Verts[vert_idx]].pos().z << std::endl;
	}

	for(int face_idx=0; face_idx<4; face_idx++){
		output_debug << "f " << (face_idx) % 4 + 1 << " " << (face_idx+1)%4 +1 << " " << (face_idx+2)%4 + 1<< std::endl; 
	}
	
	output_debug.close();

	std::cout << "finished outputting " << std::endl;

}


void cwg::TetrahedralMesh::write_obj_for_debug2(int x, int y, int z){
	int vert_size = m_storage_verts.size();
	std::vector<int> all_labels;

	for(int vert_count=0; vert_count<vert_size; vert_count++){
		if (m_storage_verts[vert_count].ID != -1){
			for(int label_index=0; label_index < m_storage_verts[vert_count].NLabels(); label_index++){
				int lb = m_storage_verts[vert_count].get_label(label_index);
				std::vector<int>::iterator it = std::lower_bound(all_labels.begin(), all_labels.end(), lb);
				if (it != all_labels.end()){
					if (*it != lb)
					all_labels.insert(it, lb);
				}
				else
					all_labels.push_back(lb);
				}
		}	
	}

	for(int label_idx=0; label_idx<all_labels.size(); label_idx++){
		int lb = all_labels[label_idx];
		string filename = "m_storage_" + std::to_string(x) + "_" + std::to_string(y) + "_" + std::to_string(z) + "_" + std::to_string(lb) + ".obj";
		ofstream output_single_label(filename);
		int tet_size = m_storage_tets.size();
		int valid_tets = 0;
		for(int tet_count=0; tet_count<tet_size; tet_count++){
			if(m_storage_tets[tet_count].ID != -1){
				if(m_storage_tets[tet_count].get_label() == lb){
					for(int vert_index=0; vert_index<4; vert_index++){
						output_single_label << "v " << std::to_string(m_storage_verts[m_storage_tets[tet_count].Verts[vert_index]].pos().x);
						output_single_label << " " << std::to_string(m_storage_verts[m_storage_tets[tet_count].Verts[vert_index]].pos().y);
						output_single_label << " " << std::to_string(m_storage_verts[m_storage_tets[tet_count].Verts[vert_index]].pos().z) << std::endl;
					}
					valid_tets ++;
				}
			}
		}

		for(int tet_count=0; tet_count < valid_tets; tet_count++){
			output_single_label << "l "  << tet_count*4+1 << " " << tet_count*4+2 << std::endl; 
			output_single_label << "l "  << tet_count*4+1 << " " << tet_count*4+3 << std::endl; 
			output_single_label << "l "  << tet_count*4+1 << " " << tet_count*4+4 << std::endl;
			output_single_label << "l "  << tet_count*4+2 << " " << tet_count*4+3 << std::endl; 
			output_single_label << "l "  << tet_count*4+2 << " " << tet_count*4+4 << std::endl; 
			output_single_label << "l "  << tet_count*4+3 << " " << tet_count*4+4 << std::endl;
		}

		output_single_label.close();
	}
}


void cwg::TetrahedralMesh::vert_adjust(){

	std::map<int, int>::iterator it;
	int vert_size = ResultV_tmp.size();
	int vert_count=0;
	for(vert_count=0; vert_count<vert_size; vert_count++){
		m_storage_verts[vert_window_map_inverse[vert_count]].pos() = Cleaver::vec3(ResultV[vert_count].pos().x, 
														ResultV[vert_count].pos().y, 
														ResultV[vert_count].pos().z);
	}

	vert_size = ResultV.size();
	for(; vert_count<vert_size; vert_count++){
		m_storage_verts[vert_window_map_inverse[vert_count]].pos() = Cleaver::vec3(ResultV[vert_count].pos().x, 
														ResultV[vert_count].pos().y, 
														ResultV[vert_count].pos().z);
		m_storage_verts[vert_window_map_inverse[vert_count]].ID = vert_window_map_inverse[vert_count];
		for (int label_count=0; label_count < ResultV[vert_count].NLabels(); label_count++)
			m_storage_verts[vert_window_map_inverse[vert_count]].assign_label(ResultV[vert_count].get_label(label_count));
		
		m_storage_verts[vert_window_map_inverse[vert_count]].keep = 0;
		
	}

	int tet_size = m_storage_tets.size();
	int result_T_size = ResultT.size();
	int m_tet_count = 0;
	for(int result_tet_count=0; result_tet_count<result_T_size; result_tet_count++){
		while (m_storage_tets[m_tet_count].ID != -1 && m_tet_count<tet_size){
			m_tet_count++;
		}

		if(m_tet_count == tet_size){
			throw exception("Full m storage tets!");
		}

		m_storage_tets[m_tet_count].ID = m_tet_count;

		for(int vert_index=0; vert_index<4; vert_index++){
			m_storage_tets[m_tet_count].Verts[vert_index] = vert_window_map_inverse[ResultT[result_tet_count].Verts[vert_index]];
		}

		m_storage_tets[m_tet_count].set_label(ResultT[result_tet_count].get_label());
	}

}


bool cwg::TetrahedralMesh::tet_boundary_check(int tet_index){
	bool is_on_boundary = true;
	is_on_boundary = (vert_boundary_check( m_storage_tets[tet_index].Verts[0])  || vert_boundary_check( m_storage_tets[tet_index].Verts[1]) ||
							vert_boundary_check( m_storage_tets[tet_index].Verts[2]) || vert_boundary_check( m_storage_tets[tet_index].Verts[3]));
	return is_on_boundary;
}


bool cwg::TetrahedralMesh::vert_boundary_check(int vert_id){

	return (m_storage_verts[vert_id].boundary || m_storage_verts[vert_id].slideBoundary || (m_storage_verts[vert_id].keep > 0));
}


bool cwg::TetrahedralMesh::result_tris_boundary_check(int vert0, int vert1, int vert2){
	bool is_on_boundary = ((ResultV[vert0].slideBoundary && ResultV[vert1].slideBoundary && ResultV[vert2].slideBoundary) ||
								( ResultV[vert0].boundary && ResultV[vert1].boundary && ResultV[vert2].boundary));
	return is_on_boundary;
}


void cwg::TetrahedralMesh::getResultV_T_window_test(Visualize& vis){

	ResultV.clear();
	ResultT.clear();
	ResultV_tmp.clear();
	ResultT_tmp.clear();
	vert_window_map.clear();
	tet_window_map.clear();
	vert_window_map_inverse.clear();

	TetVertex vert;	
	std::map<int, int>::iterator it;
	std::map<int, bool> vert_save_map;
	std::map<int, bool>::iterator vert_save_mapt_it;
	
	int vert_size = m_storage_verts.size();

	for(int vert_count=0; vert_count<vert_size; vert_count++){
		if(m_storage_verts[vert_count].ID != -1){
			vert_save_mapt_it = vert_save_map.find(vert_count);
			if(vert_save_mapt_it == vert_save_map.end()){
				vert_save_map.insert( std::pair<int, bool>(vert_count, false));
			}
		}
	}

	//Tet files
	int tet_size = m_storage_tets.size();
	for(int tet_count=0; tet_count<tet_size; tet_count++){
		if((m_storage_tets[tet_count].ID != -1) && (!tet_boundary_check(tet_count))){
			bool all_right = true;
			for(int index0=0; index0<4; index0++){
				int vert_index = m_storage_tets[tet_count].Verts[index0];
				vert_save_mapt_it = vert_save_map.find(vert_index);

				if(vert_save_mapt_it == vert_save_map.end()){
					all_right = false;
					break;
				}								 
			}

			if(all_right){
				for(int index0=0; index0<4; index0++){
					int vert_index = m_storage_tets[tet_count].Verts[index0];
					vert_save_mapt_it = vert_save_map.find(vert_index);
					vert_save_mapt_it->second = true;
				}

			}
		}
	}

	int v_id_now = ResultV.size();	

	for(int vert_count=0; vert_count<vert_size; vert_count++){
		if(m_storage_verts[vert_count].ID != -1  && vert_save_map[vert_count]){
			v_id_now = ResultV.size();
			it = vert_window_map.find(vert_count);

			if(it == vert_window_map.end()){
				vert_window_map.insert( std::pair<int, int>(vert_count, v_id_now) );
				vert_window_map_inverse.push_back(vert_count);
				vert.ID = v_id_now;

				//pos.x pos.y pos.z
				vert.pos() = Cleaver::vec3(m_storage_verts[vert_count].pos().x, 
										   m_storage_verts[vert_count].pos().y, 
										   m_storage_verts[vert_count].pos().z);

				vert.keep = m_storage_verts[vert_count].keep;
				vert.slideBoundary = m_storage_verts[vert_count].slideBoundary;
				vert.boundary = m_storage_verts[vert_count].boundary;

				int label_size = m_storage_verts[vert_count].m_labels.size();
				for(int label_index=0; label_index<label_size; label_index++)
					vert.assign_label(m_storage_verts[vert_count].m_labels[label_index]);

				ResultV.push_back(vert);
				ResultV_tmp.push_back(vert);
			}
		}
	}

	//Tet files
	int t_id_now = ResultT.size();
	tet_size = m_storage_tets.size();
	for(int tet_count=0; tet_count<tet_size; tet_count++){
		if((m_storage_tets[tet_count].ID != -1) && (!tet_boundary_check(tet_count))){
			Tetrahedron	tet;
			t_id_now = ResultT.size();

			tet.ID = t_id_now;
			tet.set_label(m_storage_tets[tet_count].get_label());

			bool all_right = true;

			for(int index0=0; index0<4; index0++){
				int vert_index = m_storage_tets[tet_count].Verts[index0];
				it = vert_window_map.find(vert_index);
				if(it != vert_window_map.end()){
						tet.Verts[index0] = it->second;
				}else{
						tet.Verts[index0] = -1;
						all_right = false;
						throw exception( "wrong here, should all be in tets" );
				}								 
			}

			if(all_right){
				ResultT.push_back(tet);				
				ResultT_tmp.push_back(tet);
				m_storage_tets[tet_count].invalidate(this);
			}
		}
	}

	//std::cout << "set result v and result t cost time: " << TimerObj.tick() << std::endl;

	int ntets = ResultT_tmp.size();
	for(int i=0; i<ntets; ++i) 
	{
		Tetrahedron* t = &ResultT_tmp[i];
		bool labeled = false;
		for(int j=0; j<4; ++j) 
		{
			int vi = t->Verts[j];
			ResultV_tmp[vi].Tets.push_back(t->ID);
		}
	}

	//std::cout << "correlate verts and tets: " << TimerObj.tick() << std::endl;
	//RenderResultVerts(vis);
	//RenderResults(vis);

}


void cwg::TetrahedralMesh::getResultV_T(){
	vert_map.clear();

	int win_inside_index =0;
	int win_x = 0, win_y = 0, win_z = 0;
	std::ifstream vertFile; //std::ofstream vertFileClean;
	std::ifstream tetFile;  //std::ofstream tetFileClean;
	
	std::string vertFileName0 = get_foldername() + "Data\\vert_";
	std::string tetFileName0  = get_foldername() + "Data\\tet_";
	std::string vertFileName;
	std::string tetFileName;

	int num_x = ceil((float)m_X / (float)(slide_X/2));
	int num_y = ceil((float)m_Y / (float)(slide_Y/2));
	int num_z = ceil((float)m_Z / (float)(slide_Z/2));

	int v_id_now = ResultV.size();	
	TetVertex vert;
	
	std::map<int, int>::iterator it;

	for(int k=0; k<num_z; k++){
		for(int i=0; i<num_x; i++){
			for(int j=0; j<num_y; j++){
				std::string tmpfilename("");
				tmpfilename += std::to_string(i);
				tmpfilename += "_";
				tmpfilename += std::to_string(j);
				tmpfilename += "_";
				tmpfilename += std::to_string(k);
				tmpfilename += ".txt";

				vertFileName = vertFileName0 + tmpfilename;
				vertFile.open(vertFileName, ios::binary);

				int x_start = (slide_X/2) * i;
				int y_start = (slide_Y/2) * j;
				int z_start = (slide_Z/2) * k;

				if(vertFile.is_open()){
					vertFile.seekg(-(sizeof(int)), vertFile.end);
					int vert_size;
					vertFile.read((char*)&vert_size, sizeof(int));
					vertFile.seekg(0, vertFile.beg);

					for(int vert_count=0; vert_count<vert_size; vert_count++){
						//idx idy idz
						int vert_i, vert_j, vert_k;
						vertFile.read((char*)&vert_i, sizeof(int));
						vertFile.read((char*)&vert_j, sizeof(int));
						vertFile.read((char*)&vert_k, sizeof(int));

						int idx = (z_start+vert_k)*(m_X+1)*(m_Y+1) + (x_start+vert_i)*(m_Y+1) + (y_start+vert_j);
						v_id_now = ResultV.size();

						it = vert_map.find(idx);
						if(it == vert_map.end()){
							vert_map.insert( std::pair<int, int>(idx, v_id_now) );
							//Make the vert and push it into the resultV!
							vert.ID = v_id_now;

							int keep_number;
							vertFile.read((char*)&keep_number, sizeof(int));

							//pos.x pos.y pos.z
							double vert_poses[3];
							vertFile.read((char*)vert_poses, sizeof(double)*3);
							vert.pos() = Cleaver::vec3(vert_poses[0], vert_poses[1], vert_poses[2]);

							//labels
							int label_size;
							vertFile.read((char*)&label_size, sizeof(int));
							int label_int;
							for(int label_index=0; label_index<label_size; label_index++){
								vertFile.read((char*)&label_int, sizeof(int));
								vert.assign_label(label_int);
							}

							ResultV.push_back(vert);
						}else{
							int keep_number;
							vertFile.read((char*)&keep_number, sizeof(int));
							double vert_poses[3];
							vertFile.read((char*)vert_poses, sizeof(double)*3);//x y z (pos)
							int label_size;
							vertFile.read((char*)&label_size, sizeof(int));
							int label_int;
							for(int label_index=0; label_index<label_size; label_index++){
								vertFile.read((char*)&label_int, sizeof(int));
							}//label
						}
					}
				}else{
					std::cout << "vert file can not open! so it can not be read!" << std::endl;
				}

				vertFile.close();

			}
		}	
	}
	
	std::cout << "Finished loading Verts!" << std::endl;
	std::cout << "size: " << vert_map.size() << std::endl;

	//Tet files
	int t_id_now = ResultT.size();
	Tetrahedron	tet;
	
	for(int k=0; k<num_z; k++){
		for(int i=0; i<num_x; i++){
			for(int j=0; j<num_y; j++){	
				std::string tmpfilename("");
				tmpfilename += std::to_string(i);
				tmpfilename += "_";
				tmpfilename += std::to_string(j);
				tmpfilename += "_";
				tmpfilename += std::to_string(k);
				tmpfilename += ".txt";
	
				tetFileName = tetFileName0 + tmpfilename;
				tetFile.open(tetFileName, ios::binary);

				int x_start = (slide_X/2) * i;
				int y_start = (slide_Y/2) * j;
				int z_start = (slide_Z/2) * k;

				if(tetFile.is_open()){
					tetFile.seekg(-(sizeof(int)), tetFile.end);
					int tet_size;
					tetFile.read((char*)&tet_size, sizeof(int));
					tetFile.seekg(0, tetFile.beg);

					for(int tet_count=0; tet_count<tet_size; tet_count++){
			
						//get the tet index
						int x, y, z, index0;
						tetFile.read((char*)&x, sizeof(int));
						tetFile.read((char*)&y, sizeof(int));
						tetFile.read((char*)&z, sizeof(int));
						tetFile.read((char*)&index0, sizeof(int));

						//tet_idx tet_idy tet_idz index0
						t_id_now = ResultT.size();
						tet.ID = t_id_now;

						//bool all_inside_cube, useless here
						int all_inside_cube_int;
						tetFile.read((char*)&all_inside_cube_int, sizeof(int));
					
						//label

						int label_int;
						tetFile.read((char*)&label_int, sizeof(int));
						tet.set_label(label_int);

						//verts
						int vert_idx_int, vert_idy_int, vert_idz_int;
						bool inside_cube_int;

						bool all_right = true;
						for(int index0=0; index0<4; index0++){
							tetFile.read((char*)&inside_cube_int, sizeof(bool));
							tetFile.read((char*)&vert_idx_int, sizeof(int));
							tetFile.read((char*)&vert_idy_int, sizeof(int));
							tetFile.read((char*)&vert_idz_int, sizeof(int));

							int vert_index = (vert_idz_int+z_start)*(m_X+1)*(m_Y+1) + (vert_idx_int+x_start)*(m_Y+1) + (vert_idy_int+y_start);
							it = vert_map.find(vert_index);
							if(it != vert_map.end()){
								tet.Verts[index0] = it->second;
							}else{
								tet.Verts[index0] = -1;
								all_right = false;
								throw exception("The vert is not in the map! something is wrong!");
							}								 
						}

						if(all_right)
							ResultT.push_back(tet);
					}
				}else{
					std::cout << "tet file can not open! so it can not be read!" << std::endl;
				}

				tetFile.close();
				//tetFileClean.open(tetFileName, ios::binary);
				//tetFileClean.close();
			}
		}
	}

}


void cwg::TetrahedralMesh::compute_component_nums(){
	
	int face_list[4][3] =  {{0, 1, 2},
							{1, 2, 3},
							{0, 1, 3},
							{0, 2, 3}};
	
	std::vector<int> labels_used;
	
	int tet_size = m_storage_tets.size();
	labels_used.push_back(0);
	
	for(int i=0; i<tet_size; i++){
		if(m_storage_tets[i].ID != -1){
			int label_this_tet = m_storage_tets[i].get_label();
			if(std::find(labels_used.begin(), labels_used.end(), label_this_tet) == labels_used.end()){
				labels_used.push_back(label_this_tet);
			}
		}	
	}
	//labels_used.erase(labels_used.begin());
	

	//labels_used.push_back(30825);

	for(int i=0; i<labels_used.size(); i++){
		std::vector<int> tets_component_counting;

		for(int j=0; j<tet_size; j++){
			if ( m_storage_tets[j].ID != -1 && m_storage_tets[j].get_label() == labels_used[i]){
				tets_component_counting.push_back(j);
			}
		}
		
		/*if(tets_component_counting.size() == compute_component_branch_size){
			continue;
		}else{
			compute_component_branch_size = tets_component_counting.size();
			if(compute_component_branch_size > 26){
				continue;
			}
		}*/

		int branch_num=0;
		map<int, bool> tets_checked;
		for(int j=0; j<tets_component_counting.size(); j++){
			tets_checked.insert(make_pair(tets_component_counting[j], false));
		}

		for(int j=0; j<tets_component_counting.size(); j++){
			if(!tets_checked[tets_component_counting[j]]){
				std::vector<int> branch_this;
				branch_this.push_back(tets_component_counting[j]);

				while(branch_this.size() > 0){
					int tet_check_now = branch_this[0];

					for(int face_index=0; face_index<4; face_index++){
						std::vector<int> shared_tets;
						std::vector<int> shared_tets_2;

						std::sort(m_storage_verts[m_storage_tets[tet_check_now].Verts[face_list[face_index][0]]].Tets.begin(), 
								  m_storage_verts[m_storage_tets[tet_check_now].Verts[face_list[face_index][0]]].Tets.end());
						
						std::sort(m_storage_verts[m_storage_tets[tet_check_now].Verts[face_list[face_index][1]]].Tets.begin(), 
								  m_storage_verts[m_storage_tets[tet_check_now].Verts[face_list[face_index][1]]].Tets.end());

						
						std::sort(m_storage_verts[m_storage_tets[tet_check_now].Verts[face_list[face_index][2]]].Tets.begin(), 
								  m_storage_verts[m_storage_tets[tet_check_now].Verts[face_list[face_index][2]]].Tets.end());

						std::set_intersection (m_storage_verts[m_storage_tets[tet_check_now].Verts[face_list[face_index][0]]].Tets.begin(), 
											   m_storage_verts[m_storage_tets[tet_check_now].Verts[face_list[face_index][0]]].Tets.end(),
												m_storage_verts[m_storage_tets[tet_check_now].Verts[face_list[face_index][1]]].Tets.begin(),
												m_storage_verts[m_storage_tets[tet_check_now].Verts[face_list[face_index][1]]].Tets.end(), 
														  back_inserter(shared_tets));

						std::set_intersection (shared_tets.begin(), shared_tets.end(), 
											   m_storage_verts[m_storage_tets[tet_check_now].Verts[face_list[face_index][2]]].Tets.begin(),
											   m_storage_verts[m_storage_tets[tet_check_now].Verts[face_list[face_index][2]]].Tets.end(), 
												back_inserter(shared_tets_2));

						for(int shared_tet_index=0; shared_tet_index<shared_tets_2.size(); shared_tet_index++){
							if(std::find(branch_this.begin(), branch_this.end(), shared_tets_2[shared_tet_index]) == branch_this.end()){
								if(m_storage_tets[shared_tets_2[shared_tet_index]].get_label() == labels_used[i]){
									if(!tets_checked[shared_tets_2[shared_tet_index]]){
										branch_this.push_back(shared_tets_2[shared_tet_index]);
									}
								}
								
							}
						}
						
					}

					tets_checked[tet_check_now] = true;
					branch_this.erase(branch_this.begin());
				}

				branch_num ++;
			}

		}

		std::cout << "branch label: "<< labels_used[i] << " branch num: " << branch_num << " tets num: " << tets_component_counting.size() << std::endl;

		if(branch_num > 1){
			int check_num;
			std::cout << "got the one that fails! " << labels_used[i] << std::endl;
			std::cin >> check_num;
		}

		storage_to_result();
		stringstream vert_sstream;
		vert_sstream << "after_improve_vert_" << labels_used[i] << ".txt";
		string vert_file_name = vert_sstream.str();

		stringstream tet_sstream;
		tet_sstream << "after_improve_tet_" << labels_used[i] << ".txt";
		string tet_file_name = tet_sstream.str();

		save_result_v_t_label(vert_file_name, tet_file_name, labels_used[i]);

		/*for(int tet_idx=0; tet_idx < ResultT.size(); tet_idx++){
			if(ResultT[tet_idx].get_label() == labels_used[i]){
				std::cout << "start outputing this tet! " << tet_idx << std::endl;
				write_obj_for_debug(tet_idx);
			}
		}*/
	}
	
}

void cwg::TetrahedralMesh::compute_component_nums_result(){
	int face_list[4][3] =  {{0, 1, 2},
							{1, 2, 3},
							{0, 1, 3},
							{0, 2, 3}};


	std::vector<int> labels_used;
	labels_used.push_back(0);
	int tet_size = ResultT.size();
	for(int i=0; i<tet_size; i++){
		int label_this_tet = ResultT[i].get_label();
		if(std::find(labels_used.begin(), labels_used.end(), label_this_tet) == labels_used.end()){
				labels_used.push_back(label_this_tet);
		}
	}

	labels_used.erase(labels_used.begin());
	for(int i=0; i<labels_used.size(); i++){
		std::vector<int> tets_component_counting;

		for(int j=0; j<tet_size; j++){
			if ( ResultT[j].get_label() == labels_used[i]){
				tets_component_counting.push_back(j);
			}
		}

		int branch_num=0;
		map<int, bool> tets_checked;
		for(int j=0; j<tets_component_counting.size(); j++){
			tets_checked.insert(make_pair(tets_component_counting[j], false));
		}

		int pushed_back_item = 0;
		for(int j=0; j<tets_component_counting.size(); j++){
			if(!tets_checked[tets_component_counting[j]]){
				std::vector<int> branch_this;
				branch_this.push_back(tets_component_counting[j]);

				while(branch_this.size() > 0){
					int tet_check_now = branch_this[0];

					for(int face_index=0; face_index<4; face_index++){
						std::vector<int> shared_tets;
						std::vector<int> shared_tets_2;

						std::sort(ResultV[ResultT[tet_check_now].Verts[face_list[face_index][0]]].Tets.begin(), 
								  ResultV[ResultT[tet_check_now].Verts[face_list[face_index][0]]].Tets.end());
						
						std::sort(ResultV[ResultT[tet_check_now].Verts[face_list[face_index][1]]].Tets.begin(), 
								  ResultV[ResultT[tet_check_now].Verts[face_list[face_index][1]]].Tets.end());

						
						std::sort(ResultV[ResultT[tet_check_now].Verts[face_list[face_index][2]]].Tets.begin(), 
								  ResultV[ResultT[tet_check_now].Verts[face_list[face_index][2]]].Tets.end());

						std::set_intersection (ResultV[ResultT[tet_check_now].Verts[face_list[face_index][0]]].Tets.begin(), 
											   ResultV[ResultT[tet_check_now].Verts[face_list[face_index][0]]].Tets.end(),
												ResultV[ResultT[tet_check_now].Verts[face_list[face_index][1]]].Tets.begin(),
												ResultV[ResultT[tet_check_now].Verts[face_list[face_index][1]]].Tets.end(), 
														  back_inserter(shared_tets));

						std::set_intersection (shared_tets.begin(), shared_tets.end(), 
											   ResultV[ResultT[tet_check_now].Verts[face_list[face_index][2]]].Tets.begin(),
											   ResultV[ResultT[tet_check_now].Verts[face_list[face_index][2]]].Tets.end(), 
												back_inserter(shared_tets_2));
						
						for(int shared_tet_index=0; shared_tet_index<shared_tets_2.size(); shared_tet_index++){
							if(std::find(branch_this.begin(), branch_this.end(), shared_tets_2[shared_tet_index]) == branch_this.end()){
								if(ResultT[shared_tets_2[shared_tet_index]].get_label() == labels_used[i]){
									if(!tets_checked[shared_tets_2[shared_tet_index]]){
										pushed_back_item ++ ;
										branch_this.push_back(shared_tets_2[shared_tet_index]);
									}
								}
							}
						}
						
					}

					tets_checked[tet_check_now] = true;
					branch_this.erase(branch_this.begin());
				}

				branch_num ++;
			}

		}

		std::cout << "tets_component_counting size: " << tets_component_counting.size() << std::endl;
		std::cout << "labels this: "<< labels_used[i] << " num of components: " << branch_num << std::endl;
		std::cout << "pushed back item num: " << pushed_back_item << std::endl;

	}
	
}

void cwg::TetrahedralMesh::storage_to_result(){
	vert_map.clear();
	ResultV.clear();
	ResultT.clear();

	std::map<int, int>::iterator it;

	for(int i=0; i<m_storage_verts.size(); i++){
		if(m_storage_verts[i].ID != -1){
			TetVertex vert;
			vert.ID = ResultV.size();
			vert.pos() = Cleaver::vec3(m_storage_verts[i].pos().x, m_storage_verts[i].pos().y, m_storage_verts[i].pos().z);
			for(int j=0; j<m_storage_verts[i].m_labels.size(); j++){
				vert.assign_label(m_storage_verts[i].m_labels[j]);
			}
			vert_map.insert(make_pair(i, ResultV.size()));
			ResultV.push_back(vert);
		}
	}

	for(int i=0; i<m_storage_tets.size(); i++){
		if(m_storage_tets[i].ID != -1){
			Tetrahedron	tet;
			tet.ID = ResultT.size();
			tet.Verts[0] = vert_map[m_storage_tets[i].Verts[0]];
			tet.Verts[1] = vert_map[m_storage_tets[i].Verts[1]];
			tet.Verts[2] = vert_map[m_storage_tets[i].Verts[2]];
			tet.Verts[3] = vert_map[m_storage_tets[i].Verts[3]];
			ResultV[tet.Verts[0]].correlate_tet(tet.ID);
			ResultV[tet.Verts[1]].correlate_tet(tet.ID);
			ResultV[tet.Verts[2]].correlate_tet(tet.ID);
			ResultV[tet.Verts[3]].correlate_tet(tet.ID);


			tet.set_label(m_storage_tets[i].get_label());

			ResultT.push_back(tet);
		}
	}


}



void cwg::TetrahedralMesh::remove_border_tets()
{
	vector<Cleaver::Tet*>& tets = m_cleaver_tet_mesh->tets;

	for (int i=0; i<tets.size(); i++)
	{
		Cleaver::Tet* t = tets[i];
		if (t->mat_label == m_labelvolume->nTsp()-1)
		{
			for (int vi=0; vi<4; vi++)
			{
				bool found;
				Cleaver::Vertex3D* v = tets[i]->verts[vi];
				int idx = find_tet_storage_index(v, t, found);	
				if (found)
					v->tets.erase(v->tets.begin()+idx);
			}
			SAFE_DELETE(tets[i]);
		}
	}
}


void cwg::TetrahedralMesh::flatten_boundary_tets()
{
	double lambda = m_labelvolume->nPaddings;
	vec3 lower = vec3(lambda, lambda, lambda);
	vec3 upper = vec3(m_labelvolume->Width()-1.0-lambda, 
		m_labelvolume->Height()-1.0-lambda, m_labelvolume->nFrames()-1.0-lambda);
	BoundingBox BoundingBox_Exact(lower, upper-lower);

	vector<Vertex3D*>& verts = m_cleaver_tet_mesh->verts;

	for (int i=0; i<verts.size(); i++)
	{
		Vertex3D* v = verts[i];
		BoundingBox_Exact.round(v->pos());
	}
}


void cwg::TetrahedralMesh::remove_degenerated_tets()
{
	vector<Cleaver::Tet*>& tets = m_cleaver_tet_mesh->tets;
	for (int i=0; i<tets.size(); i++)
	{
		Tet* t = tets[i];
		if (t)
		{
			double vol = tet_signed_volume(t);

			if (abs(vol) < 1e-5)
			{
				for (int vi=0; vi<4; vi++)
				{
					bool found;
					Vertex3D* v = tets[i]->verts[vi];
					int idx = find_tet_storage_index(v, t, found);	
					if (found)
						v->tets.erase(v->tets.begin()+idx);
				}
				SAFE_DELETE(tets[i]);
			}
		}
	}
}


void cwg::TetrahedralMesh::remove_nonreferenced_verts()
{
	vector<Cleaver::Tet*>& tets = m_cleaver_tet_mesh->tets;
	vector<Vertex3D*>& verts = m_cleaver_tet_mesh->verts;

	vector<int> verts_reference_table(verts.size(), 0);
	for (int i=0; i<tets.size(); i++)
	{
		if (!tets[i])
			continue;

		for (int j=0; j<4; j++)
		{
			Vertex3D* vij = tets[i]->verts[j];
			verts_reference_table[vij->tm_v_index]++;
		}
	}

	for (int i=0; i<verts.size(); i++)
	{
		if (verts_reference_table[i] == 0)
			SAFE_DELETE(verts[i]);
	}
}


void cwg::TetrahedralMesh::translate_model()
{
	vector<Cleaver::Vertex3D*>& verts = m_cleaver_tet_mesh->verts;
	vec3 tmp(m_labelvolume->nPaddings, m_labelvolume->nPaddings, m_labelvolume->nPaddings);
	//vec3 tmp(2,2,2);
	for (int i=0; i<verts.size(); i++)
	{
		if (!verts[i])
			continue;

		verts[i]->pos() = verts[i]->pos() - tmp;
	}
}


void cwg::TetrahedralMesh::initial_sortout()
{
	vector<Cleaver::Tet*>& tets = m_cleaver_tet_mesh->tets;
	vector<Cleaver::Vertex3D*>& verts = m_cleaver_tet_mesh->verts;

	int nverts(0), ntets(0);
	for (int i=0; i<verts.size(); i++)
	{
		if (verts[i])
			nverts++;
	}
	for (int i=0; i<tets.size(); i++)
	{
		if (tets[i])
			ntets++;
	}


	m_storage_verts.resize(nverts);
	//m_storage_verts = new TetVertex [nverts];

	//m_storage_tets = new Tetrahedron [ntets];
	m_storage_tets.resize(ntets);

	int vi(0), ti(0);
	map<int, int> IndexMap;
	for (int i=0; i<verts.size(); i++)
	{
		if (verts[i])
		{
			m_storage_verts[vi].copy_info_from_cleaver_vert(verts[i]);
			const TetVertex* v = &m_storage_verts[vi];
			IndexMap.insert( make_pair(verts[i]->tm_v_index, v->ID) );
			vi++;
		}
	}
	for (int i=0; i<tets.size(); i++)
	{
		if (tets[i])
		{
			Tetrahedron* t = &m_storage_tets[ti];
			t->set_label(tets[i]->mat_label);
			for (int j=0; j<4; j++)
			{
				int origin_id = tets[i]->verts[j]->tm_v_index;
				int now_id = IndexMap.at(origin_id);
				t->Verts[j] = now_id;
			}

			m_storage_tets[ti] = *t;
			ti++;
		}
	}

	for (int i=0; i<m_storage_tets.size(); i++)
		m_storage_tets[i].correlate_verts(this);
}


void cwg::TetrahedralMesh::build_edges_wenhua(){
	int e_onBoundary_count = 0;

	EdgeSet.clear();
	EdgeHeap.clear();

	int tet_size = m_storage_tets.size();


	int ycount = 0;
	for(int k=0; k<slide_Z; k++){
		for(int i=0; i<slide_X; i++){
			for(int j=0; j<slide_Y; j++){
				int tet_index0 = (k*slide_X*slide_Y + i*slide_Y + j)*tetn;

				for(int index=0; index<tetn; index++){
					int ti = tet_index0 +index;

					if(m_storage_tets[ti].ID != -1){
						Tetrahedron* t = &m_storage_tets[ti];

						for (int i0=0; i0<4; i0++){
							int vi = t->Verts[i0];

							for (int j0=i0+1; j0<4; j0++){
								int vj = t->Verts[j0];
								if((m_storage_verts[vi].keep < 0) || (m_storage_verts[vj].keep < 0))
									std::cout << "impossible!!! keep should never be smaller than 0!" << std::endl;

								if(!m_storage_verts[vi].boundary && !m_storage_verts[vj].boundary && (m_storage_verts[vi].keep==0) && (m_storage_verts[vj].keep==0) )
									if(!m_storage_verts[vi].slideBoundary && !m_storage_verts[vj].slideBoundary)
										EdgeSet.insert( cwg::ordered_pair(vi,vj));
							}
						}
					}
				}

			}
		}
	}

	int count = 0;	
	EdgeHeap.resize(EdgeSet.size());

	//Real calculation, initialize the matrix
	for (set<cwg::ordered_pair>::iterator it = EdgeSet.begin(); it != EdgeSet.end(); it++)
	{

		TetEdge* e = &m_storage_edges[count];
		e->onBoundary = true;
		e->ID = count;
		e->set_verts(it->First(), it->Second());
		AdjMatrix->set_entry(it->First(), it->Second(), e->ID+0.5);
		EdgeHeap.set_edge(e, count);
		++count;
	}

}


void cwg::TetrahedralMesh::release()
{
	m_storage_verts.clear();//SAFE_DELETE_ARRAY(m_storage_verts);
	m_storage_verts.shrink_to_fit();
	m_storage_tets.clear();//SAFE_DELETE_ARRAY(m_storage_tets);
	m_storage_tets.shrink_to_fit();
	m_storage_edges.clear();//SAFE_DELETE_ARRAY(m_storage_edges);
	m_storage_edges.shrink_to_fit();

	EdgeHeap.clear();
	delete_vector_elems(Tris);
	SAFE_DELETE(AdjMatrix);
}


void cwg::TetrahedralMesh::compute_verts_border_code()
{
	double lambda = 1e-5;
	vec3 lower = vec3(lambda, lambda, lambda);
	vec3 upper = vec3(m_X  - lambda , m_Y  - lambda , m_Z - lambda );
	BoundingBox BoundingBox_Smaller(lower, upper-lower);
	double BTH = lower.x;
	double W_BTH = upper.x;
	double H_BTH = upper.y;
	double F_BTH = upper.z;
	unsigned char vert_type_tmp, dimcode;
	vec3 vpos;
	for (int i=0; i<m_storage_verts.size(); i++)
	{
		vert_type_tmp = 0;
		vpos = m_storage_verts[i].pos();
		// x-dim:
		dimcode = vpos.x < BTH ? 1 : (vpos.x > W_BTH ? 2 : 0);
		vert_type_tmp |= ( dimcode << 0 );
		vert_type_tmp += ( dimcode != 0 ? (unsigned char)Cleaver::VERTEX_DLT : 0 );
		// y-dim:
		dimcode = vpos.y < BTH ? 1 : (vpos.y > H_BTH ? 2 : 0);
		vert_type_tmp |= ( dimcode << 2 );
		vert_type_tmp += ( dimcode != 0 ? (unsigned char)Cleaver::VERTEX_DLT : 0 );

		// z-dim: 
		dimcode = vpos.z < BTH ? 1 : (vpos.z > F_BTH ? 2 : 0);
		vert_type_tmp |= ( dimcode << 4 );
		vert_type_tmp += ( dimcode != 0 ? (unsigned char)Cleaver::VERTEX_DLT : 0 );

		m_storage_verts[i].set_full_border_code(vert_type_tmp);
	}
}


void cwg::TetrahedralMesh::compute_verts_Q_wenhua(TetrahedralMesh* tetmesh)
{

	int v_size = m_storage_verts.size();

	for (int i = 0; i<v_size; i++){

		if(m_storage_verts[i].ID != -1)
			m_storage_verts[i].ComputeInitialQ(tetmesh);
	}

}


void cwg::TetrahedralMesh::compute_weight_from_edge_variation()
{
	for (int i = 0 ; i < ResultV.size(); ++i )
	{
		double edge_avg = 0;
		for (int j = 0; j < ResultV[i].Tets.size(); ++j)
		{
			edge_avg += ResultT[ResultV[i].Tets[j]].avg_edge_length_sliding_result(this);

		}
		ResultV[i].setAvgWeight((1.0/ ((edge_avg/ResultV[i].Tets.size())*(edge_avg/ResultV[i].Tets.size())) )*100 );

	}
}


void cwg::TetrahedralMesh::densityField()
{
	cout << "set weight ...\t";
	TimerObj.tick();

	std::queue<int> surfV;



	for (int i=0; i<m_storage_verts.size(); ++i)
	{
		if(m_storage_verts[i].NLabels()>1)
		{
			surfV.push(i);
			m_storage_verts[i].weightLevel = 1;
		}
	}

	compute_lfs(surfV,m_storage_verts.size());

	for (int i=0; i<m_storage_verts.size(); ++i)
	{

		m_storage_verts[i].setAvgWeight(1.0/(m_storage_verts[i].weightLevel * m_storage_verts[i].weightLevel) * 100.0);
	}

	cout << "cost " << TimerObj.tick() << " sec." << endl;
}


void cwg::TetrahedralMesh::compute_lfs(std::queue<int>& surfV , int vsize)
{

	while(!surfV.empty())
	{
		int i = surfV.front();
		surfV.pop();

		if(i-1 >=0 && m_storage_verts[i-1].weightLevel == 0)
		{
			m_storage_verts[i-1].weightLevel = m_storage_verts[i].weightLevel+1;
			surfV.push(i-1);
		}
		if(i+1 < vsize && m_storage_verts[i+1].weightLevel == 0)
		{
			m_storage_verts[i+1].weightLevel = m_storage_verts[i].weightLevel+1;
			surfV.push(i+1);
		}
		if(i-(slide_Y+1) >=0 && m_storage_verts[i-(slide_Y+1)].weightLevel == 0)
		{
			m_storage_verts[i-(slide_Y+1)].weightLevel = m_storage_verts[i].weightLevel+1;
			surfV.push(i-(slide_Y+1));
		}
		if(i+(slide_Y+1) < vsize && m_storage_verts[i+(slide_Y+1)].weightLevel == 0)
		{
			m_storage_verts[i+(slide_Y+1)].weightLevel = m_storage_verts[i].weightLevel+1;
			surfV.push(i+(slide_Y+1));
		}
		if(i-(slide_X+1)*(slide_Y+1) >=0 && m_storage_verts[i-(slide_X+1)*(slide_Y+1)].weightLevel == 0)
		{
			m_storage_verts[i-(slide_X+1)*(slide_Y+1)].weightLevel = m_storage_verts[i].weightLevel+1;
			surfV.push(i-(slide_X+1)*(slide_Y+1));
		}
		if(i+(slide_X+1)*(slide_Y+1) < vsize && m_storage_verts[i+(slide_X+1)*(slide_Y+1)].weightLevel == 0)
		{
			m_storage_verts[i+(slide_X+1)*(slide_Y+1)].weightLevel = m_storage_verts[i].weightLevel+1;
			surfV.push(i+(slide_X+1)*(slide_Y+1));
		}


	}


}


void cwg::TetrahedralMesh::compute_edge_contraction()
{
	//cout << "compute edge contraction ...\t";
	for (int i= 0; i<EdgeHeap.size(); i++)
	{
		TetEdge* e = EdgeHeap.get_edge(i);
		e->ComputeContraction(this);
	}
}


void cwg::TetrahedralMesh::build_heap()
{
	EdgeHeap.random_shuffle();
}


void cwg::TetrahedralMesh::simplify( int tar_nverts )
{
	int itime(0);

	bool heap_to_bottom = false;
	int searched_edges = 0;
	double max_cost = 1000;
	int max_searched_edges = 0;

	for (itime=0; ; itime++)
	{	
		TetEdge* e = EdgeHeap.top_edge();

		double topcost = e->cost();
		if (e->ID == -1 || (searched_edges > (EdgeHeap.size()+1000)))
		{
			heap_to_bottom = true;
			break;
		}

		e->compute_around_info(this);
		e->ComputeContraction(this);

		if(e->able_to_collapse){

			if(e->Cost > max_cost)
				max_cost = e->Cost;

			EdgeHeap.pop();
			contract(e);

			searched_edges = 0;
		}else{
			max_searched_edges ++;

			e->Cost = e->Cost * (max_searched_edges + 1);

			EdgeHeap.update(e->id());
			searched_edges +=  1;
		}	
	
	}
}


void cwg::TetrahedralMesh::contract( TetEdge* e )
{
	int v0 = e->Verts[0];
	int v1 = e->Verts[1];

	assert(m_storage_verts[v0].border_type() >= m_storage_verts[v1].border_type());

	m_storage_verts[v0].pos() = m_storage_verts[v1].pos() = e->ContractPos;
	m_storage_verts[v0].update_funcval();

	const TetEdgeAroundInfo& e_around_info = e->get_around_info();
	const vector<int>& v0_private_tets = e_around_info.v0_private_tets;
	const vector<int>& v1_private_tets = e_around_info.v1_private_tets;
	const vector<int>& overlap_tets_v0 = e_around_info.v0_private_overlap_tets;
	const vector<int>& overlap_tets_v1 = e_around_info.v1_private_overlap_tets;
	const vector<int>& v0_private_nonoverlap_tets = e_around_info.v0_private_nonoverlap_tets;
	const vector<int>& v1_private_nonoverlap_tets = e_around_info.v1_private_nonoverlap_tets;
	const vector<int>& edge_tets = e_around_info.edge_tets;
	const set<int>& affected_verts = e_around_info.affected_verts;

	remove_tets(edge_tets);
	remove_edge(&m_storage_verts[v0], &m_storage_verts[v1], e);

	replace_verts_in_v1_tets(&m_storage_verts[v0], &m_storage_verts[v1]);
	replace_verts_in_v1_edges(&m_storage_verts[v0], &m_storage_verts[v1]);;

	//********************************************CHANGED HERE**********************************************
	//*************************** USED TO DELETE OVERLAP TETS V0 AND OVERLAP TETS V1 ***********************
	//*************************** CHANGED TO DELETE OVERLAP TETS V1 ONLY************************************
	
	m_storage_verts[v0].UpdateQ(m_storage_verts[v1].QMatrix());
	update_edges_of_vert(&m_storage_verts[v0]);

	m_storage_verts[v0].update_labels(this);

	m_storage_verts[v1].invalidate(); 
	
}


void cwg::TetrahedralMesh::compute_edge_tets( const TetVertex* v0, const TetVertex* v1, std::vector<int>& edge_tets ) const
{
	const vector<int>& v0_tets = v0->Tets;
	const vector<int>& v1_tets = v1->Tets;
	edge_tets.clear();
	std::set_intersection(v0_tets.begin(), v0_tets.end(), v1_tets.begin(), v1_tets.end(), 
		inserter(edge_tets, edge_tets.end()));
}


void cwg::TetrahedralMesh::compute_edge_tets( TetEdge* e, std::vector<int>& edge_tets ) const
{	
	compute_edge_tets(&m_storage_verts[ e->Verts[0] ], &m_storage_verts[ e->Verts[1] ], edge_tets);
}


void cwg::TetrahedralMesh::remove_tets( const std::vector<int>& edge_tets )
{
	for (vector<int>::const_iterator it = edge_tets.begin(); it != edge_tets.end(); it++)
	{
		Tetrahedron* t = &this->m_storage_tets[*it];

		for (int i=0; i<4; i++){
			m_storage_verts[ t->Verts[i] ].decorrelate_tet(t->ID);
		}
	
		t->ID = -1;
	}
}


void cwg::TetrahedralMesh::remove_edge( TetVertex* v0, TetVertex* v1, TetEdge* e )
{
	bool flag;
	int eij = static_cast<int>(AdjMatrix->get_entry(v0->ID, v1->ID, flag));

	if (flag)
	{
		if (eij < 0)
		{
			cout << e->ID << "\t" << e->Verts[0] << "\t" << e->Verts[1] << endl;
			cout << e->Cost << endl;
			throw exception("void cwg::TetrahedronMesh::remove_edge( TetVertex* v0, TetVertex* v1, TetEdge* e ): \
							eij < 0: IMPOSSIBLE!");
		}
		else
		{
			AdjMatrix->set_entry(v0->ID, v1->ID, -1.0);
		}
	}
	else // meaning there's no this element
	{
		throw exception("void cwg::TetrahedronMesh::remove_edge( TetVertex* v0, TetVertex* v1, TetEdge* e ): \
						there is no this edge stored in the AdjMatrix: IMPOSSIBLE!");
	}

	e->ID = -1;
	e->onBoundary = false;
	//m_invalid_num_edges++;
}


//remove the edges which are determined to be removed here.
//The numbers are calculated before and do the remove now.
void cwg::TetrahedralMesh::remove_edge_sliding( TetVertex* v0 )
{

	std::vector< std::vector<Sparse_Entry> >& entryset = AdjMatrix->get_entryset();
	std::vector<Sparse_Entry>& v0_edges = entryset[v0->ID];

	for (int j=0; j<v0_edges.size(); j++)
	{
		int vj = v0_edges[j].index; // v1's neighboring vert: vj
		double tmpval = v0_edges[j].value; 
		int eij = tmpval;

		if (tmpval >= -1e-6)
		{
			AdjMatrix->set_entry(v0->ID, vj, -1.0);
			EdgeHeap[eij]->ID = -1;
			EdgeHeap[eij]->onBoundary = false;
		}
	}



}

void cwg::TetrahedralMesh::remove_edge_sliding( TetVertex* v0, TetVertex* v1 )
{

	int ei = AdjMatrix->get_entry(v0->ID, v1->ID);
	if (ei != -1)
	{
		EdgeHeap[ei]->ID = -1;
		AdjMatrix->set_entry(v0->ID, v1->ID, -1.0);
	}



}

void cwg::TetrahedralMesh::replace_verts_in_v1_tets( TetVertex* v0, TetVertex* v1 )
{
	vector<int> v1_tets = v1->Tets;
	for (vector<int>::iterator it = v1_tets.begin(); it != v1_tets.end(); it++)
	{
		Tetrahedron* t = &this->m_storage_tets[ *it ];

		int vtid = t->find_vert(v1->ID);
		t->Verts[ vtid ] = v0->ID;
		v0->correlate_tet(t->ID);
		v1->decorrelate_tet(t->ID);
	}
}

void cwg::TetrahedralMesh::replace_verts_in_v1_edges( TetVertex* v0, TetVertex* v1)
{
	std::vector< std::vector<Sparse_Entry> >& entryset = AdjMatrix->get_entryset();
	std::vector<Sparse_Entry>& v1_edges = entryset[v1->ID];
	for (int j=0; j<v1_edges.size(); j++)
	{
		int vj = v1_edges[j].index; // v1's neighboring vert: vj
		double tmpval = v1_edges[j].value; // e1j: v1---vj
		int e1j = tmpval;

		if (tmpval >= -1e-6)
		{
			AdjMatrix->set_entry(v1->ID, vj, -1.0);

			bool flag;
			double e0j = AdjMatrix->get_entry(v0->ID, vj, flag);
			
			if (flag){ // e0j already exists, remove e1j
				EdgeHeap[e1j]->Cost = std::numeric_limits<double>::max();
				EdgeHeap.remove(EdgeHeap[e1j]->id());
			}
			else{
				AdjMatrix->set_entry(v0->ID, vj, tmpval);

				assert(EdgeHeap[e1j]->id() == e1j);	
			
				int which_v = EdgeHeap[e1j]->find_vert(v1->ID);
				EdgeHeap[e1j]->Verts[which_v] = v0->ID;
			}
		}
	}
}

void cwg::TetrahedralMesh::update_edges_of_vert( TetVertex* v0 )
{
	std::vector< std::vector<Sparse_Entry> >& entryset = AdjMatrix->get_entryset();
	std::vector<Sparse_Entry>& v0_edges = entryset[v0->ID];
	for (int j=0; j<v0_edges.size(); j++)
	{
		int vj = v0_edges[j].index;
		double tmpval = v0_edges[j].value;
		int e0j = tmpval;

		if (tmpval >= -0.5)
		{
			TetEdge* e = EdgeHeap[e0j];
			assert(e->id() == e0j);
			e->ComputeContraction(this);
			EdgeHeap.update(e->id());
		}
	}
}

void cwg::TetrahedralMesh::build_triangles()
{
	map< cwg::ordered_triple, vector<int> > triple_map;
	map< cwg::ordered_triple, vector<int> >::iterator triple_map_it;
	for (int k=0; k<ResultT.size(); k++)
	{
		Tetrahedron* t = &ResultT[k];
		if (t->ID == -1)
			throw exception("void cwg::TetrahedronMesh::build_triangles(): t->ID == -1");

		for (int i=0; i<4; i++)
		{
			int v0 = t->Verts[ (i+1)%4 ];
			int v1 = t->Verts[ (i+2)%4 ];
			int v2 = t->Verts[ (i+3)%4 ];
			cwg::ordered_triple triple(v0, v1, v2);

			triple_map[ triple ].push_back(t->ID);
		}
	}

	Tris.resize(triple_map.size(), NULL);
	int ntris = 0;
	for (triple_map_it = triple_map.begin(); triple_map_it != triple_map.end(); triple_map_it++)
	{
		TetTriangle* tri = new TetTriangle(triple_map_it->first.First(), triple_map_it->first.Second(), triple_map_it->first.Third());
		int tetj = 0;
		if (0)//triple_map_it->second.size() > 2)
		{
			cout << ntris << endl;
			cout << triple_map_it->first.First() << "\t" << triple_map_it->first.Second() << "\t" << triple_map_it->first.Third() << "\t" << endl;
			for (int j=0; j<triple_map_it->second.size(); j++)
			{
				cout << triple_map_it->second[j] << endl;
			}
			//throw exception("void cwg::TetrahedronMesh::build_triangles(): triple_map_it->second.size() > 2");
		}
		for (tetj = 0; tetj < triple_map_it->second.size(); tetj++)
			tri->Tets[tetj] = (triple_map_it->second)[tetj];

		Tris[ntris] = tri;
		ntris++;
	}

}

void cwg::TetrahedralMesh::build_internal_boundaries()
{
	map< cwg::ordered_triple, vector<int> > triple_map;
	map< cwg::ordered_triple, vector<int> >::iterator triple_map_it;
	tmp_tri_mesh.clear();
	TetTriangle::reset_nextid();

	for (int ti = 0; ti < ResultT.size(); ti++)
	{
		Tetrahedron* t = &ResultT[ ti ];
		for (int i=0; i<4; i++)
		{
			int v0 = t->Verts[ (i+1)%4 ];
			int v1 = t->Verts[ (i+2)%4 ];
			int v2 = t->Verts[ (i+3)%4 ];

			if ( (ResultV[v0].get_full_border_code() == Cleaver::VERTEX_INV && ResultV[v0].NLabels() == 1 )
				|| (ResultV[v1].get_full_border_code() == Cleaver::VERTEX_INV && ResultV[v1].NLabels() == 1 )
				|| (ResultV[v2].get_full_border_code() == Cleaver::VERTEX_INV && ResultV[v2].NLabels() == 1)){
					if (!(result_tris_boundary_check(v0, v1, v2)))
						if ((ResultV[v0].keep==0) && (ResultV[v1].keep == 0) && (ResultV[v2].keep == 0) )
							continue;
			}

			cwg::ordered_triple triple(v0, v1, v2);
			triple_map[ triple ].push_back(t->ID);
		}
	}

	int tri_count = 0;
	for (triple_map_it = triple_map.begin(); triple_map_it != triple_map.end(); triple_map_it++)
	{
		if (triple_map_it->second.size() == 2)
		{
			int tlb1 = ResultT[ (triple_map_it->second)[0] ].get_label();
			int tlb2 = ResultT[ (triple_map_it->second)[1] ].get_label();
			if (tlb1 != tlb2)
			{

				TetTriangle* tri = new TetTriangle(triple_map_it->first.First(), 
					triple_map_it->first.Second(), triple_map_it->first.Third());

				tri->ID = tri_count;
				++ tri_count;

				//for (int k=0; k<2; k++)
				//	tri->Tets[k] = (triple_map_it->second)[k];

				for (int tetj = 0; tetj < triple_map_it->second.size(); tetj++)
					tri->Tets[tetj] = (triple_map_it->second)[tetj];

				tmp_tri_mesh.push_back(tri);
			}else{
				if (result_tris_boundary_check(triple_map_it->first.First(), triple_map_it->first.Second(), triple_map_it->first.Third())){
						
					TetTriangle* tri = new TetTriangle(triple_map_it->first.First(), 
						triple_map_it->first.Second(), triple_map_it->first.Third());

					tri->ID = tri_count;
					++ tri_count;

					for (int tetj = 0; tetj < triple_map_it->second.size(); tetj++)
						tri->Tets[tetj] = (triple_map_it->second)[tetj];

					tmp_tri_mesh.push_back(tri);
				}else if((ResultV[triple_map_it->first.First()].keep > 0) || (ResultV[triple_map_it->first.Second()].keep > 0) || (ResultV[triple_map_it->first.Third()].keep > 0)){
					TetTriangle* tri = new TetTriangle(triple_map_it->first.First(), 
						triple_map_it->first.Second(), triple_map_it->first.Third());

					tri->ID = tri_count;
					++ tri_count;

					for (int tetj = 0; tetj < triple_map_it->second.size(); tetj++)
						tri->Tets[tetj] = (triple_map_it->second)[tetj];

					tmp_tri_mesh.push_back(tri);
				}

			}

		}
		else if (triple_map_it->second.size() == 1)
		{

			TetTriangle* tri = new TetTriangle(triple_map_it->first.First(), 
				triple_map_it->first.Second(), triple_map_it->first.Third());

			tri->Tets[0] = (triple_map_it->second)[0];
			tri->Tets[1] = -1;

			tri->ID = tri_count;
			++ tri_count;

			tmp_tri_mesh.push_back(tri);
		}
	}

	InternalMesh->Tris = tmp_tri_mesh;
	for (int i=0; i<InternalMesh->Tris.size(); i++)
	{
		//if(InternalMesh->Tris[i]->Tets[1] != -1 && InternalMesh->Tris[i]->Tets[0] != -1)
		if(InternalMesh->Tris[i]->Tets[0]== -1 || InternalMesh->Tris[i]->Tets[1]== -1 )	//on outer boundary
		{
			for (int j=0; j<3; j++)
			{
				int vj = InternalMesh->Tris[i]->Vert(j);			
				ResultV[vj].OuterBoundaryTris.push_back(InternalMesh->Tris[i]->ID);
			}
			//InternalMesh->Tris[i]->interSurf = false;
		}
		else 
		{
			for (int j=0; j<3; j++)
			{
				int vj = InternalMesh->Tris[i]->Vert(j);			
				ResultV[vj].correlate_boundary_tri(InternalMesh->Tris[i]->ID);
			}
			//InternalMesh->Tris[i]->interSurf = true;
		}

	}

}

void cwg::TetrahedralMesh::build_internal_boundaries_sliding_wenhua(){

	triple_map.clear();
	map< cwg::ordered_triple, vector<int> >::iterator triple_map_it;
	tmp_tri_mesh.clear();
	TetTriangle::reset_nextid();

	int tet_size = m_storage_tets.size();

	EdgeHeap.clear();
	delete AdjMatrix;
	AdjMatrix =  new Sparse_Matrix(m_storage_verts.size(), m_storage_verts.size(), SYM_BOTH, false, CRS);
	AdjMatrix->begin_fill_entry();

	for(int k=0; k<slide_Z; k++){
		for(int i=0; i<slide_X; i++){
			for(int j=0; j<slide_Y; j++){
				int tet_index0 = (k*slide_X*slide_Y + i*slide_Y + j)*tetn;

				for(int index=0; index<tetn; index++){
					int ti = tet_index0 +index;

					if(m_storage_tets[ti].ID != -1){
						for (int i0=0; i0<4; i0++){
							int v0 = m_storage_tets[ ti ].Verts[ (i0+1)%4 ];
							int v1 = m_storage_tets[ ti ].Verts[ (i0+2)%4 ];
							int v2 = m_storage_tets[ ti ].Verts[ (i0+3)%4 ];

							if (m_storage_verts[v0].NLabels() == 1 || m_storage_verts[v1].NLabels() == 1 || m_storage_verts[v2].NLabels() == 1){
								continue;
							}

							cwg::ordered_triple triple(v0, v1, v2);		
							triple_map[ triple ].push_back(m_storage_tets[ ti ].ID);
						}
					}	
					
				}


			}
		}
	}
	
	int tc=tmp_tri_mesh.size();
	int tri_start = tc;

	for (triple_map_it = triple_map.begin(); triple_map_it != triple_map.end(); triple_map_it++)
	{
		//Each triangle, has at most 2 labels, two of tets' IDs
		//Each tet, has at most 1 label, one ID
		if (triple_map_it->second.size() == 2)
		{
			int tlb1 = m_storage_tets[ (triple_map_it->second)[0] ].get_label();
			int tlb2 = m_storage_tets[ (triple_map_it->second)[1] ].get_label();
			if (tlb1 != tlb2)
			{
				TetTriangle* tri = new TetTriangle(triple_map_it->first.First(), 
					triple_map_it->first.Second(), triple_map_it->first.Third());

				tri->ID = tc;
				++tc;

				for (int tetj = 0; tetj < triple_map_it->second.size(); tetj++)
					tri->Tets[tetj] = (triple_map_it->second)[tetj];

				tmp_tri_mesh.push_back(tri);
			}
		}
	}

	InternalMesh->Tris = tmp_tri_mesh;

	int tmp_tri_mesh_size = tmp_tri_mesh.size();

	for (int i=tri_start; i<tmp_tri_mesh_size; i++)
	{
		for (int j=0; j<3; j++)
		{
			int vj = tmp_tri_mesh[i]->Vert(j);			
			m_storage_verts[vj].correlate_boundary_tri(tmp_tri_mesh[i]->ID);
		}
	}

}

int cwg::TetrahedralMesh::compute_valid_nverts()
{
	int nvalidverts = 0;
	for (int i=0; i<m_storage_verts.size(); i++)
	{
		if (m_storage_verts[i].ID != -1)
			nvalidverts++;
	}
	return nvalidverts;
}

int cwg::TetrahedralMesh::compute_valid_ntets()
{
	int nvalidtets = 0;
	for (int i=0; i<m_storage_tets.size(); i++)
	{
		if (m_storage_tets[i].ID != -1)
			nvalidtets++;
	}
	return nvalidtets;
}


int cwg::TetrahedralMesh::tetgenCDT(int odt_count, bool window_optimize, int starts_x, int starts_y, int starts_z)
{
	cout<<"*************** heap operation time :  "<<heap_time <<endl;
	//reset the vert labels
	int vert_size = ResultV.size();
	for(int i=0; i<vert_size; i++){
		ResultV[i].clear_label();
	}

	// related tet to v in result and reset vert label
	for ( int i = 0; i < ResultT.size(); ++i)
	{
		ResultT[i].ID = i;
		for ( int j = 0; j < 4; ++j)
		{
			ResultV[ResultT[i].Verts[j]].correlate_tet(i);
			ResultV[ResultT[i].Verts[j]].assign_label(ResultT[i].get_label());
		}
	}

	InternalMesh->Tris.clear();
	build_internal_boundaries();

	//set weight
	compute_weight_from_edge_variation();
	int odt_success = volumeODT_window(true, starts_x, starts_y, starts_z);	

	if (odt_success == -1){
		std::cout << "Did not apply CODT in this window!" << std::endl;
	}else{
		std::cout << "Applied CODT in this window!" << std::endl;
	}
	
	return odt_success;

}


// ======================== flip and split to improve the surface quality ===================================
bool cwg::TetrahedralMesh::improveSurfaceQ()
{

	int opt_cnt = 0;
	std::set< cwg::ordered_pair > EdgeSet;
	std::set< cwg::ordered_pair > EdgeSet_flip;
	std::vector<float> cost_flip;
	// find all edges
	for (int k = 0; k < InternalMesh->Tris.size(); ++k)
	{
		TetTriangle* tri = InternalMesh->Tris[k];

		for (int i=0; i<3; i++)
		{
			int vi = tri->Vert(i);

			for (int j=i+1; j<3; j++)
			{
				int vj = tri->Vert(j);
				EdgeSet.insert( cwg::ordered_pair(vi,vj));

			}

		}
	}

	// find all NLD edges -- need to be flipped (planar) or split
	const float PI = 3.1415926;
	// edges for splitting
	std::queue<cwg::ordered_pair> split_edges;

	std::cout<< "EdgeSet size: "<<EdgeSet.size()<<endl;
	for (set<cwg::ordered_pair>::iterator it = EdgeSet.begin(); it != EdgeSet.end(); it++)
	{
		// if the edge is feature line, skip
		if (ResultV[it->First()].NLabels() > 2 || ResultV[it->Second()].NLabels() > 2  // feature curve formed by 3 materials
			|| ResultV[it->First()].NLabels() == 2 && ResultV[it->First()].border_type() == 1  
			|| ResultV[it->Second()].NLabels() == 2 && ResultV[it->Second()].border_type() == 1	// outer face edge with 2 labels
			|| ResultV[it->First()].border_type() > 1 || ResultV[it->Second()].border_type() > 1)  // on boundary edge or corner
		{
			continue;
		}

		// find shared triangles of the edge
		std::vector<int> shared_tri(10);

		// outer boundary  -- flip
		if (ResultV[it->First()].border_type() == 1 && ResultV[it->Second()].border_type() == 1)
		{

			std::sort(ResultV[it->First()].OuterBoundaryTris.begin(), ResultV[it->First()].OuterBoundaryTris.end());
			std::sort(ResultV[it->Second()].OuterBoundaryTris.begin(), ResultV[it->Second()].OuterBoundaryTris.end());
			auto tit = std::set_intersection (ResultV[it->First()].OuterBoundaryTris.begin(), ResultV[it->First()].OuterBoundaryTris.end(),
				ResultV[it->Second()].OuterBoundaryTris.begin(), ResultV[it->Second()].OuterBoundaryTris.end(), shared_tri.begin());
			shared_tri.resize(tit-shared_tri.begin());

			if (shared_tri.size() == 2)
			{
				int new_e[2];

				float angle[2];
				for (int j = 0; j < 2; ++j)
				{
					for (int i = 0; i < 3; ++i)
					{
						int vi = InternalMesh->Tris[shared_tri[j]]->Vert(i);
						if (vi != it->First() && vi != it->Second() )
						{
							new_e[j] = vi;
							vec3 n1 = ResultV[it->First()].pos() - ResultV[vi].pos();
							vec3 n2 = ResultV[it->Second()].pos() - ResultV[vi].pos();

							n1=normalize(n1);
							n2=normalize(n2);

							angle[j] = acos(n1.dot(n2));
							break;
						}
					}

				}

				// check if the new edge is existing
				auto ite = EdgeSet.find(cwg::ordered_pair(new_e[0],new_e[1]));
				if (ite == EdgeSet.end())
				{
					// flip
					if (angle[0] + angle[1] > PI)
					{
						//cout<<shared_tri[0]<<" "<<shared_tri[1]<<endl;
						InternalMesh->Tris[shared_tri[0]]->SetVertices(it->First(), new_e[1], new_e[0]);
						InternalMesh->Tris[shared_tri[1]]->SetVertices(it->Second(), new_e[0], new_e[1]);

						// update tris of v -- delete tri
						for (int m = 0 ; m < ResultV[it->First()].OuterBoundaryTris.size(); ++m)
						{
							if (ResultV[it->First()].OuterBoundaryTris[m] == shared_tri[1])
							{
								ResultV[it->First()].OuterBoundaryTris.erase(ResultV[it->First()].OuterBoundaryTris.begin() + m);
							}
						}
						for (int m = 0 ; m < ResultV[it->Second()].OuterBoundaryTris.size(); ++m)
						{
							if (ResultV[it->Second()].OuterBoundaryTris[m] == shared_tri[0])
							{
								ResultV[it->Second()].OuterBoundaryTris.erase(ResultV[it->Second()].OuterBoundaryTris.begin() + m);
							}
						}

						// add tris of new edge's v
						bool find = false;
						for (int m = 0 ; m < ResultV[new_e[0]].OuterBoundaryTris.size(); ++m)
						{

							if (ResultV[new_e[0]].OuterBoundaryTris[m] == shared_tri[1])
							{
								find = true;
							}
						}
						if (!find)
						{
							ResultV[new_e[0]].OuterBoundaryTris.push_back(shared_tri[1]);
						}

						find = false;
						for (int m = 0 ; m < ResultV[new_e[0]].OuterBoundaryTris.size(); ++m)
						{

							if (ResultV[new_e[0]].OuterBoundaryTris[m] == shared_tri[0])
							{
								find = true;
							}
						}
						if (!find)
						{
							ResultV[new_e[0]].OuterBoundaryTris.push_back(shared_tri[0]);
						}

						find = false;
						for (int m = 0 ; m < ResultV[new_e[1]].OuterBoundaryTris.size(); ++m)
						{

							if (ResultV[new_e[1]].OuterBoundaryTris[m] == shared_tri[0])
							{
								find = true;
							}
						}
						if (!find)
						{
							ResultV[new_e[1]].OuterBoundaryTris.push_back(shared_tri[0]);
						}
						find = false;
						for (int m = 0 ; m < ResultV[new_e[1]].OuterBoundaryTris.size(); ++m)
						{

							if (ResultV[new_e[1]].OuterBoundaryTris[m] == shared_tri[1])
							{
								find = true;
							}
						}
						if (!find)
						{
							ResultV[new_e[1]].OuterBoundaryTris.push_back(shared_tri[1]);
						}

						++opt_cnt;
					}

				}

			}
		}
		else // inter-surface edge -- add to split
		{

			std::sort(ResultV[it->First()].OrigBoundaryTris.begin(), ResultV[it->First()].OrigBoundaryTris.end());
			std::sort(ResultV[it->Second()].OrigBoundaryTris.begin(), ResultV[it->Second()].OrigBoundaryTris.end());
			auto tit = std::set_intersection (ResultV[it->First()].OrigBoundaryTris.begin(), ResultV[it->First()].OrigBoundaryTris.end(),
				ResultV[it->Second()].OrigBoundaryTris.begin(), ResultV[it->Second()].OrigBoundaryTris.end(), shared_tri.begin());
			shared_tri.resize(tit-shared_tri.begin());

			if (shared_tri.size() == 2)
			{
				int new_e[2];
				float angle[2];
				vec3 norm[2];
				for (int j = 0; j < 2; ++j)
				{
					for (int i = 0; i < 3; ++i)
					{
						int vi = InternalMesh->Tris[shared_tri[j]]->Vert(i);
						if (vi != it->First() && vi != it->Second() )
						{
							new_e[j] = vi;
							vec3 n1 = ResultV[it->First()].pos() - ResultV[vi].pos();
							vec3 n2 = ResultV[it->Second()].pos() - ResultV[vi].pos();

							n1 = normalize(n1);
							n2 = normalize(n2);

							angle[j] = acos(n1.dot(n2));
							norm[j] = n1.cross(n2);
							norm[j] = normalize(norm[j]);
							break;
						}
					}

				}
				//cout<<"angles: "<<angle[0]<<" "<<angle[1]<<endl;

				if (angle[0] + angle[1] > PI)
				{
					// two triangles are coplanar

					if (norm[0].dot(norm[1]) == 1 || norm[0].dot(norm[1]) == -1)
						//if(CGAL::coplanar(p1,p2,p3,p4))
					{
						auto ite = EdgeSet.find(cwg::ordered_pair(new_e[0],new_e[1]));
						if (ite == EdgeSet.end())
						{

							InternalMesh->Tris[shared_tri[0]]->SetVertices(it->First(), new_e[0], new_e[1]);
							InternalMesh->Tris[shared_tri[1]]->SetVertices(it->Second(), new_e[1], new_e[0]);

							//==================================================

							// update tris of v -- delete tri
							for (int m = 0 ; m < ResultV[it->First()].OrigBoundaryTris.size(); ++m)
							{
								if (ResultV[it->First()].OrigBoundaryTris[m] == shared_tri[1])
								{
									ResultV[it->First()].OrigBoundaryTris.erase(ResultV[it->First()].OrigBoundaryTris.begin() + m);
								}
							}
							for (int m = 0 ; m < ResultV[it->Second()].OrigBoundaryTris.size(); ++m)
							{
								if (ResultV[it->Second()].OrigBoundaryTris[m] == shared_tri[0])
								{
									ResultV[it->Second()].OrigBoundaryTris.erase(ResultV[it->Second()].OrigBoundaryTris.begin() + m);
								}
							}
							// add tris of new edge's v
							bool find = false;
							{
								ResultV[new_e[0]].OrigBoundaryTris.push_back(shared_tri[1]);
							}

							{
								ResultV[new_e[1]].OrigBoundaryTris.push_back(shared_tri[0]);
							}

						}

						//cout<<"==========flip inter-face========"<<endl;
					}
					else  // add for split
					{
						split_edges.push(cwg::ordered_pair(it->First(), it->Second()));
						++opt_cnt;
					}
				}

			}

		}



	}




	// split edge
	cout<<"========== split ==========="<< split_edges.size()<<endl;
	int count = 0; 
	//while(false)
	while(split_edges.size() > 0 )
	{
		cwg::ordered_pair e = split_edges.front();
		split_edges.pop();

		// midpoint
		vec3 midv = (ResultV[e.First()].pos() + ResultV[e.Second()].pos())/2;
		// debug
		//cout<<ResultV[e->Verts[0]].pos() <<",   "<< ResultV[e->Verts[1]].pos() << ",  "<<midv<<endl;

		// distance
		float dist = length(midv - ResultV[e.First()].pos());
		// edge length
		float fdist = length(ResultV[e.Second()].pos() - ResultV[e.First()].pos());
		// find the nearest point s to the midpoint
		vec3 newv;
		int pow = 0;
		if (dist >= 1)
		{
			float tmp = dist;
			while (tmp >= 2)
			{
				tmp /= 2;
				++pow;
			}

			float ndist = std::pow(2,pow);

			vec3 dir = ResultV[e.Second()].pos() - ResultV[e.First()].pos();
			dir=normalize(dir);
			//cout<<fdist<<" "<<dist<<" "<<ndist<<endl;
			if (std::abs(dist - ndist) < std::abs(dist - ndist*2) && ndist < fdist && ndist > 0)
			{
				newv = ResultV[e.First()].pos() + ndist * dir;
			}
			else if (ndist*2 < fdist && ndist*2 > 0)
			{
				newv = ResultV[e.First()].pos() + ndist*2 * dir;
			}
			else
			{
				newv = midv;
			}
		}
		else 
		{
			float tmp = 1/dist;
			while (tmp >= 2)
			{
				tmp /= 2;
				--pow;
			}

			float ndist = std::pow(2,pow);
			vec3 dir = ResultV[e.Second()].pos() - ResultV[e.First()].pos();
			dir = normalize(dir);
			//cout<<fdist<<" "<<dist<<" "<<ndist<<endl;
			if (std::abs(dist - ndist) < std::abs(dist - ndist/2) && ndist < fdist && ndist >0)
			{
				newv = ResultV[e.First()].pos() + ndist * dir;
			}
			else if (ndist/2.0 < fdist && ndist > 0)
			{
				newv = ResultV[e.First()].pos() + ndist/2.0 * dir;
			}
			else
			{
				newv = midv;
			}

		}

		// add new vertex
		TetVertex nv;
		nv.setPos(newv);
		nv.m_labels = ResultV[e.First()].m_labels;
		nv.ID = ResultV.size();
		ResultV.push_back(nv);


		// update triangles
		// find shared triangles of the edge
		std::vector<int> shared_tri;

		std::sort(ResultV[e.First()].OrigBoundaryTris.begin(), ResultV[e.First()].OrigBoundaryTris.end());
		std::sort(ResultV[e.Second()].OrigBoundaryTris.begin(), ResultV[e.Second()].OrigBoundaryTris.end());

		std::set_intersection (ResultV[e.First()].OrigBoundaryTris.begin(), ResultV[e.First()].OrigBoundaryTris.end(),
			ResultV[e.Second()].OrigBoundaryTris.begin(), ResultV[e.Second()].OrigBoundaryTris.end(), back_inserter(shared_tri));

		// this edge has been split
		if (shared_tri.size() != 2)
		{
			cout<<"shared tri < 2 !!  = "<<shared_tri.size() << " "<< count <<endl;
			continue;
		}
		for (int k = 0; k < 2; ++k)
		{
			TetTriangle *t = InternalMesh->Tris[shared_tri[k]];
			int orgv[3];
			orgv[0] = t->Vert(0);
			orgv[1] = t->Vert(1);
			orgv[2] = t->Vert(2);


			for (int i=0; i < 3; ++i)
			{
				// replace e->v0
				if (orgv[i] == e.First())
				{
					int oldv1 = orgv[(i+1)%3];
					int oldv2 = orgv[(i+2)%3];
					t->SetVertices(nv.ID, oldv1, oldv2);

				}
				// replace e->v1 -- new tri
				else if (orgv[i] == e.Second())
				{
					int oldv1 = orgv[(i+1)%3];
					int oldv2 = orgv[(i+2)%3];
					TetTriangle *newt = new TetTriangle(nv.ID, oldv1, oldv2);
					newt->ID = InternalMesh->Tris.size();
					InternalMesh->Tris.push_back(newt);
					//cout<<newt->ID<<endl;

					// update ev0->origBoundaryTris
					for (int m = 0; m < ResultV[e.First()].OrigBoundaryTris.size(); ++m)
					{
						if (ResultV[e.First()].OrigBoundaryTris[m] == shared_tri[k])
						{
							ResultV[e.First()].OrigBoundaryTris[m] = newt->ID;
						}
					}
					if (orgv[(i+1)%3] != e.First())
					{
						ResultV[orgv[(i+1)%3]].OrigBoundaryTris.push_back(newt->ID);
					}
					else if (orgv[(i+2)%3] != e.First())
					{
						ResultV[orgv[(i+2)%3]].OrigBoundaryTris.push_back(newt->ID);
					}
				}
			}
		}

		// add new tri to new v
		ResultV[nv.ID].OrigBoundaryTris.clear();
		for (int i = 0 ; i < 2; ++i)
		{
			ResultV[nv.ID].OrigBoundaryTris.push_back(shared_tri[i]);
			ResultV[nv.ID].OrigBoundaryTris.push_back(InternalMesh->Tris.size()-1-i);
		}

		++count;
	}

	cout<< "Split count: " <<count<<endl;

	if (opt_cnt==0)
	{
		return false;
	}
	else
		return true;
	/*for (int i = 0; i < InternalMesh->Tris.size(); ++i)
	{
	TetTriangle *t = InternalMesh->Tris[i];
	vec3 n1 = ResultV[t->Vert(0)].pos()-ResultV[t->Vert(1)].pos();
	vec3 n2 = ResultV[t->Vert(0)].pos()-ResultV[t->Vert(2)].pos();
	n1 = normalize(n1);
	n2 = normalize(n2);
	if ((n1).dot(n2) == -1
	|| n1.dot(n2) == 1)
	{

	cout<<i<<endl;
	cout<< "////Wrong tri==== "<< t->Vert(0)<<" "<<t->Vert(1)<<" "<<t->Vert(2)<<endl;
	}
	}*/
}



void cwg::TetrahedralMesh::ODTsmoothSurfOnly()
{
	Cleaver::vec3 newpos;

	int vs = m_storage_verts.size();

	//move curve and edge first
	for (int i=0; i<vs; i++)
	{
		newpos.reset();
		double neighbor_num = 0;
		if(m_storage_verts[i].get_full_border_code() != Cleaver::VERTEX_INV )	//boundary surface
		{	
			if((int)m_storage_verts[i].border_type() == 1 && m_storage_verts[i].NLabels() > 1)	//on boundary face but among different materials
			{
				ODT_curve(i, newpos, neighbor_num);
			}
			else if((int)m_storage_verts[i].border_type() == 2 && m_storage_verts[i].NLabels() == 1)	//boundary edge
			{
				ODT_cube_edge(i, newpos, neighbor_num);


			}

		}


		if(neighbor_num != 0)
		{
			newpos /= neighbor_num;
			m_storage_verts[i].setPos(newpos);

			/*	vector<double> vol;
			//test foldover
			for (int j=0; j < Verts[i]->Tets.size(); j++)
			{
			vol.push_back(Tets[Verts[i]->Tets[j]]->compute_volume(this));
			}

			vec3 oldpos = Verts[i].pos();
			Verts[i]->setPos(newpos);
			bool foldover = false;
			for (int j=0; j < Verts[i]->Tets.size(); j++)
			{
			double new_vol = Tets[Verts[i]->Tets[j]]->compute_volume(this);
			if ( abs(new_vol) < 0.1 || new_vol * vol[j] < 0)
			{
			foldover = true;
			break;
			}
			}

			if(foldover)
			{
			Verts[i]->setPos(oldpos);
			}
			*/

		}
	}

	//smooth surface
	for (int i = 0; i<vs; i++)
	{

		newpos.reset();
		double neighbor_num = 0;
		if(m_storage_verts[i].get_full_border_code() != Cleaver::VERTEX_INV )	//boundary surface
		{

			if((int)m_storage_verts[i].border_type() == 1 && m_storage_verts[i].NLabels() == 1)	//boundary face
			{
				ODT_tri(i, newpos, neighbor_num);

			}


		}
		else if(m_storage_verts[i].NLabels() == 2)       //intersurface
		{

			//ODT_interface(i, newpos, neighbor_num);


		}

		if(neighbor_num != 0)
		{
			newpos /= neighbor_num;
			m_storage_verts[i].setPos(newpos);

			/*		vector<double> vol;
			//test foldover
			for (int j=0; j < Verts[i]->Tets.size(); j++)
			{
			vol.push_back(Tets[Verts[i]->Tets[j]]->compute_volume(this));
			}

			vec3 oldpos = Verts[i]->pos();
			Verts[i]->setPos(newpos);
			bool foldover = false;
			for (int j=0; j < Verts[i]->Tets.size(); j++)
			{
			double new_vol = Tets[Verts[i]->Tets[j]]->compute_volume(this);
			if ( abs(new_vol) < 0.1 || new_vol * vol[j] < 0)
			{
			foldover = true;
			break;
			}
			}

			if(foldover)
			{
			Verts[i]->setPos(oldpos);
			}
			*/
		}
	}


}

void cwg::TetrahedralMesh::ODTsmooth(bool first_smooth)
{
	//cout<<"ODT smooth === "<<endl;
	Cleaver::vec3 newpos;

	int vs = ResultV.size();
	double precision_vol = 0.1;
	//move curve and edge first
	if (first_smooth)
	{
		for (int i=0; i<vs; i++)
		{
			newpos.reset();
			double neighbor_num = 0;
			if(ResultV[i].get_full_border_code() != Cleaver::VERTEX_INV )	//boundary surface
			{	
				if((int)ResultV[i].border_type() == 1 && ResultV[i].NLabels() > 1)	//on boundary face but among different materials
				{
					ODT_curve(i, newpos, neighbor_num);
				}
				else if((int)ResultV[i].border_type() == 2 && ResultV[i].NLabels() == 1)	//boundary edge
				{
					ODT_cube_edge(i, newpos, neighbor_num);
				}

			}
			else if(ResultV[i].NLabels() == 3)       //non-manifold curve
			{

				ODT_curve(i, newpos, neighbor_num);

			}

			if(neighbor_num != 0)
			{
				newpos /= neighbor_num;

				vector<double> vol;
				//test foldover
				for (int j=0; j < ResultV[i].Tets.size(); j++)
				{
					vol.push_back(ResultT[ResultV[i].Tets[j]].compute_volume_sliding_result(this));
				}

				vec3 oldpos = ResultV[i].pos();
				ResultV[i].setPos(newpos);
				bool foldover = false;
				for (int j=0; j < ResultV[i].Tets.size(); j++)
				{
					double new_vol = ResultT[ResultV[i].Tets[j]].compute_volume_sliding_result(this);
					int flag_after = isgn(new_vol, precision_vol);
					int flag_before = isgn(vol[j], precision_vol);
					if (flag_after == 0 || flag_after != flag_before)
						//if ( abs(new_vol) < 0.1 || new_vol * vol[j] < 0)
					{
						foldover = true;
						break;
					}
				}

				if(foldover)
				{
					ResultV[i].setPos(oldpos);
				}

			}
		}
	}


	//smooth surface
	for (int i = 0; i<vs; i++)
	{

		newpos.reset();
		double neighbor_num = 0;
		if(ResultV[i].get_full_border_code() != Cleaver::VERTEX_INV )	//boundary surface
		{

			if((int)ResultV[i].border_type() == 1 && ResultV[i].NLabels() == 1)	//boundary face
			{
				ODT_tri(i, newpos, neighbor_num);	

			}


		}
		else if(ResultV[i].NLabels() == 2 && ResultV[i].get_full_border_code() == Cleaver::VERTEX_INV)       //intersurface
		{

			//ODT_interface(i, newpos, neighbor_num);


		}

		if(neighbor_num != 0)
		{
			newpos /= neighbor_num;


			vector<double> vol;
			//test foldover
			for (int j=0; j < ResultV[i].Tets.size(); j++)
			{
				vol.push_back(ResultT[ResultV[i].Tets[j]].compute_volume_sliding_result(this));
			}

			vec3 oldpos = ResultV[i].pos();
			ResultV[i].setPos(newpos);
			bool foldover = false;
			for (int j=0; j < ResultV[i].Tets.size(); j++)
			{
				double new_vol = ResultT[ResultV[i].Tets[j]].compute_volume_sliding_result(this);
				int flag_after = isgn(new_vol, precision_vol);
				int flag_before = isgn(vol[j], precision_vol);
				if (flag_after == 0 || flag_after != flag_before)
					//if ( abs(new_vol) < 0.1 || new_vol * vol[j] < 0)
				{
					foldover = true;
					break;
				}
			}

			if(foldover)
			{
				ResultV[i].setPos(oldpos);
			}

		}
	}

	//smooth tet
	for (int i = 0; i<vs; i++)
	{

		newpos.reset();
		double neighbor_num = 0;
		if(ResultV[i].get_full_border_code() == Cleaver::VERTEX_INV && ResultV[i].NLabels() == 1)	//boundary surface
		{
			if (first_smooth)
			{
				ODT_tet_no_weight(i, newpos, neighbor_num);
			}
			else
			{
				ODT_tet(i, newpos, neighbor_num);
			}
		}

		//cout<<i<<endl;
		//cout<<neighbor_num<<endl;
		if(neighbor_num != 0)
		{
			newpos /= neighbor_num;

			vector<double> vol;
			//test foldover
			for (int j=0; j < ResultV[i].Tets.size(); j++)
			{
				vol.push_back(ResultT[ResultV[i].Tets[j]].compute_volume_sliding_result(this));
			}

			vec3 oldpos = ResultV[i].pos();
			ResultV[i].setPos(newpos);
			bool foldover = false;
			for (int j=0; j < ResultV[i].Tets.size(); j++)
			{
				double new_vol = ResultT[ResultV[i].Tets[j]].compute_volume_sliding_result(this);
				int flag_after = isgn(new_vol, precision_vol);
				int flag_before = isgn(vol[j], precision_vol);
				if (flag_after == 0 || flag_after != flag_before)
					//if ( abs(new_vol) < 0.1 || new_vol * vol[j] < 0)
				{
					foldover = true;
					break;
				}
			}

			if(foldover)
			{
				ResultV[i].setPos(oldpos);
			}

		}
	}

}

void cwg::TetrahedralMesh::ODT_tet(int i, Cleaver::vec3 &newpos, double &neighbor_num)
{
	//each tet contains v
	for (int j=0; j < ResultV[i].Tets.size(); j++)
	{
		if(ResultT[ResultV[i].Tets[j]].ID == -1)
			continue;

		CGAL::Point_3<Kernel> v[4];
		double avgWeight = 0;
		for (int k=0; k < 4; k++)
		{
			int idx = ResultT[ResultV[i].Tets[j]].Verts[k];
			CGAL::Point_3<Kernel> p(ResultV[idx].pos().x,ResultV[idx].pos().y,ResultV[idx].pos().z);
			v[k] = p;

			avgWeight += ResultV[idx].weight();
		}

		CGAL::Tetrahedron_3<Kernel> t(v[0], v[1],v[2],v[3]);
		//get volume and circumcircle's center
		//double tetV = tet_vol(Tets[Verts[i]->Tets[j]]->Verts[0], Tets[Verts[i]->Tets[j]]->Verts[1],Tets[Verts[i]->Tets[j]]->Verts[2],Tets[Verts[i]->Tets[j]]->Verts[3]);
		double tetV = CGAL::volume(v[0], v[1],v[2],v[3]);
		/*if(tetV==0)
		{
		cout<<"vol=====  "<<tetV<<endl;
		system("pause");
		}
		*/
		CGAL::Point_3<Kernel> circ_c = CGAL::circumcenter(t);

		Cleaver::vec3 circ_c2;
		circ_c2.x = circ_c.x();
		circ_c2.y = circ_c.y();
		circ_c2.z = circ_c.z();

		newpos += tetV * circ_c2 * (avgWeight/4.0);

		neighbor_num += tetV * (avgWeight/4.0);
	}		
}

void cwg::TetrahedralMesh::ODT_tet_no_weight(int i, Cleaver::vec3 &newpos, double &neighbor_num)
{
	//each tet contains v
	for (int j=0; j < ResultV[i].Tets.size(); j++)
	{
		if(ResultT[ResultV[i].Tets[j]].ID == -1)
			continue;

		CGAL::Point_3<Kernel> v[4];
		double avgWeight = 0;
		for (int k=0; k < 4; k++)
		{
			int idx = ResultT[ResultV[i].Tets[j]].Verts[k];
			CGAL::Point_3<Kernel> p(ResultV[idx].pos().x,ResultV[idx].pos().y,ResultV[idx].pos().z);
			v[k] = p;

			avgWeight += ResultV[idx].weight();
		}

		CGAL::Tetrahedron_3<Kernel> t(v[0], v[1],v[2],v[3]);
		//get volume and circumcircle's center
		//double tetV = tet_vol(Tets[Verts[i]->Tets[j]]->Verts[0], Tets[Verts[i]->Tets[j]]->Verts[1],Tets[Verts[i]->Tets[j]]->Verts[2],Tets[Verts[i]->Tets[j]]->Verts[3]);
		double tetV = CGAL::volume(v[0], v[1],v[2],v[3]);
		/*if(tetV==0)
		{
		cout<<"vol=====  "<<tetV<<endl;
		system("pause");
		}
		*/
		CGAL::Point_3<Kernel> circ_c = CGAL::circumcenter(t);

		Cleaver::vec3 circ_c2;
		circ_c2.x = circ_c.x();
		circ_c2.y = circ_c.y();
		circ_c2.z = circ_c.z();

		newpos += tetV * circ_c2;// * (avgWeight/4.0);

		neighbor_num += tetV;// * (avgWeight/4.0);
	}		
}
void cwg::TetrahedralMesh::ODT_tri(int i, Cleaver::vec3 &newpos, double &neighbor_num)
{

	for (int p=0; p<ResultV[i].OuterBoundaryTris.size(); ++p)
	{
		CGAL::Point_3<Kernel> cgal_v[3];
		int tri_v[3];
		for (int j=0; j<3; ++j)
		{
			int idx = InternalMesh->Tris[ResultV[i].OuterBoundaryTris[p]]->Vert(j);
			CGAL::Point_3<Kernel> pi(ResultV[idx].pos().x,ResultV[idx].pos().y,ResultV[idx].pos().z);
			cgal_v[j] = pi;
			tri_v[j] = idx;
		}
		double triArea = tri_area(ResultV[tri_v[0]].pos(), ResultV[tri_v[1]].pos(),ResultV[tri_v[2]].pos());

		double avgWeight = (ResultV[tri_v[0]].weight() + ResultV[tri_v[1]].weight() + ResultV[tri_v[2]].weight())/3.0;

		/*if(triArea<=0)
		{
		cout<<"tri_area1: "<<triArea<<endl;
		system("pause");
		}
		*/
		CGAL::Point_3<Kernel> circ_c = CGAL::circumcenter<Kernel> (cgal_v[0],cgal_v[1],cgal_v[2] );
		Cleaver::vec3 circ_c2;
		circ_c2.x = circ_c.x();
		circ_c2.y = circ_c.y();
		circ_c2.z = circ_c.z();

		newpos += triArea * circ_c2 * avgWeight;

		neighbor_num += triArea * avgWeight;
	}




}

void cwg::TetrahedralMesh::ODT_interface(int i, Cleaver::vec3 &newpos, double &neighbor_num)
{
	//newpos += ResultV[i].pos();
	//neighbor_num += ResultV[i].weight();
	for (int p=0; p<ResultV[i].OrigBoundaryTris.size(); ++p)
	{
		CGAL::Point_3<Kernel> cgal_v[3];
		int tri_v[3];

		for (int j=0; j<3; ++j)
		{
			int idx = InternalMesh->Tris[ResultV[i].OrigBoundaryTris[p]]->Vert(j);
			CGAL::Point_3<Kernel> pi(ResultV[idx].pos().x,ResultV[idx].pos().y,ResultV[idx].pos().z);
			cgal_v[j] = pi;
			tri_v[j] = idx;


		}
		double triArea = tri_area(ResultV[tri_v[0]].pos(), ResultV[tri_v[1]].pos(),ResultV[tri_v[2]].pos());

		double avgWeight = (ResultV[tri_v[0]].weight() + ResultV[tri_v[1]].weight() + ResultV[tri_v[2]].weight())/3.0;


		CGAL::Point_3<Kernel> circ_c = CGAL::circumcenter<Kernel> (cgal_v[0],cgal_v[1],cgal_v[2] );
		Cleaver::vec3 circ_c2;
		circ_c2.x = circ_c.x();
		circ_c2.y = circ_c.y();
		circ_c2.z = circ_c.z();



		newpos += triArea * circ_c2 * avgWeight;

		neighbor_num += triArea * avgWeight;


	}
}


void cwg::TetrahedralMesh::ODT_cube_edge(int i, Cleaver::vec3 &newpos, double &neighbor_num)
{
	//newpos += Verts[i]->pos() * Verts[i]->weight();
	//neighbor_num += Verts[i]->weight();
	newpos += ResultV[i].pos();
	neighbor_num += 1;
	//each tet contains v
	for (int j=0; j < ResultV[i].Tets.size(); j++)
	{

		//newpos += Verts[i]->pos();
		//neighbor_num += 1;
		//find neighbors
		for (int k=0; k < 4; k++)
		{

			int idx = ResultT[ResultV[i].Tets[j]].Verts[k];

			if(ResultV[idx].ID == -1)
				continue;

			if (idx != i && (ResultV[idx].get_full_border_code() == ResultV[i].get_full_border_code()))
			{

				newpos += ResultV[idx].pos();// * Verts[idx]->weight();
				neighbor_num += 1;
				//neighbor_num += Verts[idx]->weight();
			}
			else if(ResultV[idx].border_type()==3)
			{
				unsigned char intersection = ResultV[idx].get_border_code() & ResultV[i].get_border_code();
				if(intersection)
				{
					int comp = 1, count = 0;
					for (int a=0; a<6; ++a)
					{
						if(intersection & comp)
							++count;
						a<<1;
					}
					if(count == 2)
					{
						newpos += ResultV[idx].pos();//* Verts[idx]->weight();
						neighbor_num += 1;
						//neighbor_num += Verts[idx]->weight();
					}
				}
			}

		}
	}


}
void cwg::TetrahedralMesh::ODT_curve(int i, Cleaver::vec3 &newpos, double &neighbor_num)
{

	//newpos += Verts[i]->pos() * Verts[i]->weight();
	//neighbor_num +=  Verts[i]->weight();
	newpos += ResultV[i].pos();
	neighbor_num += 1;
	//each tet contains v
	for (int j=0; j < ResultV[i].Tets.size(); j++)
	{

		//find neighbors
		for (int k=0; k < 4; k++)
		{
			int idx = ResultT[ResultV[i].Tets[j]].Verts[k];

			if(ResultV[idx].ID == -1)
				continue;


			if (idx != i && ResultV[idx].m_labels == ResultV[i].m_labels 
				&& ((int)ResultV[idx].get_border_code() & (int)ResultV[i].get_border_code()) != 0)
			{

				std::vector<int> sharedtri;
				sort(ResultV[i].OrigBoundaryTris.begin(),ResultV[i].OrigBoundaryTris.end());
				sort(ResultV[idx].OrigBoundaryTris.begin(),ResultV[idx].OrigBoundaryTris.end());

				set_intersection(ResultV[i].OrigBoundaryTris.begin(),ResultV[i].OrigBoundaryTris.end(), 
					ResultV[idx].OrigBoundaryTris.begin(),ResultV[idx].OrigBoundaryTris.end(),back_inserter(sharedtri) );

				if(sharedtri.size()>0)
				{			
					newpos += ResultV[idx].pos() ;//* Verts[idx]->weight();
					neighbor_num += 1;
					//neighbor_num += Verts[idx]->weight();
				}

			}


		}
	}


}

void cwg::TetrahedralMesh::private_tets( TetVertex* v0, TetVertex* v1, std::vector<int>& v0pvtets, std::vector<int>& v1pvtets )
{
	const vector<int>& v0tets = v0->Tets;
	const vector<int>& v1tets = v1->Tets;
	set_difference(v0tets.begin(), v0tets.end(), v1tets.begin(), v1tets.end(), inserter(v0pvtets, v0pvtets.end()));
	set_difference(v1tets.begin(), v1tets.end(), v0tets.begin(), v0tets.end(), inserter(v1pvtets, v1pvtets.end()));
}

void cwg::TetrahedralMesh::resize( const Cleaver::vec3& sxsysz )
{
	for (int i=0; i<m_storage_verts.size(); i++)
	{
		m_storage_verts[i].pos().x *= sxsysz.x;
		m_storage_verts[i].pos().y *= sxsysz.y;
		m_storage_verts[i].pos().z *= sxsysz.z;
	}
}

void cwg::TetrahedralMesh::translate( const Cleaver::vec3& v)
{
	for (int i=0; i<ResultV.size(); i++)
	{
		ResultV[i].pos().x += v.x;
		ResultV[i].pos().y += v.y;
		ResultV[i].pos().z += v.z;
	}
}

void cwg::TetrahedralMesh::copy_to( cwg::TetrahedralMesh* temp_mesh ) const
{
	temp_mesh->m_storage_tets = this->m_storage_tets;
	temp_mesh->m_storage_verts = this->m_storage_verts;
	temp_mesh->m_storage_verts.resize(this->m_storage_verts.size());
	temp_mesh->m_storage_tets.resize(this->m_storage_tets.size());

	for (int i=0; i<this->m_storage_verts.size(); i++)
		temp_mesh->m_storage_verts[i] = temp_mesh->m_storage_verts[i];

	for (int i=0; i<this->m_storage_tets.size(); i++)
		temp_mesh->m_storage_tets[i] = temp_mesh->m_storage_tets[i];
}

void cwg::TetrahedralMesh::copy_submesh_to( cwg::TetrahedralMesh* temp_mesh, int lb ) const
{
	if(temp_mesh->m_storage_verts.size() > 0)
		temp_mesh->release();
	const vector<int>& tet_ids_of_cur_lb = m_submesh_hashtable.at(lb);
	map<int, int> VIndexMap;
	int new_id = 0;
	for (int i=0; i<tet_ids_of_cur_lb.size(); i++)
	{
		const Tetrahedron* t = &ResultT[ tet_ids_of_cur_lb[i] ];
		for (int j=0; j<4; j++)
		{
			int id = ResultV[ t->Verts[j] ].ID;
			if(id == -1)
				cout<<"split wrong";

			if ( VIndexMap.find(id) == VIndexMap.end() )
			{
				VIndexMap.insert( make_pair(id, new_id) );
				new_id++;
			}
		}
	}
	temp_mesh->ResultV.resize(new_id);
	temp_mesh->ResultT.resize(tet_ids_of_cur_lb.size());

	for (map<int,int>::iterator it = VIndexMap.begin(); it != VIndexMap.end(); it++)
	{
		int old_id = it->first;
		int new_id = it->second;
		temp_mesh->ResultV[new_id] = ResultV[old_id];
		temp_mesh->ResultV[new_id].invalidate();
		temp_mesh->ResultV[new_id].ID = new_id;
		temp_mesh->ResultV[new_id] = temp_mesh->ResultV[new_id];
	}
	int count = 0;
	for (int i=0; i<tet_ids_of_cur_lb.size(); i++)
	{
		temp_mesh->ResultT[count] = ResultT[ tet_ids_of_cur_lb[i] ];
		for (int j=0; j<4; j++)
		{
			int old_id = temp_mesh->ResultT[count].Verts[j];
			int new_id = VIndexMap.at(old_id);
			temp_mesh->ResultT[count].Verts[j] = new_id;
		}
		temp_mesh->ResultT[count].ID = count;
		temp_mesh->ResultT[count].correlate_verts_sliding_result(temp_mesh);
		temp_mesh->ResultT[count] = temp_mesh->ResultT[count];
		count++;
	}
}

void cwg::TetrahedralMesh::split_mesh( const std::string& foldername )
{
	for (int i=0; i<ResultT.size(); i++)
	{
		int lb = ResultT[i].get_label();
		if(lb != 0){
			m_submesh_hashtable[lb].push_back(i);
		}	
	}

	std::cout << "m_submesh_size" << m_submesh_hashtable.size() << std::endl;
	char str[300];
	for (map<int,vector<int>>::iterator it = m_submesh_hashtable.begin(); 
		it != m_submesh_hashtable.end(); it++)
	{
		std::cout << "label: " << it->first << std::endl;

		TetrahedralMesh* temp_mesh = new TetrahedralMesh;
		copy_submesh_to(temp_mesh, it->first);

		temp_mesh->build_triangles();
		mkdir(m_foldername.c_str());
		mkdir((m_foldername+foldername).c_str());
		sprintf(str, "surf.%04d.ply", it->first);
		temp_mesh->save_txt_surf_ply(m_foldername+foldername+str);
		sprintf(str, "vol.%04d.ply", it->first);
		temp_mesh->save_txt_vol_ply(m_foldername+foldername+str);
		//sprintf(str, "%03d.bin.tetm", it->first);
		//temp_mesh->save_bin_tetm(m_foldername+foldername+str);
		SAFE_DELETE(temp_mesh);
	}
}

void cwg::TetrahedralMesh::test(double val)
{
	int nedges = 0;
	int nedges2 = 0;
	char str[300];
	sprintf(str, "debug%lf.test.txt", val);
	ofstream ofile(str);
	sprintf(str, "debug-map%lf.test.txt", val);
	ofstream ofile2(str);
	map<ordered_pair, int> VertexMap, VertexMap2, D1, D2;

	for (int i=0; i<EdgeHeap.size(); i++)
	{
		TetEdge* e = EdgeHeap.get_edge(i);

		if (e->id() == -1)
			cout << "void cwg::TetrahedronMesh::test(double val) error!" << endl;
		VertexMap.insert( pair<ordered_pair,int>(ordered_pair(e->Verts[0], e->Verts[1]), e->id()) );
		nedges++;		
	}

	const std::vector< vector< Sparse_Entry> >& entryset = AdjMatrix->get_entryset();
	for (int i=0; i<entryset.size(); i++)
	{
		for (int j=0; j<entryset[i].size(); j++)
		{
			int vj = entryset[i][j].index;
			double tmpeij = entryset[i][j].value;
			if (tmpeij > -1e-6)
			{
				ordered_pair tmp(i, vj);
				VertexMap2.insert( pair<ordered_pair,int>(ordered_pair(i, vj), tmpeij) );
				nedges2++;
			}
		}
	}

	set_difference(VertexMap.begin(), VertexMap.end(), VertexMap2.begin(), VertexMap2.end(), inserter(D1, D1.begin()));
	set_difference(VertexMap2.begin(), VertexMap2.end(), VertexMap.begin(), VertexMap.end(), inserter(D2, D2.begin()));

	for (auto it=D1.begin(); it != D1.end(); it++)
	{
		ofile2 << "D1: " << (it->first).First() << "\t" << (it->first).Second() << "\t" << it->second << endl;
	}
	for (auto it=D2.begin(); it != D2.end(); it++)
	{
		ofile2 << "D2: " << (it->first).First() << "\t" << (it->first).Second() << "\t" << it->second << endl;
	}

	ofile << "Edges: " << nedges2 << endl << endl;

	ofile.close();
	ofile2.close();

	if (nedges*2 != nedges2 || D1.size() != 0 || D2.size() != 0)
	{
		cout << nedges*2 << "\t" << nedges2 << endl;
		char str[300];
		sprintf(str, "val=%lf:nedges != nedges2", val);
		throw exception(str);
		SYSTEMPAUSE;
	}
	else
		cout << "pass" << endl;
}

void cwg::TetrahedralMesh::test2( double val )
{
	map< cwg::ordered_triple, vector<int> > triple_map;
	map< cwg::ordered_triple, vector<int> >::iterator triple_map_it;
	for (int k=0; k<m_storage_tets.size(); k++)
	{
		Tetrahedron* t = &m_storage_tets[k];
		if (t->ID == -1)
			continue;

		for (int i=0; i<4; i++)
		{
			int v0 = t->Verts[ (i+1)%4 ];
			int v1 = t->Verts[ (i+2)%4 ];
			int v2 = t->Verts[ (i+3)%4 ];
			cwg::ordered_triple triple(v0, v1, v2);

			triple_map_it = triple_map.find( triple );
			if (triple_map_it != triple_map.end())
				triple_map_it->second.push_back(t->ID);
			else
				triple_map.insert(make_pair(triple, vector<int>(1,t->ID)));
		}
	}

	vector<TetTriangle*> TmpTris(triple_map.size(), NULL);
	int ntris = 0;
	for (triple_map_it = triple_map.begin(); triple_map_it != triple_map.end(); triple_map_it++)
	{
		TetTriangle* tri = new TetTriangle(triple_map_it->first.First(), triple_map_it->first.Second(), triple_map_it->first.Third());
		int tetj = 0;
		if (triple_map_it->second.size() > 2)
		{
			cout << ntris << endl;
			cout << triple_map_it->first.First() << "\t" << triple_map_it->first.Second() << "\t" << triple_map_it->first.Third() << "\t" << endl;
			for (int j=0; j<triple_map_it->second.size(); j++)
			{
				cout << triple_map_it->second[j] << endl;
			}
			throw exception("void cwg::TetrahedronMesh::build_triangles(): triple_map_it->second.size() > 2");
		}
		for (tetj = 0; tetj < triple_map_it->second.size(); tetj++)
			tri->Tets[tetj] = (triple_map_it->second)[tetj];

		TmpTris[ntris] = tri;
		ntris++;
	}
	delete_vector_elems(TmpTris);
}

void cwg::TetrahedralMesh::reset_id()
{
	TetVertex::reset_nextid();
	TetEdge::reset_nextid();
	TetTriangle::reset_nextid();
	Tetrahedron::reset_nextid();
}


struct vec3comp
{
	bool operator() (const Cleaver::vec3& left, const Cleaver::vec3& right) const
	{
		const double eps = 1e-3;
		bool flag = false;
		if (left.z < right.z - eps)
			return true;
		else
		{
			// 1. equal z:
			if (left.z < right.z + eps)
			{
				// check y
				if (left.y < right.y - eps)
					return true;
				else
				{
					if (left.y < right.y + eps)
					{
						// check x
						if (left.x < right.x - eps)
							return true;
						else
							return false;								
					}
					else
						return false;
				}
			}
			else // 2. left.z > right.z + eps
				return false;
		}
	}
};


void cwg::TetrahedralMesh::merge( const std::vector<TetrahedralMesh*>& submeshes )
{
	map< Cleaver::vec3, int, vec3comp > BoundaryPointsMap;
	map< Cleaver::vec3, int, vec3comp >::iterator it;

	vector< map<int,int> > VIndexMapVec;
	int glbid = 0;
	int num_internal_points = 0;
	int num_tets = 0;
	for (int i=0; i<submeshes.size(); i++)
	{
		TetrahedralMesh* tmp = submeshes[i];
		num_tets += tmp->m_storage_tets.size();
		for (int vi=0; vi<tmp->m_storage_verts.size(); vi++)
		{
			unsigned char fbcode = tmp->m_storage_verts[vi].get_full_border_code();
			if (fbcode == Cleaver::VERTEX_INV)
				num_internal_points++;
			else
			{
				Cleaver::vec3 tmppos = tmp->m_storage_verts[vi].pos();
				it = BoundaryPointsMap.find(tmppos);
				if (it == BoundaryPointsMap.end())
					BoundaryPointsMap.insert( make_pair(tmp->m_storage_verts[vi].pos(), glbid++) );
			}
		}
	}

	int total_points = BoundaryPointsMap.size() + num_internal_points;
	m_storage_verts.resize(total_points);
	m_storage_tets.resize(num_tets);

	int count = 0;
	for (int i=0; i<submeshes.size(); i++)
	{
		TetrahedralMesh* tmp = submeshes[i];
		map<int, int> VIndexMap;
		for (int vi=0; vi<tmp->m_storage_verts.size(); vi++)
		{
			unsigned char fbcode = tmp->m_storage_verts[vi].get_full_border_code();
			if (fbcode == Cleaver::VERTEX_INV)
			{
				int storage_pos = BoundaryPointsMap.size() + count;
				m_storage_verts[storage_pos].pos() = tmp->m_storage_verts[vi].pos();
				m_storage_verts[storage_pos].ID = storage_pos;
				VIndexMap.insert( make_pair(tmp->m_storage_verts[vi].ID, storage_pos) );
				count++;
			}
			else
			{
				it = BoundaryPointsMap.find(tmp->m_storage_verts[vi].pos());
				if (it == BoundaryPointsMap.end())
					throw exception("void cwg::TetrahedralMesh::merge( const std::vector<TetrahedralMesh*>& submeshes ):\
									BoundaryPointsMap not found!");
				else
				{
					glbid = it->second;
					// Submesh: i, Vert: vi, ==> glbid
					m_storage_verts[glbid].pos() = tmp->m_storage_verts[vi].pos();
					m_storage_verts[glbid].ID = glbid;
					VIndexMap.insert( make_pair(tmp->m_storage_verts[vi].ID, glbid) );
				}
			}
		}
		VIndexMapVec.push_back(VIndexMap);
	}

	count = 0;
	for (int i=0; i<submeshes.size(); i++)
	{
		TetrahedralMesh* tmp = submeshes[i];
		for (int ti=0; ti<tmp->m_storage_tets.size(); ti++)
		{
			m_storage_tets[count].ID = count;
			m_storage_tets[count].set_label(tmp->m_storage_tets[ti].get_label());
			for (int j=0; j<4; j++)
			{
				int tivj = tmp->m_storage_tets[ti].Verts[j];
				int tivj_mapped_id = VIndexMapVec[i].at(tivj);
				m_storage_tets[count].Verts[j] = tivj_mapped_id;
			}
			count++;
		}
	}

	this->m_storage_verts.resize( m_storage_verts.size() );
	this->m_storage_tets.resize( m_storage_tets.size() );
	for (int i=0; i<m_storage_verts.size(); i++)
		this->m_storage_verts[i] = m_storage_verts[i];
	for (int i=0; i<m_storage_tets.size(); i++)
		this->m_storage_tets[i] = m_storage_tets[i];

	for (int i=0; i<this->m_storage_tets.size(); i++)
		this->m_storage_tets[i].correlate_verts(this);

	double xmin, xmax, ymin, ymax, zmin, zmax;
	xmin = ymin = zmin = numeric_limits<double>::max();
	xmax = ymax = zmax = numeric_limits<double>::min();
	for (int i=0; i<this->m_storage_verts.size(); i++)
	{
		this->m_storage_verts[i].update_funcval();
		xmin = std::min(xmin, this->m_storage_verts[i].pos().x);
		xmax = std::max(xmax, this->m_storage_verts[i].pos().x);
		ymin = std::min(ymin, this->m_storage_verts[i].pos().y);
		ymax = std::max(ymax, this->m_storage_verts[i].pos().y);
		zmin = std::min(zmin, this->m_storage_verts[i].pos().z);
		zmax = std::max(zmax, this->m_storage_verts[i].pos().z);
		this->m_storage_verts[i].set_full_border_code(Cleaver::VERTEX_UDV);
	}

	m_X = xmax + 1;
	m_Y = ymax + 1;
	m_Z = zmax + 1;

	compute_verts_border_code();
}

int cwg::TetrahedralMesh::volumeODT_window(bool last, int starts_x, int starts_y, int starts_z)
{
	tetgenio tetin, tetout, addin;
	// build tetgenio from boundary

	tetin.numberofpoints = ResultV.size() ;
	tetin.pointlist = new REAL[tetin.numberofpoints * 3] ;
	tetin.pointmarkerlist = new int[tetin.numberofpoints] ;
	tetin.numberoffacets = InternalMesh->Tris.size();

	tetin.facetlist = new tetgenio::facet[tetin.numberoffacets] ;
	tetin.numberofregions = 0; //assume there is no hole.
	tetin.numberofholes = 0;

	// vertex
	for(int i=0; i<ResultV.size(); ++i) {
		tetin.pointlist[i*3]   = ResultV[i].pos().x ; //vit->point()[0];
		tetin.pointlist[i*3+1] = ResultV[i].pos().y ; //vit->point()[1];
		tetin.pointlist[i*3+2] = ResultV[i].pos().z ; //it->point()[2];

		if(ResultV[i].NLabels()>1 || ResultV[i].border_type() > 0 || ResultV[i].boundary || ResultV[i].slideBoundary || ResultV[i].keep > 0)
		{
			tetin.pointmarkerlist[i] = -1;
		}
		else
		{
			tetin.pointmarkerlist[i] = ResultV[i].get_label(0); 
		}

	}

	// face
	tetgenio::facet *f;
	tetgenio::polygon *p;

	for(unsigned int i = 0; i < InternalMesh->Tris.size(); i++) 
	{
		f = &tetin.facetlist[i]; 
		tetin.init(f);      
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[1];
		p = &f->polygonlist[0];
		tetin.init(p);
		p->numberofvertices = 3; 
		// Allocate memory for face vertices
		p->vertexlist = new int[p->numberofvertices];
		p->vertexlist[0] = InternalMesh->Tris[i]->Vert(0);
		p->vertexlist[1] = InternalMesh->Tris[i]->Vert(1);
		p->vertexlist[2] = InternalMesh->Tris[i]->Vert(2);
	}

	::tetrahedralize("p/0.000001YQMO2/3", &tetin, &tetout); // "pYnQ" is the parameter for Tetgen 

	int npoints = tetout.numberofpoints;
	int nedges = tetout.numberofedges ;
	int nfacets = tetout.numberoffacets ;
	int ntris = tetout.numberoftrifaces ;
	int ntets = tetout.numberoftetrahedra ;

	// update vertex
	//std::cout<<ResultV.size()<<" "<<npoints<<endl;

	if (npoints > ResultV.size()){
		return -1;
	}
	for (int i=0; i<ResultV.size(); ++i)
	{

		Cleaver::vec3 p;
		p.x = tetout.pointlist[i*3];
		p.y = tetout.pointlist[i*3+1];
		p.z = tetout.pointlist[i*3+2];

		ResultV[i].setPos(p);
		ResultV[i].Tets.clear();

	}

	int ii=ResultV.size();
	int orig_points_size = ResultV.size();

	// update tet
	ResultT.clear();
	ResultT.resize(ntets);

	// tet cannot defined label
	std::vector<int> failtet;

	for(int i=0; i<ntets; ++i) 
	{
		Tetrahedron* t = &ResultT[i];
		bool labeled = false;
		bool all_old = true;
		for(int j=0; j<4; ++j) 
		{
			int vi = tetout.tetrahedronlist[i*4+j] ;
			t->Verts[j] = vi ;
			t->ID = i;
			if(ResultV[vi].NLabels() == 1)
			{
				t->set_label(ResultV[vi].get_label(0));
				labeled = true;
			}

			if (vert_window_map_inverse[vi] < 0)
				all_old = false;
			ResultV[vi].Tets.push_back(t->ID);
		}

		if(!all_old){
			for(int j=0; j<4; ++j) {
				int vi = tetout.tetrahedronlist[i*4+j] ;
			
				if(vert_window_map_inverse[vi] < 0){
					int vert_count;
					int vert_count_x, vert_count_y, vert_count_z;

					for(int adjust_vert_index=0; adjust_vert_index<3; adjust_vert_index++){
						if( tetout.tetrahedronlist[i*4+(j+adjust_vert_index)%4] < orig_points_size){
							vert_count = vert_window_map_inverse[tetout.tetrahedronlist[i*4+(j+adjust_vert_index)%4]];
							vert_count_z = (vert_count / slide_X / slide_Y);
							vert_count_x = (vert_count % (slide_X * slide_Y)) / slide_Y;
							vert_count_y = (vert_count % (slide_X * slide_Y)) % slide_Y;

							if(vert_count_z > 0){
								vert_count_z -= 1;
								vert_count = vert_count_z * (slide_X * slide_Y) + vert_count_x * slide_Y + vert_count_y;
								if( m_storage_verts[vert_count].ID == -1){
									vert_window_map_inverse[vi] = vert_count;
									break;
								}else{
									if(vert_count_x > 0){
										vert_count_x -= 1;
										vert_count = vert_count_z * (slide_X * slide_Y) + vert_count_x * slide_Y + vert_count_y;
										if( m_storage_verts[vert_count].ID == -1){
											vert_window_map_inverse[vi] = vert_count;
											break;
										}else{
											if(vert_count_y > 0){
												vert_count_y -= 1;
												vert_count = vert_count_z * (slide_X * slide_Y) + vert_count_x * slide_Y + vert_count_y;
												if( m_storage_verts[vert_count].ID == -1){
													vert_window_map_inverse[vi] = vert_count;
													break;
												}
											}

										}
									}

								}
							}
						}			
					}
				}
			}
		}

		if(last && !labeled)
		{
			t->set_label(-5);
			failtet.push_back(i);
		}

		ResultT[i] = *t;
	}


	for(int check_index=orig_points_size; check_index < npoints; check_index++)
		if(vert_window_map_inverse[check_index] < 0)
			throw exception("added verts does not get correct id");

	if(last)
	{
		// tet cannot defined label
		int failcount = 0;

		//set label for boundary tet
		for (int i = 0; i < failtet.size(); ++i )
		{
			int surfcount = 0;
			bool labeled = false;

			sort(ResultV[ResultT[failtet[i]].Verts[0]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[0]].m_labels.end());
			sort(ResultV[ResultT[failtet[i]].Verts[1]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[1]].m_labels.end());
			sort(ResultV[ResultT[failtet[i]].Verts[2]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[2]].m_labels.end());
			sort(ResultV[ResultT[failtet[i]].Verts[3]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[3]].m_labels.end());

			std::vector<int> sharedlabel_tmp;
			std::set_intersection(ResultV[ResultT[failtet[i]].Verts[0]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[0]].m_labels.end(), 
				ResultV[ResultT[failtet[i]].Verts[1]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[1]].m_labels.end(),back_inserter(sharedlabel_tmp) );

			std::vector<int> sharedlabel_tmp2;
			std::set_intersection(ResultV[ResultT[failtet[i]].Verts[2]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[2]].m_labels.end(), 
				ResultV[ResultT[failtet[i]].Verts[3]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[3]].m_labels.end(),back_inserter(sharedlabel_tmp2) );

			std::vector<int> sharedlabel;
			std::set_intersection(sharedlabel_tmp.begin(),sharedlabel_tmp.end(), 
				sharedlabel_tmp2.begin(),sharedlabel_tmp2.end(),back_inserter(sharedlabel) );

			if(sharedlabel.size() == 1)
			{
				ResultT[failtet[i]].set_label(sharedlabel[0]);
				labeled = true;

			}
			else
			{
				//std::cout << "debug here" << std::endl;
				float center_x=0, center_y=0, center_z =0;
				for(int vert_index=0; vert_index<4; vert_index++){
					center_x += ResultV[ResultT[failtet[i]].Verts[vert_index]].pos().x;
					center_y += ResultV[ResultT[failtet[i]].Verts[vert_index]].pos().y;
					center_z += ResultV[ResultT[failtet[i]].Verts[vert_index]].pos().z;

				}
				center_x /= 4.;
				center_y /= 4.;
				center_z /= 4.;

				Cleaver::vec3 center;
				center.x = center_x; center.y = center_y; center.z = center_z;

				if(!labeled){
					int while_count = 0; bool first_try = true;
					while(!labeled && while_count<20){
						if(!first_try){
							center_x=0, center_y=0, center_z =0;
							float alphas[4] = {0., 0., 0., 0.};
							for(int alpha_index=0; alpha_index<4; alpha_index++){
								int rand_num = rand();
								while(rand_num ==0 || rand_num == RAND_MAX){
									rand_num = rand();
								}
								alphas[alpha_index] = static_cast <float> (rand_num) / static_cast <float> (RAND_MAX);
							}

							float weights[4] = {0., 0., 0., 0.};
							weights[0] = pow(sin(alphas[0]*2*3.14), 2);
							weights[1] = (1 - weights[0]) * pow(sin(alphas[1]*2*3.14), 2);
							weights[2] = (1 - weights[0]) * pow(cos(alphas[1]*2*3.14), 2) * pow(sin(alphas[2]*2*3.14), 2);
							weights[3] = (1 - weights[0]) * pow(cos(alphas[1]*2*3.14), 2) * pow(cos(alphas[2]*2*3.14), 2);

							for(int vert_index=0; vert_index<4; vert_index++){
								center_x += weights[vert_index] * ResultV[ResultT[failtet[i]].Verts[vert_index]].pos().x;
								center_y += weights[vert_index] * ResultV[ResultT[failtet[i]].Verts[vert_index]].pos().y;
								center_z += weights[vert_index] * ResultV[ResultT[failtet[i]].Verts[vert_index]].pos().z;
							}
						}

						int tet_size = ResultT_tmp.size();
						for(int tet_count=0; tet_count<tet_size; tet_count++){
							bool inside = ResultT_tmp[tet_count].inside_tet_check(center_x, center_y, center_z, &ResultV_tmp);
							if(inside){
								int tl = ResultT_tmp[tet_count].get_label();
								ResultT[failtet[i]].set_label(tl);
								labeled=true;
								break;
							}
						}
						while_count ++;
						first_try = false;
					}
				}

				if(!labeled)
				{
					int label_size = sharedlabel.size();
					ofstream file_out0("debug3.obj");		
					if (!file_out0)
						throw exception("Can not open debug file!");
					file_out0 << "v " << ResultV[ResultT[failtet[i]].Verts[0]].pos().x << " " << ResultV[ResultT[failtet[i]].Verts[0]].pos().y << " " << ResultV[ResultT[failtet[i]].Verts[0]].pos().z << std::endl;
					file_out0 << "v " << ResultV[ResultT[failtet[i]].Verts[1]].pos().x << " " << ResultV[ResultT[failtet[i]].Verts[1]].pos().y << " " << ResultV[ResultT[failtet[i]].Verts[1]].pos().z << std::endl;
					file_out0 << "v " << ResultV[ResultT[failtet[i]].Verts[2]].pos().x << " " << ResultV[ResultT[failtet[i]].Verts[2]].pos().y << " " << ResultV[ResultT[failtet[i]].Verts[2]].pos().z << std::endl;
					file_out0 << "v " << ResultV[ResultT[failtet[i]].Verts[3]].pos().x << " " << ResultV[ResultT[failtet[i]].Verts[3]].pos().y << " " << ResultV[ResultT[failtet[i]].Verts[3]].pos().z << std::endl;
					file_out0 << "v " << center_x << " " << center_y << " " << center_z << std::endl;
					file_out0 << "l " << 1 << " " << 2 << std::endl; file_out0 << "l " << 1 << " " << 3 << std::endl; file_out0 << "l " << 1 << " " << 4 << std::endl;
					file_out0 << "l " << 2 << " " << 3 << std::endl; file_out0 << "l " << 2 << " " << 4 << std::endl; file_out0 << "l " << 3 << " " << 4 << std::endl;		
					file_out0.close();


					int tet_size = ResultT_tmp.size();
					
					for(int label_index=0; label_index<label_size; label_index++){
						string file_out_name_0 = "debug_" + std::to_string(label_index) + ".obj";
						ofstream file_out_0(file_out_name_0);
						if (!file_out_0)
							throw exception("Can not open debug file!");
						int label_0_count = 0;

						for(int tet_count=0; tet_count<tet_size; tet_count++){
							if(ResultT_tmp[tet_count].get_label() == sharedlabel[label_index]){
								file_out_0 << "v " << ResultV_tmp[ResultT_tmp[tet_count].Verts[0]].pos().x << " " << ResultV_tmp[ResultT_tmp[tet_count].Verts[0]].pos().y << " " << ResultV_tmp[ResultT_tmp[tet_count].Verts[0]].pos().z << std::endl;
								file_out_0 << "v " << ResultV_tmp[ResultT_tmp[tet_count].Verts[1]].pos().x << " " << ResultV_tmp[ResultT_tmp[tet_count].Verts[1]].pos().y << " " << ResultV_tmp[ResultT_tmp[tet_count].Verts[1]].pos().z << std::endl;
								file_out_0 << "v " << ResultV_tmp[ResultT_tmp[tet_count].Verts[2]].pos().x << " " << ResultV_tmp[ResultT_tmp[tet_count].Verts[2]].pos().y << " " << ResultV_tmp[ResultT_tmp[tet_count].Verts[2]].pos().z << std::endl;
								file_out_0 << "v " << ResultV_tmp[ResultT_tmp[tet_count].Verts[3]].pos().x << " " << ResultV_tmp[ResultT_tmp[tet_count].Verts[3]].pos().y << " " << ResultV_tmp[ResultT_tmp[tet_count].Verts[3]].pos().z << std::endl;
								label_0_count ++;
							}
						}
						for(int index =0; index<label_0_count; index++){
							file_out_0 << "l " << index*4+1 << " " << index*4+2 << std::endl; file_out_0 << "l " << index*4+1 << " " << index*4+3 << std::endl; file_out_0 << "l " << index*4+1 << " " << index*4+4 << std::endl;
							file_out_0 << "l " << index*4+2 << " " << index*4+3 << std::endl; file_out_0 << "l " << index*4+2 << " " << index*4+4 << std::endl; file_out_0 << "l " << index*4+3 << " " << index*4+4 << std::endl;		
						}

						file_out_0.close();

					}
				}

			}
		}
	}

	return 0;
}


void cwg::TetrahedralMesh::volumeODT(bool last)
{
	tetgenio tetin, tetout, addin;
	// build tetgenio from boundary

	tetin.numberofpoints = ResultV.size() ;
	tetin.pointlist = new REAL[tetin.numberofpoints * 3] ;
	tetin.pointmarkerlist = new int[tetin.numberofpoints] ;
	tetin.numberoffacets = InternalMesh->Tris.size();

	tetin.facetlist = new tetgenio::facet[tetin.numberoffacets] ;
	//assumen there is no hole ! If go as I expected, there would be holes.
	tetin.numberofregions = 0; //assume there is no hole.
	tetin.numberofholes = 0;

	// vertex
	for(int i=0; i<ResultV.size(); ++i) {
		tetin.pointlist[i*3]   = ResultV[i].pos().x ; //vit->point()[0];
		tetin.pointlist[i*3+1] = ResultV[i].pos().y ; //vit->point()[1];
		tetin.pointlist[i*3+2] = ResultV[i].pos().z ; //it->point()[2];

		if(ResultV[i].NLabels()>1 || ResultV[i].border_type() > 0)
		{
			tetin.pointmarkerlist[i] = -1;
		}
		else
		{
			tetin.pointmarkerlist[i] = ResultV[i].get_label(0); 
		}

	}

	// face
	tetgenio::facet *f;
	tetgenio::polygon *p;

	int fcount = 0;
	for(unsigned int i = 0; i < InternalMesh->Tris.size(); i++) 
	{
		f = &tetin.facetlist[i]; 
		tetin.init(f);      
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[1];
		p = &f->polygonlist[0];
		tetin.init(p);
		p->numberofvertices = 3; 
		// Allocate memory for face vertices
		p->vertexlist = new int[p->numberofvertices];
		p->vertexlist[0] = InternalMesh->Tris[i]->Vert(0);
		p->vertexlist[1] = InternalMesh->Tris[i]->Vert(1);
		p->vertexlist[2] = InternalMesh->Tris[i]->Vert(2);
	}

	::tetrahedralize("pYQ", &tetin, &tetout); // "pYnQ" is the parameter for Tetgen 

	int npoints = tetout.numberofpoints;
	int nedges = tetout.numberofedges ;
	int nfacets = tetout.numberoffacets ;
	int ntris = tetout.numberoftrifaces ;
	int ntets = tetout.numberoftetrahedra ;

	// update vertex
	cout<<ResultV.size()<<" "<<npoints<<endl;

	if (npoints > ResultV.size()){
	}
	for (int i=0; i<ResultV.size(); ++i)
	{

		Cleaver::vec3 p;
		p.x = tetout.pointlist[i*3];
		p.y = tetout.pointlist[i*3+1];
		p.z = tetout.pointlist[i*3+2];

		ResultV[i].setPos(p);
		ResultV[i].Tets.clear();

	}

	int ii=ResultV.size();

	if(npoints > ResultV.size())
	{
		cout<<"tetgen add v"<<endl;

		ResultV.resize(npoints);
		for (; ii<npoints; ++ii)
		{

			Cleaver::vec3 p;
			p.x = tetout.pointlist[ii*3];
			p.y = tetout.pointlist[ii*3+1];
			p.z = tetout.pointlist[ii*3+2];

			ResultV[ii].ID = ii;
			ResultV[ii].m_labels.clear();
			ResultV[ii].assign_label(tetout.pointmarkerlist[ii]);
			ResultV[ii].setPos(p);
			ResultV[ii].Tets.clear();


		}
	}

	// update tet

	ResultT.clear();
	ResultT.resize(ntets);

	// tet cannot defined label
	std::vector<int> failtet;

	for(int i=0; i<ntets; ++i) 
	{
		Tetrahedron* t = &ResultT[i];
		bool labeled = false;
		for(int j=0; j<4; ++j) 
		{
			int vi = tetout.tetrahedronlist[i*4+j] ;
			t->Verts[j] = vi ;
			t->ID = i;
			if(ResultV[vi].NLabels() == 1)
			{
				t->set_label(ResultV[vi].get_label(0));
				labeled = true;
			}


			ResultV[vi].Tets.push_back(t->ID);
		}

		if(last && !labeled)
		{
			t->set_label(-5);
			failtet.push_back(i);
		}

		ResultT[i] = *t;


	}

	if(last)
	{
		int failcount = 0;
		int count = failtet.size();
		//set label for boundary tet
		for (int i = 0; i < failtet.size(); ++i )
		{
			int surfcount = 0;

			bool labeled = false;


			sort(ResultV[ResultT[failtet[i]].Verts[0]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[0]].m_labels.end());
			sort(ResultV[ResultT[failtet[i]].Verts[1]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[1]].m_labels.end());
			sort(ResultV[ResultT[failtet[i]].Verts[2]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[2]].m_labels.end());
			sort(ResultV[ResultT[failtet[i]].Verts[3]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[3]].m_labels.end());

			std::vector<int> sharedlabel_tmp;
			std::set_intersection(ResultV[ResultT[failtet[i]].Verts[0]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[0]].m_labels.end(), 
				ResultV[ResultT[failtet[i]].Verts[1]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[1]].m_labels.end(),back_inserter(sharedlabel_tmp) );

			std::vector<int> sharedlabel_tmp2;
			std::set_intersection(ResultV[ResultT[failtet[i]].Verts[2]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[2]].m_labels.end(), 
				ResultV[ResultT[failtet[i]].Verts[3]].m_labels.begin(),ResultV[ResultT[failtet[i]].Verts[3]].m_labels.end(),back_inserter(sharedlabel_tmp2) );

			std::vector<int> sharedlabel;
			std::set_intersection(sharedlabel_tmp.begin(),sharedlabel_tmp.end(), 
				sharedlabel_tmp2.begin(),sharedlabel_tmp2.end(),back_inserter(sharedlabel) );


			if(sharedlabel.size() == 1)
			{
				ResultT[failtet[i]].set_label(sharedlabel[0]);
				labeled = true;

			}
			else
			{
				std::set<int> surfRelateLabel;	//surface triangle adjacent tet labels
				std::vector<int> candidateLabel;
				std::vector<std::vector<int>> candidateLabelset;	//candidate label from surface triangle


				std::vector<int> sharedtet;
				std::vector<int> sharedtet_tmp;

				std::vector<int> sharedtri;
				std::vector<int> sharedtri_tmp;
				sort(ResultV[ResultT[failtet[i]].Verts[0]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[0]].OrigBoundaryTris.end());
				sort(ResultV[ResultT[failtet[i]].Verts[1]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[1]].OrigBoundaryTris.end());
				sort(ResultV[ResultT[failtet[i]].Verts[2]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[2]].OrigBoundaryTris.end());
				sort(ResultV[ResultT[failtet[i]].Verts[3]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[3]].OrigBoundaryTris.end());

				sort(ResultV[ResultT[failtet[i]].Verts[0]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[0]].Tets.end());
				sort(ResultV[ResultT[failtet[i]].Verts[1]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[1]].Tets.end());
				sort(ResultV[ResultT[failtet[i]].Verts[2]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[2]].Tets.end());
				sort(ResultV[ResultT[failtet[i]].Verts[3]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[3]].Tets.end());



				std::set_intersection(ResultV[ResultT[failtet[i]].Verts[0]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[0]].OrigBoundaryTris.end(), 
					ResultV[ResultT[failtet[i]].Verts[1]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[1]].OrigBoundaryTris.end(),back_inserter(sharedtri_tmp) );
				std::set_intersection(sharedtri_tmp.begin(),sharedtri_tmp.end(), 
					ResultV[ResultT[failtet[i]].Verts[2]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[2]].OrigBoundaryTris.end(),back_inserter(sharedtri) );

				//===================
				std::set_intersection(ResultV[ResultT[failtet[i]].Verts[0]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[0]].Tets.end(), 
					ResultV[ResultT[failtet[i]].Verts[1]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[1]].Tets.end(),back_inserter(sharedtet_tmp) );
				std::set_intersection(sharedtet_tmp.begin(),sharedtet_tmp.end(), 
					ResultV[ResultT[failtet[i]].Verts[2]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[2]].Tets.end(),back_inserter(sharedtet) );

				//std::cout<<"** shared tet num 1: "<<sharedtet.size()<<endl;

				for(int p = 0; p < sharedtet.size(); ++p)
				{
					//std::cout<<"shared tet: "<<Tets[sharedtet[p]]->get_label()<<endl;
					int tl = ResultT[sharedtet[p]].get_label();
					if (tl >=0)
					{
						if(sharedtri.size() == 0)
						{
							ResultT[failtet[i]].set_label(tl);
							labeled = true;
							break;
						}
						else
						{
							surfRelateLabel.insert(tl);
							surfcount++;

						}

					}
				}

				sharedtri_tmp.clear();
				sharedtri.clear();

				std::set_intersection(ResultV[ResultT[failtet[i]].Verts[0]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[0]].OrigBoundaryTris.end(), 
					ResultV[ResultT[failtet[i]].Verts[2]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[2]].OrigBoundaryTris.end(),back_inserter(sharedtri_tmp) );
				std::set_intersection(sharedtri_tmp.begin(),sharedtri_tmp.end(), 
					ResultV[ResultT[failtet[i]].Verts[3]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[3]].OrigBoundaryTris.end(),back_inserter(sharedtri) );


				//================
				sharedtet_tmp.clear();
				sharedtet.clear();

				std::set_intersection(ResultV[ResultT[failtet[i]].Verts[0]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[0]].Tets.end(), 
					ResultV[ResultT[failtet[i]].Verts[3]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[3]].Tets.end(),back_inserter(sharedtet_tmp) );
				std::set_intersection(sharedtet_tmp.begin(),sharedtet_tmp.end(), 
					ResultV[ResultT[failtet[i]].Verts[2]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[2]].Tets.end(),back_inserter(sharedtet) );

				if(!labeled)
				{
					for(int p = 0; p < sharedtet.size(); ++p)
					{

						int tl = ResultT[sharedtet[p]].get_label();
						if (tl>=0)
						{
							if(sharedtri.size() == 0)
							{
								ResultT[failtet[i]].set_label(tl);
								labeled = true;
								break;
							}
							else
							{
								surfRelateLabel.insert(tl);
								surfcount++;
							}

						}
					}
				}

				sharedtri_tmp.clear();
				sharedtri.clear();
				std::set_intersection(ResultV[ResultT[failtet[i]].Verts[1]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[1]].OrigBoundaryTris.end(), 
					ResultV[ResultT[failtet[i]].Verts[2]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[2]].OrigBoundaryTris.end(),back_inserter(sharedtri_tmp) );
				std::set_intersection(sharedtri_tmp.begin(),sharedtri_tmp.end(), 
					ResultV[ResultT[failtet[i]].Verts[3]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[3]].OrigBoundaryTris.end(),back_inserter(sharedtri) );

				//================
				sharedtet_tmp.clear();
				sharedtet.clear();

				std::set_intersection(ResultV[ResultT[failtet[i]].Verts[1]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[1]].Tets.end(), 
					ResultV[ResultT[failtet[i]].Verts[3]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[3]].Tets.end(),back_inserter(sharedtet_tmp) );
				std::set_intersection(sharedtet_tmp.begin(),sharedtet_tmp.end(), 						
					ResultV[ResultT[failtet[i]].Verts[2]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[2]].Tets.end(),back_inserter(sharedtet) );

				//std::cout<<"** shared tet num 3: "<<sharedtet.size()<<endl;
				if(!labeled)
				{
					for(int p = 0; p < sharedtet.size(); ++p)
					{
						//std::cout<<"shared tet: "<<Tets[sharedtet[p]]->get_label()<<endl;

						int tl = ResultT[sharedtet[p]].get_label();
						if (tl>=0)
						{
							if(sharedtri.size() == 0)
							{
								ResultT[failtet[i]].set_label(tl);
								labeled = true;
								break;
							}
							else
							{
								surfRelateLabel.insert(tl);
								surfcount++;
							}

						}
					}
				}



				sharedtri_tmp.clear();
				sharedtri.clear();
				std::set_intersection(ResultV[ResultT[failtet[i]].Verts[1]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[1]].OrigBoundaryTris.end(), 
					ResultV[ResultT[failtet[i]].Verts[0]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[0]].OrigBoundaryTris.end(),back_inserter(sharedtri_tmp) );
				std::set_intersection(sharedtri_tmp.begin(),sharedtri_tmp.end(), 
					ResultV[ResultT[failtet[i]].Verts[3]].OrigBoundaryTris.begin(),ResultV[ResultT[failtet[i]].Verts[3]].OrigBoundaryTris.end(),back_inserter(sharedtri) );


				//================
				sharedtet_tmp.clear();
				sharedtet.clear();

				std::set_intersection(ResultV[ResultT[failtet[i]].Verts[0]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[0]].Tets.end(), 
					ResultV[ResultT[failtet[i]].Verts[3]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[3]].Tets.end(),back_inserter(sharedtet_tmp) );
				std::set_intersection(sharedtet_tmp.begin(),sharedtet_tmp.end(), 
					ResultV[ResultT[failtet[i]].Verts[1]].Tets.begin(),ResultV[ResultT[failtet[i]].Verts[1]].Tets.end(),back_inserter(sharedtet) );

				//std::cout<<"** shared tet num 4: "<<sharedtet.size()<<endl;
				if(!labeled)
				{
					for(int p = 0; p < sharedtet.size(); ++p)
					{
						//std::cout<<"shared tet: "<<Tets[sharedtet[p]]->get_label()<<endl;

						int tl = ResultT[sharedtet[p]].get_label();
						if (tl>=0)
						{
							if(sharedtri.size() == 0)
							{
								ResultT[failtet[i]].set_label(tl);
								labeled = true;
								break;
							}
							else
							{
								surfRelateLabel.insert(tl);
								surfcount++;
							}

						}
					}
				}

				if(!labeled)
				{

					for (int q = 0; q < sharedlabel.size(); q++ )
					{
						//cout<<sharedlabel[q]<<endl;

						std::set<int>::iterator it;
						bool find = true;

						//cout<<"surface related label size: "<< surfRelateLabel.size()<<endl;
						for ( it = surfRelateLabel.begin(); it != surfRelateLabel.end(); it++)
						{
							//cout<<"tl:  "<<*it<<endl;
							if(sharedlabel[q] == *it)
							{

								find = false;
								break;
							}
						}
						if(find && surfRelateLabel.size() >0)
						{

							ResultT[failtet[i]].set_label(sharedlabel[q]);
							labeled = true;


							break;
						}

					}

				}
			}


			if(!labeled)
			{
				failtet.push_back(failtet[i]);

			}


		}

		cout<<"fail: "<<failcount<<endl;
	}


}
