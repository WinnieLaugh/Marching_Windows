#ifndef _TETRAHEDRON_MESH_H_
#define _TETRAHEDRON_MESH_H_

#ifdef SIMPLIFICATION3D_EXPORTS
#define SIMPLIFICATION3D_API __declspec(dllexport)
#else
#define SIMPLIFICATION3D_API __declspec(dllimport)
#endif

#include <cwgUtilities.h>

#include <Visualize.h>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>

#include "Cleaverlib/Cleaver.h"
#include "TspVideo.h"
#include "Tetrahedron.h"
#include "TetEdgeHeap.h"
#include "TetCube.h"
#include "TriangularMesh.h"

#include <CGAL/global_functions_circular_kernel_2.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#define TETLIBRARY 
#include <tetgen.h>

//Visualize include
#include <Vert.h>

extern SIMPLIFICATION3D_API bool glb_boundary_fixed;
const int index_size = 2*2*2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;


namespace cwg
{
	class SIMPLIFICATION3D_API TetrahedralMesh
	{
	public:
		TetrahedralMesh();
		~TetrahedralMesh();

		void								tetrahedralize_cube();
		void								release();
		void								copy_to(TetrahedralMesh* temp_mesh) const;
		void								copy_submesh_to(TetrahedralMesh* temp_mesh, int lb) const;
		void								reset_id();
		void								merge(const std::vector<TetrahedralMesh*>& submeshes);

		//set weight -- get density field based on volume points
		void								densityField();		
		void								compute_lfs(std::queue<int>& surfV, int vsize);	//return the coefficient of surface v, coefficient*lfs = weight
		
		//========= compute weight to smooth from edge variation ==========
		void								compute_weight_from_edge_variation();

		// Simplification:
		void								build_internal_boundaries();
		void								compute_edge_contraction();
		void								build_heap();
		void								simplify(int tar_nverts);
		void								contract(TetEdge* e);

		
		// odt smooth
		void								ODTsmooth(bool first_smooth);
		void								ODTsmoothSurfOnly();
		void								ODT_tet(int i, Cleaver::vec3 &newpos, double &neighbor_num);
		void								ODT_tet_no_weight(int i, Cleaver::vec3 &newpos, double &neighbor_num);
		void								ODT_tri(int i, Cleaver::vec3 &newpos, double &neighbor_num);
		void								ODT_interface(int i, Cleaver::vec3 &newpos, double &neighbor_num);
		void								ODT_cube_edge(int i, Cleaver::vec3 &newpos, double &neighbor_num);
		void								ODT_curve(int i, Cleaver::vec3 &newpos, double &neighbor_num);

		//using tetgen to do constrained DT
		int									tetgenCDT(int odt_count, bool window_optimize, int starts_x, int starts_y, int starts_z);

		bool								improveSurfaceQ();	// flip and split surface edge to improve quality

		//void								volumeODT(const std::vector<int> surface, const std::vector<int> vertices, int label, std::vector<Tetrahedron*> &newTet);
		void								volumeODT(bool last);
		int								    volumeODT_window(bool last, int starts_x, int starts_y, int starts_z);

		void								split_mesh(const std::string& foldername);
		// Geo Processing:
		void								build_triangles();
		inline void							resize(double x, double y, double z) { resize(Cleaver::vec3(x,y,z)); }
		void								resize(const Cleaver::vec3& sxsysz);
		void								translate(const Cleaver::vec3& v);

		// I/O:
		bool								save_txt_tetm(const std::string& filename) const;
		bool								load_txt_tetm(const std::string& filename);
		bool								save_bin_tetm(const std::string& filename) const;
		bool								save_result_v_t(const std::string& filename_verts, const std::string& filename_tets);
		bool								save_result_v_t_label(const std::string& filename_verts, const std::string& filename_tets, int label_this);

		
		bool								load_bin_tetm(const std::string& filename);
		bool								load_bin_tetm_to_result(const std::string& filename);

		bool								load_txt_obj(const std::string& filename);

		bool								load_mesh(const std::string& filename);
		bool								load_mesh_surf(const std::string& filename);	//load surface

		bool								save_txt_vol_ply(const std::string& filename) const;
		bool								save_txt_vol_ply_without0_weights(const std::string& filename) const;
		bool								save_txt_vol_ply_quality(const std::string& filename, float min_angle, float max_angle) const;
		bool								save_txt_surf_ply(const std::string& filename) const;
		bool								save_obj_quality(const std::string& filename) const;
		
		bool								save_txt_vol_stellar_mesh(const std::string& filename) const;
		bool								save_txt_vol_ply_weight(const std::string& filename) const;
		bool								save_ascii_cgal(const std::string& filename);
		bool								save_ascii_cgal_from_result(const std::string& filename);
		bool								save_mesh_from_result(const std::string& filename);
		bool								save_mesh(const std::string& filename);

		inline void							set_volume_data(LabelVolume* labelvol) { m_labelvolume = labelvol; }
		inline void							set_input_filename(const std::string& filename) { m_foldername = filedir(filename); mkdir(m_foldername.c_str()); }
		inline void							set_input_foldername(const std::string& foldername) { m_foldername = foldername+"mesh\\"; mkdir(m_foldername.c_str()); }
		inline void							get_size(int sz[3]) const { sz[0] = m_X; sz[1] = m_Y; sz[2] = m_Z; }
		inline int							get_tar_verts() const { return m_tar_verts; }
		inline std::string					get_foldername() const { return m_foldername; }
		inline double						get_final_percentage() const { return m_final_percentage; }
		
		// Test Code:
		void								test(double val);
		void								test2(double val);

		//============== sliding ===============================
		void								init_memory(int slidex, int slidey, int slidez);
		void								clear_memory();

		//==================================wenhua test sliding================================
		void								walls(Visualize& vis, int startx, int starty, int startz, int endx, int endy, int endz);
		void								sliding_boundaries(Visualize& vis, int startx, int starty, int startz, int endx, int endy, int endz);
		void								sliding_tet_labels(Visualize& vis, const std::string& foldername , const std::vector<std::string>& filelist, int startx, int starty, int startz, int endx, int endy, int endz);
		void								loadVolumeData_sliding_wenhua(Visualize& vis, const std::string& foldername , const std::vector<std::string>& filelist, int startx, int starty, int startz, int endx, int endy, int endz);
		void								saveResultData_wenhua(int data_buffer_startx, int data_buffer_starty, int data_buffer_startz, int buffer_win_x, int buffer_win_y, int buffer_win_z);
		void								reset_win_verts(int win_inside_index, int startx, int starty, int startz);
		void								reset_all();

		//load from file, right, text
		void								loadData_fromFile_wenhua(const std::string& foldername , const std::vector<std::string>& filelist, int Data_cube_start_x, int Data_cube_start_y, int Data_cube_start_z, int win_inside_index);
		void								loadData_fromRight_wenhua(int startx, int starty, int startz);
		void								loadData_fromText_wenhua(int win_x, int win_y, int win_z, int win_inside_index);

		//load to buffer
		void								saveData(int data_cube_start_x, int data_cube_start_y, int data_cube_start_z);
		void								loadTet_tmp(int start_x, int start_y, int start_z);
			
		void								loadData_betweenCubes_wenhua();
		void								keep_Adjust(bool x_end);
		void								tet_Adjust();
		void								vert_adjust();
		bool								tet_boundary_check(int tet_index);
		bool								vert_boundary_check(int vert_id);
		bool								result_tris_boundary_check(int vert0, int vert1, int vert2);
		int									find_available_tet_id(int orig_tet_x, int orig_tet_y, int orig_tet_z);
		int									find_available_tmp_tet_id(int orig_tmp_tet_x, int orig_tmp_tet_y, int orig_tmp_tet_z);

		int									start_index[index_size][4];
		int									find_adjust_tet_id(int vert[4][3]);
		////////////////////////////*********** DEBUG ************************////////////////////////
		void								getResultV_T();
		bool								to_debug_bool;
		void								compute_component_nums();
		void								compute_component_nums_result();
		int									compute_component_branch_size;
		void								getResultV_T_window_test(Visualize& vis);
		void								storage_to_result();
		void								write_obj_for_debug(int tet_id);
		void								write_obj_for_debug2(int x, int y, int z);
		////////////////////////////*********** DEBUG ************************////////////////////////
		
		std::map<int, int>					vert_map;
		std::map<int, int>					vert_window_map;
		std::map<int, int>					tet_window_map;

		std::vector<int>					vert_window_map_inverse;
		bool								keep_in_previous;
		int									result_v_id;
		int									result_t_id;
		bool								render_mode;
		
		std::vector<Tetrahedron>			m_storage_tets_tmp; // to store the tetradral which are not fully inside the cube, to see whether other cubes have it's lost verts.
		
		//===================================format=============================================
		void								build_internal_boundaries_sliding_wenhua();

		//===================================wenhua test simplification based on YAting's work================
		void								build_edges_wenhua();
		void								compute_verts_Q_wenhua(TetrahedralMesh* tetmesh);

		//======= sliding sim  start point ===========
		
		int									tri_start;	//new slide start from this surface triangle

		int									result_t_start;	// for updating result t's v
		int									step;	//sliding window - step size
		int									steps[3];
		int									load_vidx;	//the v index to load v
		int									save_vidx;	//save v
		int									load_tidx;	//the t index to load v
		int									save_tidx;	//save t
		int									slide_X, slide_Y, slide_Z;
		int									slide_move_vsize; //the number of update v during move
		int									slide_move_tsize; //the number of update t during move
		int									m_X, m_Y, m_Z;	//volume data size
		int									win;	//the number of window
		int									direction;	//window move direction
		bool								process_all;	//if process all at once

		int									tetn; //the number of tets of each cube

		std::vector<TetVertex>				ResultV;	//store the simplified result v
		std::vector<Tetrahedron>			ResultT;	//store the simplified result v

		std::vector<TetVertex>				ResultV_tmp; // store the tet id of each vert before ODT simplification
		std::vector<Tetrahedron>			ResultT_tmp; // store the vert id before ODT simplification
		//============ edge variation thr===============
		double								edge_var_thr;

		// Main Data Structures:
		std::vector<Tetrahedron>			m_storage_tets;
		std::vector<TetVertex>				m_storage_verts;
		std::vector<TetTriangle*>			Tris;
		TriangularMesh*						InternalMesh;
		std::set< cwg::ordered_pair >		EdgeSet;
		std::map< cwg::ordered_triple, std::vector<int> > triple_map;
		

		Sparse_Matrix*						AdjMatrix;
		cwg::TetEdgeHeap					EdgeHeap;
		
		Timer								TimerObj;
		Timer								heap_opt_time;	//heap operation time
		double								heap_time;


		//***********************************wenhua****************************
		int									walls_begin[3]; //to record if the cube to load is on the boundary of the whole data, i.e. startx == 0
		int									walls_end[3]; //to record if the cube to load is at the end of the big data, i.e. endx == m_Y
		bool								directions[3];	
		
		//********************************************************************************

	private:
		void								remove_border_tets();
		void								flatten_boundary_tets();
		void								remove_degenerated_tets();
		void								remove_nonreferenced_verts();
		void								translate_model();
		void								initial_sortout();
		void								compute_verts_border_code();

		void								compute_edge_tets(const TetVertex* v0, const TetVertex* v1, std::vector<int>& edge_tets) const;
		void								compute_edge_tets(TetEdge* e, std::vector<int>& edge_tets) const;
		void								remove_tets(const std::vector<int>& edge_tets);
		void								remove_edge(TetVertex* v0, TetVertex* v1, TetEdge* e);
		void								remove_edge_sliding(TetVertex* v0);
		void								remove_edge_sliding(TetVertex* v0, TetVertex* v1);
		void								replace_verts_in_v1_tets(TetVertex* v0, TetVertex* v1);
		void								replace_verts_in_v1_edges(TetVertex* v0, TetVertex* v1);
		void								update_edges_of_vert(TetVertex* v);
		int									compute_valid_nverts();
		int									compute_valid_ntets();
		void								private_tets(TetVertex* v0, TetVertex* v1, std::vector<int>& v0pvtets, std::vector<int>& v1pvtets);

		// functions may be removed:
		void								contract_old(TetEdge* e);
		void								tetrahedronlize(int i, int j, int k, TetCube& tmp_cube, std::vector<Tetrahedron>::iterator& it, int& tID);
		

	private:
		Cleaver::TetMesh*					m_cleaver_tet_mesh;		
		LabelVolume*						m_labelvolume;

		std::string							m_foldername;
		int									m_tar_verts;
		double								m_final_percentage;
		
		std::map<int, std::vector<int>>		m_submesh_hashtable;

	private:
		std::vector<TetTriangle>			m_storage_tris; //copy_to
		std::vector<TetEdge>				m_storage_edges;
		std::vector<TetTriangle*>			tmp_tri_mesh;

		std::vector<int>					last_slideboundary_v;	//the last line of the slide boundary
	
	};


	void									write_ply_file_head(std::ofstream& ofile, int vertex_count, int face_count);
	
}

#endif