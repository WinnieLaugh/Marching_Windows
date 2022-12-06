#include "stdafx.h"
#include "simplification3d.h"
#include "Utils/MemorySurveillance.h"
using namespace std;
#include "H5Cpp.h"
using namespace H5;
#include <iostream>

cwg::Simplification3D::Simplification3D()
{
	m_labelvolume = NULL;
	m_tetmesh = NULL;
	is_double_padded_volume = false;
	cleaver_sparsity = 1.0;
	render_mode =  false;
}


cwg::Simplification3D::~Simplification3D()
{
	SAFE_DELETE(m_labelvolume);
	SAFE_DELETE(m_tetmesh);
}


void cwg::Simplification3D::GenerateSingleLabelData( int h, int w, int nframes )
{
	m_labelvolume = new cwg::PaddedTSPVideo();
	m_labelvolume->GenerateSimpleTSPVideo(h, w, nframes);
}


/*
Function to intialize the sliding.
*/
void cwg::Simplification3D::Init_sliding_wenhua(){
	
	m_tetmesh->TimerObj.tick();
	m_tetmesh->build_internal_boundaries_sliding_wenhua();

	m_tetmesh->build_edges_wenhua();

	m_tetmesh->compute_verts_Q_wenhua(m_tetmesh);
	m_tetmesh->compute_edge_contraction();
	m_tetmesh->build_heap();
}


void cwg::Simplification3D::Simplify( int tar_nverts )
{
	m_tetmesh->simplify(tar_nverts);
}


void cwg::Simplification3D::ImproveQuality(Visualize &vis,int odt_count)
{
	if (odt_count>0)
	{
		m_tetmesh->tetgenCDT(odt_count, false, 0, 0, 0);
	}

}


int cwg::Simplification3D::ImproveQuality_window(int odt_count, int starts_x, int  starts_y, int starts_z)
{

	return m_tetmesh->tetgenCDT(odt_count, true, starts_x, starts_y, starts_z);

}


void cwg::Simplification3D::LoadVolumeData( const std::string& foldername )
{
	if (is_double_padded_volume)
		m_labelvolume = new cwg::PaddedTSPVideo(4);
	else
		m_labelvolume = new cwg::PaddedTSPVideo(2);
	m_labelvolume->Load_from_Folder(foldername);

	m_tetmesh = new TetrahedralMesh();
	m_tetmesh->set_volume_data(m_labelvolume);
	m_tetmesh->set_input_foldername(foldername);

	if (foldername.find("block") != std::string::npos)
	{
		string vec_filename = foldername + "shiftvec.txt";
		ifstream ifile(vec_filename);
		if (ifile)
		{
			int x, y, z;
			ifile >> x >> y >> z;
			shift_vec = Cleaver::vec3(x,y,z);
		}
	}

	//m_tetmesh->TimerObj.start();
}


void cwg::Simplification3D::LoadMeshData( const std::string& filename )
{


	m_tetmesh = new TetrahedralMesh();
	m_tetmesh->set_volume_data(NULL);
	m_tetmesh->set_input_filename(filename);
	m_tetmesh->TimerObj.start();

	cout << "Loading Mesh File... \t";
	m_tetmesh->load_bin_tetm(filename);
	cout << "cost " << m_tetmesh->TimerObj.tick() << " sec." << endl;
}


void cwg::Simplification3D::ExportData()
{

	m_tetmesh->build_triangles();

	m_tetmesh->translate(shift_vec);

	char str[300];
	string foldername = m_tetmesh->get_foldername();

	//cout<<"====================="<<foldername<<endl;

	cout<<"Vertex number:  "<<m_tetmesh->ResultV.size()<<endl;

	sprintf(str, "final_%3.4lf_v%06d_t%06d.tetm.ply", m_tetmesh->get_final_percentage(),
		m_tetmesh->ResultV.size(), m_tetmesh->ResultT.size());
	cout << "Now Saving... May be Slow!" << endl;
	m_tetmesh->TimerObj.tick();
	m_tetmesh->save_bin_tetm(foldername+str+".tetm");
	m_tetmesh->save_txt_vol_ply_weight(foldername+str);

	sprintf(str, "final_%3.4lf_v%06d_t%06d_no0.tetm.ply", m_tetmesh->get_final_percentage(),
		m_tetmesh->ResultV.size(), m_tetmesh->ResultT.size());
	m_tetmesh->save_txt_vol_ply_without0_weights(foldername+str);


	sprintf(str, "final_%3.4lf_v%06d_t%06d.surf.ply", m_tetmesh->get_final_percentage(),
		m_tetmesh->ResultV.size(), m_tetmesh->ResultT.size());
	m_tetmesh->save_txt_surf_ply(foldername+str);
	cout << "Saving File: ";
	cout << "cost " << m_tetmesh->TimerObj.tick() << " sec." << endl;

	m_tetmesh->split_mesh("submesh\\");
	PRINT_TIME;
	cout << "DONE!" << endl;


}


void cwg::Simplification3D::ExportData(const std::string& filename)
{

	m_tetmesh->build_triangles();

	m_tetmesh->translate(shift_vec);

	char str[300];

	m_tetmesh->save_txt_vol_ply_weight(filename+".ply");
	m_tetmesh->save_bin_tetm(filename+".tetm");

	//m_tetmesh->save_ascii_cgal_from_result (filename+".ascii.cgal");

	cout << "Saving File: ";
	cout << "cost " << m_tetmesh->TimerObj.tick() << " sec." << endl;

	//m_tetmesh->split_mesh("submesh\\");
	PRINT_TIME;
	cout << "DONE!" << endl;


}


void cwg::Simplification3D::DensityField( const std::string& foldername )
{
	//m_tetmesh->TimerObj.start();

	if (is_double_padded_volume)
		m_labelvolume = new cwg::PaddedTSPVideo(4);
	else
		m_labelvolume = new cwg::PaddedTSPVideo(2);
	m_labelvolume->Load_from_Folder(foldername);

	m_tetmesh = new TetrahedralMesh();
	m_tetmesh->set_volume_data(m_labelvolume);
	m_tetmesh->set_input_foldername(foldername);

	if (foldername.find("block") != std::string::npos)
	{
		string vec_filename = foldername + "shiftvec.txt";
		ifstream ifile(vec_filename);
		if (ifile)
		{
			int x, y, z;
			ifile >> x >> y >> z;
			shift_vec = Cleaver::vec3(x,y,z);
		}
	}

	//m_tetmesh->TimerObj.start();
}


void cwg::Simplification3D::SlidingSimplifyMultiProcess (Visualize &vis, const std::string& foldername, TetrahedralMesh& m_tetmesh, double tar_percentage, const int slides[3], const int sizes[3])
{
	int half_slides[3] = {slides[0]/2, slides[1]/2, slides[2]/2};

	int starts[3] ={0,0,0};

	int targetv = (slides[0]+1)*( slides[1]+1)*(slides[2]+1)*tar_percentage;
	vector<string> filelist = show_file_list(foldername, "svdata");

	for(int i=0; i<3; i++){
		m_tetmesh.steps[i] = slides[i]/2;
	}	

	string& data_folder_name = foldername + "\\mesh\\Data";
	mkdir(data_folder_name.c_str());

	m_tetmesh.compute_component_branch_size = 0;
	int odt_success = -1;

	time_t simplify_time_start, simplify_time_end, io_time_start, io_time_end;
	double simplify_time=0, io_time=0;

	while(starts[2] + half_slides[2]< sizes[2]){ //z 
		if(starts[2] + half_slides[2] + m_tetmesh.steps[2] > sizes[2])
			m_tetmesh.steps[2] = sizes[2] - half_slides[2] - starts[2];
		
		//a new face
		starts[1] = 0;
		while(starts[1] + half_slides[1]< sizes[1]){ // y
			if(starts[1] + half_slides[1] + m_tetmesh.steps[1] > sizes[1])
				m_tetmesh.steps[1] = sizes[1] - half_slides[1]- starts[1];
			
			//a new line
			starts[0] = 0;
			while(starts[0] + half_slides[0]<sizes[0]){ //x
				if(starts[0] + half_slides[0] + m_tetmesh.steps[0] > sizes[0])
					m_tetmesh.steps[0] = sizes[0] - half_slides[0] -starts[0];

				//do the initialization
				std::cout<<"=====================  win ===================  \n ";			
				m_tetmesh.slide_move_tsize = (half_slides[0]+m_tetmesh.steps[0]) * (half_slides[1]+m_tetmesh.steps[1]) * (half_slides[2]+m_tetmesh.steps[2]) * m_tetmesh.tetn;
				m_tetmesh.slide_move_vsize = (half_slides[0]+m_tetmesh.steps[0] + 1) * (half_slides[1]+m_tetmesh.steps[1] + 1) * (half_slides[2]+m_tetmesh.steps[2] + 1);
				targetv = m_tetmesh.slide_move_vsize*tar_percentage;

				int endx_now = starts[0] + half_slides[0] + m_tetmesh.steps[0];
				int endy_now = starts[1] + half_slides[1] + m_tetmesh.steps[1];
				int endz_now = starts[2] + half_slides[2] + m_tetmesh.steps[2];

				time(&io_time_start);
				
				m_tetmesh.loadVolumeData_sliding_wenhua(vis, foldername, filelist, starts[0], starts[1], starts[2], endx_now, endy_now, endz_now);
				Init_sliding_wenhua();

				Simplify(targetv);

				//m_tetmesh.getResultV_T_window_test(vis);
				//odt_success = ImproveQuality_window(odt_count, starts[0], starts[1], starts[2]);

				//if (odt_success == 0)
				//	m_tetmesh.vert_adjust();

				m_tetmesh.tet_Adjust();
				m_tetmesh.keep_Adjust( endx_now == m_tetmesh.m_X );
				time(&simplify_time_end);
				simplify_time += difftime(simplify_time_end, simplify_time_start);

				time(&io_time_start);
				m_tetmesh.saveData(starts[0], starts[1], starts[2]);
				time(&io_time_end);
				io_time += difftime(io_time_end, io_time_start);

				if(endx_now == m_tetmesh.m_X)
					m_tetmesh.reset_all();

				starts[0] += m_tetmesh.steps[0];
				m_tetmesh.steps[0] = slides[0]/2;
			}

			starts[1] += m_tetmesh.steps[1];
			m_tetmesh.steps[1] = slides[1]/2;
		}
		starts[2] += m_tetmesh.steps[2];
		m_tetmesh.steps[2] = slides[2]/2;
	}

	ofstream ofile("time_splits.txt");
	ofile << "io time used: "<< io_time << std::endl;
	ofile << "simplify time used: "<< simplify_time << std::endl;
	ofile.close();

}


void cwg::Simplification3D::SlidingSimplify(Visualize &vis, const std::string& foldername, double tar_percentage, double edge_var_thr, int slidex, int slidey, int slidez)
{
	std::cout << "start sliding !" << std::endl;

	m_tetmesh = new TetrahedralMesh();
	m_tetmesh->set_input_foldername(foldername);
	m_tetmesh->TimerObj.start();

	m_tetmesh->edge_var_thr = edge_var_thr;
	m_tetmesh->tetn = 5;

	vector<string> filelist = show_file_list(foldername, "svdata");

	string& filename = foldername+filelist[0];
	ifstream ifile;
	ifile.open( filename.c_str() , ios::binary );

	int iRows = 0;
	int iCols = 0;
	ifile.read( (char*)&iRows , sizeof(int) );
	ifile.read( (char*)&iCols , sizeof(int) );

	ifile.close();
	m_tetmesh->m_X = iRows;
	std::cout << "x: " << iRows << std::endl;
	m_tetmesh->m_Y = iCols;
	std::cout << "y: " << iCols << std::endl;
	m_tetmesh->m_Z = filelist.size();
	std::cout << "z: " << filelist.size() << std::endl;

	//along the max axis to sliding
	int x_start = 0, y_start = 0, z_start = 0;
	int dir_cnt = 0; //direction count
	m_tetmesh->directions[0] = true;
	m_tetmesh->directions[1] = true;
	m_tetmesh->directions[2] = true;

	m_tetmesh->init_memory(slidex, slidey, slidez);
	m_tetmesh->slide_X = slidex;
	m_tetmesh->slide_Y = slidey;
	m_tetmesh->slide_Z = slidez;

	m_tetmesh->win = 0;	//the number of window movement
	m_tetmesh->direction = 0;
	int targetv = (slidex+1)*(slidey+1)*(slidez+1)*tar_percentage;

	int slides[3] = {slidex, slidey, slidez};
	int sizes[3] = {iRows, iCols, filelist.size()};
	
	SlidingSimplifyMultiProcess (vis, foldername, *m_tetmesh, tar_percentage, slides, sizes);

	std::cout << "=============================================================" << std::endl;
	std::cout << "===============Finished Sliding! Now combining===============" << std::endl;
	std::cout << "=============================================================" << std::endl;

	m_tetmesh->clear_memory();
	
	m_tetmesh->getResultV_T();

}


void cwg::Simplification3D::ResetID()
{
	m_tetmesh->reset_id();
}