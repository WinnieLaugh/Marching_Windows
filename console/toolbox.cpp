#include "toolbox.h"
#include <process.h>

using namespace std;
using namespace cv;
using namespace cwg;

extern SIMPLIFICATION3D_API bool glb_boundary_fixed;
Visualize vis;

void convert2svdata(int argc, char** argv)
{
	string filename = argv[1];
	string foldername = argv[2];

	ifstream ifile;
	ifile.open( filename.c_str() , ios::binary );
	if (!ifile.is_open())
		return;

	vector<cv::Mat> labelvolume;

	int Z = -1;
	int Y = -1;
	int X = -1;
	ifile.read( (char*)&X, sizeof(int) );
	ifile.read( (char*)&Y, sizeof(int) );
	ifile.read( (char*)&Z, sizeof(int) );

	

	labelvolume.resize(Z);
	int PixelsPerFrame = X*Y;
	for (int k=0; k<Z; k++)
	{
		Mat kMatrix(Y, X, CV_32SC1);
		char* pbuff = (char*) kMatrix.data;
		ifile.read( pbuff, sizeof(int)*PixelsPerFrame );
		labelvolume[k] = kMatrix;
	}

	ifile.close();

	mkdir(foldername.c_str());
	char str[300];
	double ratio = 1.0;
	for (int i=0; i<labelvolume.size(); i++)
	{
		Mat tmp;
		cv::resize(labelvolume[i], tmp, cv::Size(), ratio, ratio, INTER_NEAREST);
		sprintf(str, "frame%04d.svdata", i);
		writemat(foldername+str, tmp);
	}
}



void split_whole_data_into_blocks( int argc, char** argv )
{
	try
	{
		string foldername = argv[1];
		int X = atoi(argv[2]);
		int Y = atoi(argv[3]);
		int Z = atoi(argv[4]);

		TSPVideo the_full_tsp_video;
		the_full_tsp_video.Load_from_Folder(foldername);

		int W = the_full_tsp_video.Width();
		int H = the_full_tsp_video.Height();
		int F = the_full_tsp_video.nFrames();

		char str[300];
		string output_foldername;

		ofstream ofileBlocks(foldername+"BlockList.txt");
		for (int k=0; k<=F/Z; k++)
		{
			int ks = k * Z;
			int ke = min( (k+1)*Z, F );

			if (ke <= ks)
				continue;

			for (int i=0; i<=H/Y; i++)
			{
				int is = i * Y;
				int ie = min( (i+1)*Y, H );

				if (ie <= is)
					continue;

				for (int j=0; j<=W/X; j++)
				{
					int js = j * X;
					int je = min( (j+1)*X, W );

					if (je <= js)
						continue;

					sprintf(str, "block%02d_%02d_%02d_%04d_%04d_%04d\\", k, i, j, k*Z, i*Y, j*X);
					ofileBlocks << str << endl;
					output_foldername = foldername + str;
					mkdir(output_foldername.c_str());
					ofstream ofile(output_foldername+"shiftvec.txt");
					ofile << j*X << " " << i*Y << " " << k*Z;
					ofile.close();
					save_tsp_blocks(the_full_tsp_video, ks, ke, is, ie, js, je, output_foldername);
				}
			}
		}
		ofileBlocks.close();
	}
	catch (exception e)
	{
		cout << e.what() << endl;
	}
}

void save_tsp_blocks( const cwg::TSPVideo& video, int ks, int ke, int is, int ie, int js, int je, const std::string& output_foldername )
{
	char str[300];
	for (int k=ks; k<ke; k++)
	{
		cv::Mat tmp = video.frame(k);
		cv::Mat tmppart = tmp(cv::Range(is,ie), cv::Range(js,je));

		sprintf(str, "%03d.svdata", k);
		string filename = output_foldername + str;
		writemat(filename, tmppart);
	}
}

void tetmesh_simplify( int argc, char** argv )
{
	//input parameter: foldername, sim_percentage, odt_count, edge variation thr, slide size_x, slide size_y, slide size_z(0 means full size)
	// slide size should be even, easy to move
	try 
	{	
		string foldername = argv[1];
		glb_boundary_fixed = false;
		double percentage = atof(argv[2]);

		int tar_verts = 0;

		int odt_count = atoi(argv[3]);

		Simplification3D simpobj;

		double edge_var_thr = atof(argv[4]);

		int slidex = atof(argv[5]);
		int slidey = atof(argv[6]);
		int slidez = atof(argv[7]);

		time_t rawtime;
		struct tm * timeinfo;
		char time_output[20];
		ofstream ofile("time_cost.txt");

		time(&rawtime);
		timeinfo = localtime(&rawtime);
		printf("before simplify: %s", asctime(timeinfo));	
		strftime(time_output, 20, "%d/%m/%Y %H:%M:%S", timeinfo);
	    ofile << "before simplify: "<< time_output << std::endl;
		simpobj.odt_count = odt_count;

		vis.setFoldername(foldername);
		simpobj.SlidingSimplify(vis, foldername,percentage, edge_var_thr, slidex, slidey, slidez);

		time(&rawtime);	
		timeinfo = localtime(&rawtime);
		printf("After improvement time:: %s", asctime(timeinfo));
		strftime(time_output, 20, "%d/%m/%Y %H:%M:%S", timeinfo);
	    ofile << "after_improve: "<< time_output << std::endl;
		simpobj.ExportData("after_improve");

		simpobj.ExportData();
		simpobj.ResetID();

	}
	catch (exception e)
	{
		cout << e.what() << endl;
	}
}
 
void improve_tetm( int argc, char** argv )
{
	try 
	{
		
		string filename = argv[1];
		int odt_count = atoi(argv[2]);

		Simplification3D simpobj;
		
		simpobj.LoadMeshData(filename);
		simpobj.ImproveQuality(vis, odt_count);
		simpobj.ExportData(filename);
		

	}
	catch (exception e)
	{
		cout << e.what() << endl;
	}
}

void merge_tetm( int argc, char** argv )
{
	string foldername = argv[1];
	ifstream ifile(foldername+"BlockList.txt");
	if (!ifile)
		return;

	char str[300];
	vector<TetrahedralMesh*> blockmeshes;
	while (ifile.good())
	{
		ifile.getline(str, 256);
		string tmp(str);
		if (tmp.size() == 0)
			break;
		string subfoldername = foldername + str + "mesh\\";
		vector<string> filelists = show_file_list(subfoldername, "tetm");
		string filename_tetm = subfoldername+filelists[0];

		TetrahedralMesh* tetmesh = new TetrahedralMesh();
		tetmesh->load_bin_tetm(filename_tetm);

		blockmeshes.push_back(tetmesh);
	}

	TetrahedralMesh* wholemesh = new TetrahedralMesh();
	wholemesh->set_input_foldername(foldername);
	wholemesh->merge(blockmeshes);
	wholemesh->save_bin_tetm(foldername+"wholemesh.tetm");

	wholemesh->build_triangles();
	wholemesh->save_txt_vol_ply(foldername+"wholemesh.tetm.vol.ply");
	wholemesh->save_txt_surf_ply(foldername+"wholemesh.tetm.surf.ply");
	wholemesh->split_mesh("submesh\\");

	for (int i=0; i<blockmeshes.size(); i++)
		SAFE_DELETE(blockmeshes[i]);
	SAFE_DELETE(wholemesh);
}

void resize_tetm( int argc, char** argv )
{
	if (argc != 6)
		return;

	string inputmesh_filename = argv[1];
	string outputmesh_filename = argv[2];
	double x = atof(argv[3]);
	double y = atof(argv[4]);
	double z = atof(argv[5]);
	TetrahedralMesh *tm = new TetrahedralMesh();
	tm->load_bin_tetm(inputmesh_filename);
	tm->resize(x,y,z);
	tm->save_bin_tetm(outputmesh_filename);
	tm->build_triangles();
	tm->save_txt_surf_ply(outputmesh_filename+".surf.ply");
	tm->save_txt_vol_ply(outputmesh_filename+".vol.ply");
	tm->split_mesh("submesh\\");

	SAFE_DELETE(tm);
}
