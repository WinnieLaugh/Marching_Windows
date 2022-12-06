#include "toolbox.h"


using namespace std;
using namespace cwg;
using namespace cv;

#include <string>
#include "H5Cpp.h"

using namespace H5;

void print_help()
{
	cout << "How-to-Run: SimpTetMesh.exe PathName TargetPecent [P2(0)/P4(1)] [Density]" << endl;
}

void ChangeColor_density(char* inputFile, char* objFile)
{
	//read
	ifstream fin(inputFile);
	if (!fin)
	{
		exit(1);
	}

	ofstream fout(objFile);
	if (!fout)
	{
		exit(2);
	}

	string s="";
	int num;
	double vx,vy,vz;
	double r,g,b;
	int label;
	int vIdx[4];

	double maxx=100,maxy=100,maxz=15;

	while(fin>>s)
	{

		if(s=="v")
		{
			fin>>vx;
			fin>>vy;
			fin>>vz;

			fin>>r;
			fin>>g;
			fin>>b;

			if(r == 0)
			{
				fout<<"v "<<vx<<" "<<vy<<" "<<vz<<" "<<0<<" "<<1<<" "<<g<<endl;
			}
			else if(b == 0)
			{
				fout<<"v "<<vx<<" "<<vy<<" "<<vz<<" "<<1<<" "<<1-r<<" "<<0<<endl;
			}

		}
		else
		{
			fout<<s<<" ";
			getline(fin,s);
			fout<<s<<endl;
		}
	}

	fin.close();
	fout.close();
}

void ChangeColor_metro(char* inputFile, char* objFile)
{
	//read
	ifstream fin(inputFile);
	if (!fin)
	{
		exit(1);
	}

	ofstream fout(objFile);
	if (!fout)
	{
		exit(2);
	}

	string s="";
	int num;
	double vx,vy,vz;
	double r,g,b;
	int label;
	int vIdx[4];

	double maxx=100,maxy=100,maxz=15;

	while(fin>>s)
	{

		if(s=="v")
		{
			fin>>vx;
			fin>>vy;
			fin>>vz;

			fin>>r;
			fin>>g;
			fin>>b;


			fout<<"v "<<vx<<" "<<vy<<" "<<vz<<" "<<b<<" "<<g<<" "<<r<<endl;
		}
		else
		{
			fout<<s<<" ";
			getline(fin,s);
			fout<<s<<endl;
		}
	}

	fin.close();
	fout.close();
}


void TetQuality_R(string filename){

	std::cout << "check quality Tet path: " << filename << std::endl; 

	TetrahedralMesh *tet = new TetrahedralMesh();
	std::cout << filename << std::endl;
	tet->load_bin_tetm(filename.c_str());
	std::cout << "finished loading" << std::endl;

	//get all ratio
	cwg::Tetrahedron t;
	t.compute_sphere_ratio(tet);

	double r_ratio;
	double minr=2,maxr=0;
	double PI=3.141592653;
	double regular_angle = 7*PI/18.0;	//arccos (1/3)  70 degree,   radian = 1.22173
	double sum = 0;
	double *ratio = new double[200000000];

	int sliver = 0;
	int invalid = 0;
	for (int i=0;i<tet->m_storage_tets.size(); i++)
	{
		cwg::Tetrahedron *t = &tet->m_storage_tets.at(i);
		r_ratio = t->get_ratio();

		if (r_ratio < 0.1 )
		{
			++sliver;
		}
		ratio[i] = r_ratio;

		if(r_ratio > -0.1 && r_ratio < 2.1){
			sum += r_ratio;
		}else{
			invalid ++;
		}

		if(r_ratio<minr)
			minr=r_ratio;
		if(r_ratio>maxr)
			maxr=r_ratio;

	}

	float range = maxr - minr;

	std::cout << "num of verts: " << tet->m_storage_verts.size() << std::endl;
	std::cout << "Sliver: %" << (double)sliver/tet->m_storage_tets.size() <<std::endl;
	std::cout << "Min ratio: " << minr << std::endl;
	std::cout << "Average ratio" << sum/ (tet->m_storage_tets.size()-invalid)<<std::endl;
	std::cout << "invalid: " << invalid << std::endl;

	int stat[200];
	for(int i=0; i<200;i++)
	{
		stat[i]=0;
	}
	for (int i=0; i<tet->m_storage_tets.size(); i++)
	{
		for (int j=0; j<200; ++j)
		{
			if(ratio[i] >= minr+range/200*j && ratio[i] < minr+range/200*(j+1))
			{
				stat[j]++;
				break;
			}
		}

	}

	tet->build_triangles();
	//tet->save_txt_vol_ply_quality(outfile2.c_str(),mina,maxa);
	//tet->save_obj_quality(outfile2);
}

void readsvdata(const std::string& foldername)
{
	vector<string> filelist = show_file_list(foldername, "svdata");
	int label;
	//for (int k = 0; k < filelist.size(); ++k)
	{
		string& filename = foldername+filelist[1];
		ifstream ifile;
		ifile.open( filename.c_str() , ios::binary );
		//ifile.seekg((3 )*sizeof(int));
		int iRows = 0;
		int iCols = 0;
		ifile.read( (char*)&iRows , sizeof(int) );
		ifile.read( (char*)&iCols , sizeof(int) );
		for (int i = 0; i < iRows; ++i)
		{
			for (int j = 0; j < iCols; ++j)
			{
				ifile.read( (char*)&label , sizeof(int) );
				cout<<label<<" ";

			}
			cout<<endl;
		}


	}
}


void main(int argc, char** argv)
{
	time_t rawtime;
	struct tm * timeinfo;
	char time_output[20];
	ofstream ofile("time.txt");

	time(&rawtime);
	
	timeinfo = localtime(&rawtime);
	printf("Current local time and data: %s", asctime(timeinfo));
	strftime(time_output, 20, "%d/%m/%Y %H:%M:%S", timeinfo);
	std::cout << time_output << std::endl;
	ofile << "start time: " << time_output << std::endl;

	string pathname = argv[1];
	tetmesh_simplify(argc, argv);

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current local time and data: %s", asctime(timeinfo));
	strftime(time_output, 20, "%d/%m/%Y %H:%M:%S", timeinfo);
	ofile << "end time: " << time_output << std::endl;
	ofile.close();
	
	TetQuality_R("after_improve.tetm");

	TetrahedralMesh *tet = new TetrahedralMesh();
	std::cout << "after_improve.tetm" << std::endl;
	tet->load_bin_tetm_to_result("after_improve.tetm");
	tet->save_mesh_from_result("after_improve_vert.mesh");

	int test;
	std::cout << "finished! " << std::endl;
	std::cin >> test;

	system("pause");

}

