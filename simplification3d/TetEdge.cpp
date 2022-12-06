#include "stdafx.h"
#include "TetEdge.h"
#include "TetrahedralMesh.h"
#include "Tetrahedron.h"

using namespace std;
using namespace cv;
int cwg::TetEdge::ms_nextid = 0;
int debug_num=0;
double x_0 = 3.0; double y_0 = 4.0; double z_0 = 0.0;
double x_1 = 2.0; double y_1 = 4.0; double z_1 = 0.0;

cwg::TetEdgeAroundInfo cwg::TetEdge::ms_edge_around_info;
extern SIMPLIFICATION3D_API bool glb_boundary_fixed;

void cwg::TetEdge::interface_weight(TetrahedralMesh* tetmesh, std::set<int>& interV)
{
	

	TetVertex* v0 = &tetmesh->m_storage_verts[ this->Verts[0] ];
	TetVertex* v1 = &tetmesh->m_storage_verts[ this->Verts[1] ];

	if (v0->NLabels() > 1 && v1->NLabels() > 1)// || v0->get_full_border_code() != Cleaver::VERTEX_INV)
	{
		if (v0->border_type() < v1->border_type())
		{
			std::swap(v0, v1);
			std::swap(Verts[0], Verts[1]);
		}

		cwg::SymMat4 Q4x4 = v0->QMatrix() + v1->QMatrix();

		//double diagval = Q4x4.diagval;
		double diagval = Q4x4(0,0);
		double fval = Q4x4(3,3);
		//Cleaver::vec3 q3x1(Q4x4.x, Q4x4.y, Q4x4.z);
		cv::Vec3d q3x1;
		for (int i=0; i<3; i++)
			q3x1(i) = Q4x4(i,3);

		Cleaver::vec3 q3x1_2;
		q3x1_2.x = q3x1[0];
		q3x1_2.y = q3x1[1];
		q3x1_2.z = q3x1[2];

		cv::Matx33d Q3x3;
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				Q3x3(i,j) = Q4x4(i,j);

	

		if (v0->NLabels() > 1 && v1->NLabels() > 1)
		{

			qem_pos(Q3x3, q3x1, v0, v1);

			cv::Vec3d conPos;
			conPos[0] = ContractPos.x;
			conPos[1] = ContractPos.y;
			conPos[2] = ContractPos.z;
			

			cv::Scalar rst = conPos.t()*Q3x3*conPos + conPos.t()*(q3x1)*2.0;
			this->Cost = rst[0] + fval;

			tetmesh->m_storage_verts[ this->Verts[0] ].setAvgWeight(MAX(1,this->Cost*100));
			tetmesh->m_storage_verts[ this->Verts[1] ].setAvgWeight(MAX(1,this->Cost*100));

			interV.insert(this->Verts[0]);
			interV.insert(this->Verts[1]);

			tetmesh->m_storage_verts[ this->Verts[0] ].weightLevel = 1;
			tetmesh->m_storage_verts[ this->Verts[1] ].weightLevel = 1;
													  
			tetmesh->m_storage_verts[ this->Verts[0] ].weighted = true;
			tetmesh->m_storage_verts[ this->Verts[1] ].weighted = true;

		}

			
	}
	
		
	
}

//
void cwg::TetEdge::set_debug_bool(bool to_debug_bool){
	debug_bool = to_debug_bool;
}

void cwg::TetEdge::ComputeContraction( TetrahedralMesh* tetmesh )
{

	TetVertex* v0 = &tetmesh->m_storage_verts[ this->Verts[0] ];
	TetVertex* v1 = &tetmesh->m_storage_verts[ this->Verts[1] ];

	if (v0->border_type() < v1->border_type())
	{
		std::swap(v0, v1);
		std::swap(Verts[0], Verts[1]);
	}

	cwg::SymMat4 Q4x4 ;

	v0->ComputeInitialQ();
	v1->ComputeInitialQ();

	Q4x4 = v0->QMatrix() + v1->QMatrix();

	double diagval = Q4x4(0,0);
	double fval = Q4x4(3,3);
	cv::Vec3d q3x1;
	for (int i=0; i<3; i++)
		q3x1(i) = Q4x4(i,3);

	Cleaver::vec3 q3x1_2;
	q3x1_2.x = q3x1[0];
	q3x1_2.y = q3x1[1];
	q3x1_2.z = q3x1[2];

	cv::Matx33d Q3x3;
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			Q3x3(i,j) = Q4x4(i,j);

	able_to_collapse = true;
	bool special_case = false;
	{
		if (v0->get_full_border_code() == Cleaver::VERTEX_INV)
		{
			if (v1->get_full_border_code() == Cleaver::VERTEX_INV)
				
			{
				// Case 1: Both order equals 1
				if (v0->NLabels() == 1 && v1->NLabels() == 1)
				{			
					qem_pos(diagval, q3x1_2);
				}//v0 order equals 1, and v1 >1
				else if (v0->NLabels() == 1 && v1->NLabels() > 1)
				{
					ContractPos = v1->pos();
				}//v1 order equals 1, and v0 >1
				else if (v0->NLabels() > 1 && v1->NLabels() == 1)
				{
					ContractPos = v0->pos();
				}
				else
				{
					special_case = true;
					if(v0->NLabels() > v1->NLabels())
						ContractPos = v0->pos();
					else if(v0->NLabels() < v1->NLabels())
						ContractPos = v1->pos();
					else
						qem_pos(Q3x3, q3x1, v0, v1);
				}
			}
			else // IMPORSSIBLE
				throw exception("void cwg::TetEdge::ComputeContraction( TetrahedronMesh* tetmesh ):\
								v0 < v1!");
		}
		else // v0 is external boundary
		{   
			if (v1->get_full_border_code() == Cleaver::VERTEX_INV)			
			{   //v1 order equals 1
				if (v1->NLabels() == 1)
				{
					ContractPos = v0->pos();			
				} // v1 order > 1
				else{
					if(v0->NLabels() > 1){
						ContractPos = v0->pos();
						special_case = true;
					}else{
						able_to_collapse = false;
					}
					
				}
			}
			else
			{
				able_to_collapse = !glb_boundary_fixed;  
				if (able_to_collapse)
				{
					unsigned char union_code, intersection_code;
					union_code = v0->get_border_code() | v1->get_border_code();
					intersection_code = v0->get_border_code() & v1->get_border_code();
					if (intersection_code)
					{
						//same face or edge
						if (v0->get_border_code() == union_code && v1->get_border_code() == union_code)
						{

							if (v0->NLabels() == 1 && v1->NLabels() == 1)
							{
								ContractPos = (v0->pos() + v1->pos())*0.5;
							}
							else if (v0->NLabels() == 1 && v1->NLabels() > 1)
							{
								ContractPos = v1->pos();		
							}
							else if (v0->NLabels() > 1 && v1->NLabels() == 1)
							{
								ContractPos = v0->pos();			
							}
							else
							{
								special_case = true;
								if(v0->NLabels() > v1->NLabels())
									ContractPos = v0->pos();
								else if(v0->NLabels() < v1->NLabels())
									ContractPos = v1->pos();
								else
									ContractPos = (v0->pos() + v1->pos())*0.5;
							}
						}
						else if (v0->get_border_code() == union_code && v1->get_border_code() != union_code)
						{
							if (v1->NLabels() == 1)
							{
								this->ContractPos = v0->pos();						
							}
							else
								able_to_collapse = false;
						}
						else if (v1->get_border_code() == union_code && v0->get_border_code() != union_code)
							throw exception("void cwg::TetEdge::ComputeContraction( TetrahedronMesh* tetmesh ):\
											unioncode: v0 < v1!");
						else
						{
							able_to_collapse = false;
						}
					}
					else
					{
						able_to_collapse = false;
					}
				}
			}
		}
	}

	if (able_to_collapse)
	{
		cv::Vec3d conPos;
		conPos[0] = ContractPos.x;
		conPos[1] = ContractPos.y;
		conPos[2] = ContractPos.z;
		
		cv::Scalar rst = conPos.t()*Q3x3*conPos + conPos.t()*(q3x1)*2.0;
		this->Cost = rst[0] + fval;

		this->compute_around_info(tetmesh);

		if(special_case){
			//check case 3.a
			int edge_order = get_order(tetmesh);
			if ((edge_order < v0->NLabels()) && (edge_order < v1->NLabels())){
				this->Cost = 1000;
				able_to_collapse = false;
			}
			else{//check case 3.b
				bool dual_edges_check = check_dual_edges(tetmesh);
				if(!dual_edges_check){
					this->Cost = 1000;
					able_to_collapse = false;
				}
			}
		}

		//check case 4, the quality check
		if(able_to_collapse){
			bool folderover = this->validate_folderover(tetmesh);		
			if (!folderover){
				this->Cost = 1000;
				able_to_collapse = false;
			}
		}

	}
	else{
		this->Cost = 1000;
	}
		
}

void cwg::TetEdge::qem_pos( const cv::Matx33d& Q3x3, const cv::Vec3d& q3x1, TetVertex* v0, TetVertex* v1 )
{
	//Matx33d vpos;
	double detQ3x3 = cv::determinant(Q3x3);
	if (abs(detQ3x3) < 1e-5)
		this->ContractPos = ( v0->pos() + v1->pos() ) * 0.5;
	else
	{
		cv::Vec3d vpos;
		bool flag = cv::solve( Q3x3, -q3x1, vpos);
		if (!flag)
			this->ContractPos = ( v0->pos() + v1->pos() ) * 0.5;
		else
		{
			ContractPos.x = vpos(0);
			ContractPos.y = vpos(1);
			ContractPos.z = vpos(2);
		}
	}
}

int cwg::TetEdge::find_vert( int vid ) const
{
	for (int i=0; i<2; i++)
	{
		if (Verts[i] == vid)
			return i;
	}
	return -1;
}

void cwg::TetEdge::compute_around_info( const TetrahedralMesh* tetmesh )
{
	ms_edge_around_info.clear();

	std::vector<int>* ptr_edge_tets                     = & ms_edge_around_info.edge_tets                 ;
	std::vector<int>* ptr_v0_private_tets				= & ms_edge_around_info.v0_private_tets           ;
	std::vector<int>* ptr_v1_private_tets				= & ms_edge_around_info.v1_private_tets           ;
	std::vector<int>* ptr_v0_private_overlap_tets		= & ms_edge_around_info.v0_private_overlap_tets   ;
	std::vector<int>* ptr_v0_private_nonoverlap_tets	= & ms_edge_around_info.v0_private_nonoverlap_tets;
	std::vector<int>* ptr_v1_private_overlap_tets		= & ms_edge_around_info.v1_private_overlap_tets   ;
	std::vector<int>* ptr_v1_private_nonoverlap_tets	= & ms_edge_around_info.v1_private_nonoverlap_tets;

	std::set<int>* ptr_affected_verts = & ms_edge_around_info.affected_verts;

	const TetVertex* v0 = &tetmesh->m_storage_verts[ this->Verts[0] ];
	const TetVertex* v1 = &tetmesh->m_storage_verts[ this->Verts[1] ];

	if (v0->border_type() < v1->border_type())
		throw exception("void cwg::TetEdge::compute_around_info( const TetrahedralMesh* tetmesh ): v0 < v1");

	const vector<int>& v0tets = v0->Tets;
	const vector<int>& v1tets = v1->Tets;

	set_intersection(v0tets.begin(), v0tets.end(), v1tets.begin(), v1tets.end(), inserter(*ptr_edge_tets, ptr_edge_tets->end()));
	set_difference(v0tets.begin(), v0tets.end(), v1tets.begin(), v1tets.end(), inserter(*ptr_v0_private_tets, ptr_v0_private_tets->end()));
	set_difference(v1tets.begin(), v1tets.end(), v0tets.begin(), v0tets.end(), inserter(*ptr_v1_private_tets, ptr_v1_private_tets->end()));

	vector<int> titjverts;
	vector<cwg::ordered_triple> v0_faces(ptr_v0_private_tets->size());
	for (int i=0; i<ptr_v0_private_tets->size(); i++)
	{
		const Tetrahedron* v0ti = &tetmesh->m_storage_tets[ (*ptr_v0_private_tets)[i] ];
		v0_faces[i] = v0ti->get_opposite_face(v0->ID);
	}

	vector<cwg::ordered_triple> v1_faces(ptr_v1_private_tets->size());
	for (int i=0; i<ptr_v1_private_tets->size(); i++)
	{
		const Tetrahedron* v1tj = &tetmesh->m_storage_tets[ (*ptr_v1_private_tets)[i] ];
		v1_faces[i] = v1tj->get_opposite_face(v1->ID);
	}

	for (int i=0; i<v0_faces.size(); i++)
	{
		for (int j=0; j<v1_faces.size(); j++)
		{
			if (v0_faces[i] == v1_faces[j])
			{
				//I added the following line, the if sentence below.
				ptr_v0_private_overlap_tets->push_back((*ptr_v0_private_tets)[i]);
				
				int t1 = (*ptr_v1_private_tets)[j];
				vector<int>::iterator ittmp = std::lower_bound(ptr_v1_private_overlap_tets->begin(), 
																ptr_v1_private_overlap_tets->end(), t1);
				if ( ittmp != ptr_v1_private_overlap_tets->end() )
				{
					if (*ittmp == t1){
						int x0 = (this->Verts[0] % ((tetmesh->slide_X+1) * (tetmesh->slide_Y+1))) / (tetmesh->slide_Y+1);
						int y0 = (this->Verts[0] % ((tetmesh->slide_X+1) * (tetmesh->slide_Y+1))) % (tetmesh->slide_Y+1);
						int z0 = (this->Verts[0] / ((tetmesh->slide_X+1) * (tetmesh->slide_Y+1)));

						int x1 = (this->Verts[1] % ((tetmesh->slide_X+1) * (tetmesh->slide_Y+1))) / (tetmesh->slide_Y+1);
						int y1 = (this->Verts[1] % ((tetmesh->slide_X+1) * (tetmesh->slide_Y+1))) % (tetmesh->slide_Y+1);
						int z1 = (this->Verts[1] / ((tetmesh->slide_X+1) * (tetmesh->slide_Y+1)));

						std::cout << "vert0: " << x0 << "\t" << y0 << "\t" << z0 << std::endl; 
						std::cout << "vert1: " << x1 << "\t" << y1 << "\t" << z1 << std::endl;
						throw exception("validate_folderover: *ittmp == *it1");
					}
						
					else
						ptr_v1_private_overlap_tets->insert(ittmp, t1);
				}
				else
					ptr_v1_private_overlap_tets->push_back(t1);

				ptr_affected_verts->insert(v0_faces[i].First());
				ptr_affected_verts->insert(v0_faces[i].Second());
				ptr_affected_verts->insert(v0_faces[i].Third());

			}
		}
	}

	set_difference(ptr_v0_private_tets->begin(), ptr_v0_private_tets->end(), 
		ptr_v0_private_overlap_tets->begin(), ptr_v0_private_overlap_tets->end(), 
		inserter(*ptr_v0_private_nonoverlap_tets, ptr_v0_private_nonoverlap_tets->end()));

	set_difference(ptr_v1_private_tets->begin(), ptr_v1_private_tets->end(), 
		ptr_v1_private_overlap_tets->begin(), ptr_v1_private_overlap_tets->end(), 
		inserter(*ptr_v1_private_nonoverlap_tets, ptr_v1_private_nonoverlap_tets->end()));

	ptr_affected_verts->insert(v0->ID);
}


bool cwg::TetEdge::check_labels_debug(TetrahedralMesh* tetmesh, int label_id){

	bool has_label = false;
	if((std::find(tetmesh->m_storage_verts[Verts[0]].m_labels.begin(), 
				  tetmesh->m_storage_verts[Verts[0]].m_labels.end(), label_id) != tetmesh->m_storage_verts[Verts[0]].m_labels.end()) ||
			    (std::find(tetmesh->m_storage_verts[Verts[1]].m_labels.begin(), 
						  tetmesh->m_storage_verts[Verts[1]].m_labels.end(), label_id) != tetmesh->m_storage_verts[Verts[1]].m_labels.end())){					  
							  has_label = true;
			}

	return has_label;
}


bool cwg::TetEdge::check_tet_debug(TetrahedralMesh* tetmesh, int tet_id){

	bool is_related=false;
	if (tetmesh->m_storage_tets[tet_id].ID == -1)
		return false;

	for(int vert_idx=0; vert_idx<4; vert_idx++){
		if (tetmesh->m_storage_tets[tet_id].Verts[vert_idx] == this->Verts[0] || 
			tetmesh->m_storage_tets[tet_id].Verts[vert_idx] == this->Verts[1]){
				is_related = true;
				break;
		}
	}

	return is_related;
}



bool cwg::TetEdge::validate_folderover(TetrahedralMesh* tetmesh)
{
	// v1 ===> v0, v0's border_type must >= v1's border_type
	int v0 = Verts[0];
	int v1 = Verts[1];


	const TetEdgeAroundInfo& e_around_info = get_around_info();
	const vector<int>& v0_private_tets = e_around_info.v0_private_tets;
	const vector<int>& v1_private_tets = e_around_info.v1_private_tets;
	const vector<int>& overlap_tets_v0 = e_around_info.v0_private_overlap_tets;
	const vector<int>& overlap_tets_v1 = e_around_info.v1_private_overlap_tets;
	const vector<int>& v0_private_nonoverlap_tets = e_around_info.v0_private_nonoverlap_tets;
	const vector<int>& v1_private_nonoverlap_tets = e_around_info.v1_private_nonoverlap_tets;
	const vector<int>& edge_tets = e_around_info.edge_tets;

	bool pass = true;

	// 0. Reserve original positions:
	Cleaver::vec3 v0origpos = tetmesh->m_storage_verts[v0].pos();
	Cleaver::vec3 v1origpos = tetmesh->m_storage_verts[v1].pos();

	if (overlap_tets_v0.size() > 0 || overlap_tets_v1.size() > 0){
		pass = false;
	}

	if(pass){
		// 1. Judge TetVolumes:
		vector<double> tetvolumes(v0_private_nonoverlap_tets.size()+v1_private_nonoverlap_tets.size());
		int k = 0;
		for (vector<int>::const_iterator it = v0_private_nonoverlap_tets.begin(); it != v0_private_nonoverlap_tets.end(); it++)
			tetvolumes[k++] = tetmesh->m_storage_tets[ *it ].compute_volume(tetmesh);			


		for (vector<int>::const_iterator it = v1_private_nonoverlap_tets.begin(); it != v1_private_nonoverlap_tets.end(); it++)
			tetvolumes[k++] = tetmesh->m_storage_tets[ *it ].compute_volume(tetmesh);

		// ***********************************************************************
		tetmesh->m_storage_verts[v0].pos() = tetmesh->m_storage_verts[v1].pos() = ContractPos;
		// ***********************************************************************

		k = 0;
		const double precision_vol = 0.1;
		const double precision_edgevariance = tetmesh->edge_var_thr;
		const double precision_mda = 0.0;//CV_PI / 180.0 * 180.0;

		if(pass)
		{
			for (vector<int>::const_iterator it = v0_private_nonoverlap_tets.begin(); pass && it != v0_private_nonoverlap_tets.end(); it++)
			{
				double tetvolume_after = tetmesh->m_storage_tets[ *it ].compute_volume(tetmesh);
				int flag_after = isgn(tetvolume_after, precision_vol);
				int flag_before = isgn(tetvolumes[k], precision_vol);
				if (flag_after == 0 || flag_after != flag_before){
					pass = false;
					break;
				}
				k++;
			}
		}

		if(pass)
		{
			for (vector<int>::const_iterator it = v1_private_nonoverlap_tets.begin(); pass && it != v1_private_nonoverlap_tets.end(); it++)
			{
				double tetvolume_after = tetmesh->m_storage_tets[ *it ].compute_volume(tetmesh);
				int flag_after = isgn(tetvolume_after, precision_vol);
				int flag_before = isgn(tetvolumes[k], precision_vol);
				if (flag_after == 0 || flag_after != flag_before){
					pass = false;
					break;
				}
				k++;
			}
		}

		//=============== compute the edge variation of each tet =======================
		if(pass)
		{
			for (vector<int>::const_iterator it = v0_private_nonoverlap_tets.begin(); pass && it != v0_private_nonoverlap_tets.end(); it++)
			{
				double edge_avg=1;
				double tetvariance = tetmesh->m_storage_tets[ *it ].variance(tetmesh, edge_avg);

				if( tetvariance > precision_edgevariance * edge_avg * edge_avg)
				{
					pass = false;
					break;
				}		

			}
		}

		if(pass){
			for (vector<int>::const_iterator it = v1_private_nonoverlap_tets.begin(); pass && it != v1_private_nonoverlap_tets.end(); it++){
				double edge_avg=1;
				double tetvariance = tetmesh->m_storage_tets[ *it ].variance(tetmesh,edge_avg);

				if( tetvariance > precision_edgevariance * edge_avg * edge_avg){
					pass = false;
					break;
				}		
			}
		}
	}

	tetmesh->m_storage_verts[v0].pos() = v0origpos;
	tetmesh->m_storage_verts[v1].pos() = v1origpos;

	return pass;
}

int cwg::TetEdge::get_order( const TetrahedralMesh* tetmesh ){
	std::vector<int>* ptr_edge_tets = & ms_edge_around_info.edge_tets;

	std::vector<int> edge_labels;

	for (int i=0; i<ptr_edge_tets->size(); i++)
	{
		const Tetrahedron* v0ti = &tetmesh->m_storage_tets[ (*ptr_edge_tets)[i] ];
		int lb = v0ti->get_label();
		std::vector<int>::iterator it = std::lower_bound(edge_labels.begin(), edge_labels.end(), lb);

		if (it != edge_labels.end()){
			if (*it != lb)
				edge_labels.insert(it, lb);
		}else
			edge_labels.push_back(lb);
	}

	return edge_labels.size();
}

int cwg::TetEdge::get_edge_order( const TetrahedralMesh* tetmesh, int v0, int v1){
	const TetVertex* vertex0 = &tetmesh->m_storage_verts[v0];
	const TetVertex* vertex1 = &tetmesh->m_storage_verts[v1];

	std::vector<int>* ptr_edge_tets;

	const vector<int>& v0tets = vertex0->Tets;
	const vector<int>& v1tets = vertex1->Tets;

	set_intersection(v0tets.begin(), v0tets.end(), v1tets.begin(), v1tets.end(), inserter(*ptr_edge_tets, ptr_edge_tets->end()));

	std::vector<int> edge_labels;

	for (int i=0; i<ptr_edge_tets->size(); i++)
	{
		const Tetrahedron* v0ti = &tetmesh->m_storage_tets[ (*ptr_edge_tets)[i] ];
		int lb = v0ti->get_label();
		std::vector<int>::iterator it = std::lower_bound(edge_labels.begin(), edge_labels.end(), lb);

		if (it != edge_labels.end()){
			if (*it != lb)
				edge_labels.insert(it, lb);
		}else
			edge_labels.push_back(lb);
	}

	return edge_labels.size();
}


bool cwg::TetEdge::check_dual_edges(const TetrahedralMesh* tetmesh){
	std::vector<int>* ptr_edge_tets = & ms_edge_around_info.edge_tets;

	bool dual_edges_pass = true;

	const TetVertex* v0 = &tetmesh->m_storage_verts[ this->Verts[0] ];
	const TetVertex* v1 = &tetmesh->m_storage_verts[ this->Verts[1] ];

	bool border_edge = (v0->get_full_border_code() != Cleaver::VERTEX_INV) && (v1->get_full_border_code() != Cleaver::VERTEX_INV);

	for (int i=0; i<ptr_edge_tets->size(); i++)
	{
		const Tetrahedron* ti = &tetmesh->m_storage_tets[ (*ptr_edge_tets)[i] ];
		bool isolated = true;

		int other_verts[2];
		int other_vert_index=0;

		for(int vert_i= 0; vert_i < 4; vert_i++){
			if( (ti->Verts[vert_i] != Verts[0]) && (ti->Verts[vert_i] != Verts[1]) ){
				other_verts[other_vert_index] = ti->Verts[vert_i];
				other_vert_index += 1;
			}

			if (other_vert_index >= 2)
				break;
		}

		const TetVertex* v2 = &tetmesh->m_storage_verts[ other_verts[0] ];
		const TetVertex* v3 = &tetmesh->m_storage_verts[ other_verts[1] ];

		const vector<int>& v0tets = v0->Tets;
		const vector<int>& v1tets = v1->Tets;
		const vector<int>& v2tets = v2->Tets;
		const vector<int>& v3tets = v3->Tets;
		
		std::vector<int> face_tets_v0_v2_v3;
		std::vector<int> tmp_face_tets;
		//get the tets of v0-v2-v3, where v2 and v3 are the two other verts
		set_intersection(v2tets.begin(), v2tets.end(), v3tets.begin(), v3tets.end(), back_inserter(tmp_face_tets));
		set_intersection(tmp_face_tets.begin(), tmp_face_tets.end(), v0tets.begin(), v0tets.end(), back_inserter(face_tets_v0_v2_v3));
		//get the tets of v1-v2-v3, where v2 and v3 are the verts of this edge.
		std::vector<int> face_tets_v1_v2_v3;
		set_intersection(tmp_face_tets.begin(), tmp_face_tets.end(), v1tets.begin(), v1tets.end(), back_inserter(face_tets_v1_v2_v3));

		//check if face v0-v2-v3 has dual edge, if face v1-v2-v3 has dual edge
		bool has_dual_edge_023=false, has_volume_dual_edge_023=false; 
		bool has_dual_edge_123=false, has_volume_dual_edge_123=false;

		if(face_tets_v0_v2_v3.size() == 2)
			 has_dual_edge_023 = (tetmesh->m_storage_tets[face_tets_v0_v2_v3[0]].get_label() == tetmesh->m_storage_tets[face_tets_v0_v2_v3[1]].get_label());
		if((v0->get_full_border_code() != Cleaver::VERTEX_INV) && 
			(v2->get_full_border_code() != Cleaver::VERTEX_INV) && 
			(v3->get_full_border_code() != Cleaver::VERTEX_INV))
		    has_volume_dual_edge_023 = true;

		if(face_tets_v1_v2_v3.size() == 2)
			 has_dual_edge_123 = (tetmesh->m_storage_tets[face_tets_v1_v2_v3[0]].get_label() == tetmesh->m_storage_tets[face_tets_v1_v2_v3[1]].get_label());
		if((v1->get_full_border_code() != Cleaver::VERTEX_INV) && 
			(v2->get_full_border_code() != Cleaver::VERTEX_INV) && 
			(v3->get_full_border_code() != Cleaver::VERTEX_INV))
			has_volume_dual_edge_123 = true;

		tmp_face_tets.clear();
		//get the tets of v0-v1-v2, where v0 and v1 are the verts of this edge.
		std::vector<int> face_tets_v0_v1_v2;
		set_intersection(v0tets.begin(), v0tets.end(), v1tets.begin(), v1tets.end(), back_inserter(tmp_face_tets));
		set_intersection(tmp_face_tets.begin(), tmp_face_tets.end(), v2tets.begin(), v2tets.end(), back_inserter(face_tets_v0_v1_v2));
		//get the tets of v0-v1-v3, where v0 and v1 are the verts of this edge.
		std::vector<int> face_tets_v0_v1_v3;
		set_intersection(tmp_face_tets.begin(), tmp_face_tets.end(), v3tets.begin(), v3tets.end(), back_inserter(face_tets_v0_v1_v3));

		bool has_dual_edge_012=false, has_volume_dual_edge_012=false;
		//check material-dual-edge of v0-v1-v2
		if(face_tets_v0_v1_v2.size() == 2)
			has_dual_edge_012 = tetmesh->m_storage_tets[face_tets_v0_v1_v2[0]].get_label() == tetmesh->m_storage_tets[face_tets_v0_v1_v2[1]].get_label();
		//check volume-dual-edge of v0-v1-v2
		if((v0->get_full_border_code() != Cleaver::VERTEX_INV) &&
			(v1->get_full_border_code() != Cleaver::VERTEX_INV) && 
			(v2->get_full_border_code() != Cleaver::VERTEX_INV))
			 has_volume_dual_edge_012 = true;

		bool has_dual_edge_013=false, has_volume_dual_edge_013=false;
		//check material-dual-edge of v0-v1-v3
		if(face_tets_v0_v1_v3.size() == 2)
			has_dual_edge_013 = tetmesh->m_storage_tets[face_tets_v0_v1_v3[0]].get_label() == tetmesh->m_storage_tets[face_tets_v0_v1_v3[1]].get_label();
		//check volume-dual-edge of v0-v1-v3
		if((v0->get_full_border_code() != Cleaver::VERTEX_INV) &&
			(v1->get_full_border_code() != Cleaver::VERTEX_INV) &&
			(v3->get_full_border_code() != Cleaver::VERTEX_INV))
			has_volume_dual_edge_013 = true;

		bool has_dual_face_02_12 = false, has_dual_face_03_13=false;
		bool has_volume_dual_face_02_12 = false, has_volume_dual_face_03_13=false;

		//check material dual face of edge 02 and edge 12
		has_dual_face_02_12 = check_dual_faces(tetmesh, this->Verts[0], other_verts[0]) || check_dual_faces(tetmesh, this->Verts[1], other_verts[0]);
		//check material dual face of edge 03 and edge 13
		has_dual_face_03_13 = check_dual_faces(tetmesh, this->Verts[0], other_verts[1]) || check_dual_faces(tetmesh, this->Verts[1], other_verts[1]);
		//check volume dual face of edge 02 and edge 12
		has_volume_dual_face_02_12 = check_volume_dual_faces(tetmesh, this->Verts[0], other_verts[0]) || 
												check_volume_dual_faces(tetmesh, this->Verts[1], other_verts[0]);
		//check volume dual face of edge 03 and edge 13
		has_volume_dual_face_03_13 = check_volume_dual_faces(tetmesh, this->Verts[0], other_verts[1]) || 
												check_volume_dual_faces(tetmesh, this->Verts[1], other_verts[1]);

		// material connection filtering
		if( has_dual_edge_023 || has_dual_edge_123){
			isolated = false;

			if(has_dual_edge_012)
				if(!has_dual_face_02_12)
					return false;

			if(has_dual_edge_013)
				if(!has_dual_face_03_13)
					return false;

		}

		//volume connection filtering
		if(has_volume_dual_edge_012)
			if(!has_volume_dual_face_02_12)
				return false;
			
		if(has_volume_dual_edge_013)
			if(!has_volume_dual_face_03_13)
				return false;


		if(isolated)
			return !isolated;
	}
	
	
	return true;
}

bool cwg::TetEdge::check_dual_faces( const TetrahedralMesh* tetmesh, int v0, int v1){
		
		const TetVertex* vertex0 = &tetmesh->m_storage_verts[ v0 ];
		const TetVertex* vertex1 = &tetmesh->m_storage_verts[ v1 ];

		if (!((vertex0->get_full_border_code() == Cleaver::VERTEX_INV) || (vertex1->get_full_border_code() == Cleaver::VERTEX_INV) ))
			return false;

		const vector<int>& v0tets = vertex0->Tets;
		const vector<int>& v1tets = vertex1->Tets;
		
		//get the tets of v0-v1 to see if it is a inner edge(has dual face)
		std::vector<int> v0_v1_tets;
		set_intersection(v0tets.begin(), v0tets.end(), v1tets.begin(), v1tets.end(), back_inserter(v0_v1_tets));

		int label_check = tetmesh->m_storage_tets[v0_v1_tets[0]].get_label();
		bool has_dual_face = true;
		for(int tet_index=1; tet_index < v0_v1_tets.size(); tet_index++){
			if( tetmesh->m_storage_tets[v0_v1_tets[tet_index]].get_label() != label_check){
				has_dual_face = false;
				break;
			}
		}

		return has_dual_face;
		
}


bool cwg::TetEdge::check_volume_dual_faces( const TetrahedralMesh* tetmesh, int v0, int v1){
		
		const TetVertex* vertex0 = &tetmesh->m_storage_verts[ v0 ];
		const TetVertex* vertex1 = &tetmesh->m_storage_verts[ v1 ];

		const vector<int>& v0tets = vertex0->Tets;
		const vector<int>& v1tets = vertex1->Tets;
		
		//get the tets of v0-v1 to see if it is a inner edge(has dual face)
		std::vector<int> v0_v1_tets;
		set_intersection(v0tets.begin(), v0tets.end(), v1tets.begin(), v1tets.end(), back_inserter(v0_v1_tets));

		int label_check = tetmesh->m_storage_tets[v0_v1_tets[0]].get_label();
		bool has_dual_face = true;
		for(int tet_index=1; tet_index < v0_v1_tets.size(); tet_index++){
			if( tetmesh->m_storage_tets[v0_v1_tets[tet_index]].get_label() != label_check){
				has_dual_face = false;
				break;
			}
		}

		return has_dual_face;
		
}

void cwg::TetEdge::debug_obj(const TetrahedralMesh* tetmesh, string filename, int tet_index){

	ofstream output_file(filename);
	
	for (int vert_i =0; vert_i<4; vert_i++){
		output_file << "v " << tetmesh->m_storage_verts[tetmesh->m_storage_tets[tet_index].Verts[vert_i]].pos().x << " " <<
			                  tetmesh->m_storage_verts[tetmesh->m_storage_tets[tet_index].Verts[vert_i]].pos().y << " " <<
			                  tetmesh->m_storage_verts[tetmesh->m_storage_tets[tet_index].Verts[vert_i]].pos().z << " " << std::endl;

	}

	for (int vert_i = 0; vert_i < 4; vert_i++){
		for(int vert_j = 0; vert_j < vert_i; vert_j++){
			output_file << "l " << vert_i+1 << " " << vert_j+1 << std::endl;
		}
	}

	for(int face_idx=0; face_idx<4; face_idx++){
		output_file << "f " << (face_idx) % 4 + 1 << " " << (face_idx+1)%4 +1 << " " << (face_idx+2)%4 + 1<< std::endl; 
	}

	output_file.close();
}

void cwg::TetEdge::debug_edge(const TetrahedralMesh* tetmesh, string filename, int vert_0, int vert_1){

	ofstream output_file(filename);
	

	output_file << "v " << tetmesh->m_storage_verts[vert_0].pos().x << " " <<
			               tetmesh->m_storage_verts[vert_0].pos().y << " " <<
			               tetmesh->m_storage_verts[vert_0].pos().z << " " << std::endl;

	output_file << "v " << tetmesh->m_storage_verts[vert_1].pos().x << " " <<
			               tetmesh->m_storage_verts[vert_1].pos().y << " " <<
			               tetmesh->m_storage_verts[vert_1].pos().z << " " << std::endl;

	output_file << "l " << 1 << " " << 2 << std::endl;

}

bool cwg::TetEdge::check_edge_pos(const TetrahedralMesh* tetmesh, double x_0, double y_0, double z_0, double x_1, double y_1, double z_1){

	double diff0 = abs(tetmesh->m_storage_verts[Verts[0]].pos().x - x_0) + abs(tetmesh->m_storage_verts[Verts[0]].pos().y - y_0) +abs(tetmesh->m_storage_verts[Verts[0]].pos().z - z_0);
	double diff1 = abs(tetmesh->m_storage_verts[Verts[1]].pos().x - x_1) + abs(tetmesh->m_storage_verts[Verts[1]].pos().y - y_1) +abs(tetmesh->m_storage_verts[Verts[1]].pos().z - z_1);
	double diff3 = abs(tetmesh->m_storage_verts[Verts[1]].pos().x - x_0) + abs(tetmesh->m_storage_verts[Verts[1]].pos().y - y_0) +abs(tetmesh->m_storage_verts[Verts[1]].pos().z - z_0);
	double diff2 = abs(tetmesh->m_storage_verts[Verts[0]].pos().x - x_1) + abs(tetmesh->m_storage_verts[Verts[0]].pos().y - y_1) +abs(tetmesh->m_storage_verts[Verts[0]].pos().z - z_1);

	if((diff0<0.1 && diff1<0.1) || (diff2<0.1 && diff3<0.1))
		return true;
	else
		return false;

}