#include "stdafx.h"
#include "TetCube.h"

using namespace std;
using namespace cv;

void cwg::TetCube::ComputeTetrahedra(std::vector<Tetrahedron>::iterator& it)
{
	it->Verts[0] = Verts[0];
	it->Verts[1] = Verts[2];
	it->Verts[2] = Verts[3];
	it->Verts[3] = Verts[6];
	it->set_label(m_label);
	it++;

	it->Verts[0] = Verts[0];
	it->Verts[1] = Verts[3];
	it->Verts[2] = Verts[4];
	it->Verts[3] = Verts[6];
	it->set_label(m_label);
	it++;

	it->Verts[0] = Verts[0];
	it->Verts[1] = Verts[1];
	it->Verts[2] = Verts[3];
	it->Verts[3] = Verts[4];
	it->set_label(m_label);
	it++;

	it->Verts[0] = Verts[7];
	it->Verts[1] = Verts[1];
	it->Verts[2] = Verts[4];
	it->Verts[3] = Verts[5];
	it->set_label(m_label);
	it++;

	it->Verts[0] = Verts[7];
	it->Verts[1] = Verts[3];
	it->Verts[2] = Verts[4];
	it->Verts[3] = Verts[6];
	it->set_label(m_label);
	it++;

	it->Verts[0] = Verts[7];
	it->Verts[1] = Verts[1];
	it->Verts[2] = Verts[3];
	it->Verts[3] = Verts[4];
	it->set_label(m_label);
	it++;
}

void cwg::TetCube::ComputeTetrahedra(std::vector<Tetrahedron>::iterator& it, int &tID)
{
	it->Verts[0] = Verts[0];
	it->Verts[1] = Verts[2];
	it->Verts[2] = Verts[3];
	it->Verts[3] = Verts[6];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;

	it->Verts[0] = Verts[0];
	it->Verts[1] = Verts[3];
	it->Verts[2] = Verts[4];
	it->Verts[3] = Verts[6];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;

	it->Verts[0] = Verts[0];
	it->Verts[1] = Verts[1];
	it->Verts[2] = Verts[3];
	it->Verts[3] = Verts[4];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;

	it->Verts[0] = Verts[7];
	it->Verts[1] = Verts[1];
	it->Verts[2] = Verts[4];
	it->Verts[3] = Verts[5];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;

	it->Verts[0] = Verts[7];
	it->Verts[1] = Verts[3];
	it->Verts[2] = Verts[4];
	it->Verts[3] = Verts[6];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;

	it->Verts[0] = Verts[7];
	it->Verts[1] = Verts[1];
	it->Verts[2] = Verts[3];
	it->Verts[3] = Verts[4];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;
}

void cwg::TetCube::ComputeTetrahedra2(std::vector<Tetrahedron>::iterator& it, int &tID)
{
	it->Verts[0] = Verts[0];
	it->Verts[1] = Verts[1];
	it->Verts[2] = Verts[2];
	it->Verts[3] = Verts[4];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;

	it->Verts[0] = Verts[1];
	it->Verts[1] = Verts[4];
	it->Verts[2] = Verts[5];
	it->Verts[3] = Verts[7];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;

	it->Verts[0] = Verts[2];
	it->Verts[1] = Verts[4];
	it->Verts[2] = Verts[6];
	it->Verts[3] = Verts[7];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;

	it->Verts[0] = Verts[1];
	it->Verts[1] = Verts[2];
	it->Verts[2] = Verts[3];
	it->Verts[3] = Verts[7];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;

	it->Verts[0] = Verts[1];
	it->Verts[1] = Verts[2];
	it->Verts[2] = Verts[4];
	it->Verts[3] = Verts[7];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;

	
}

void cwg::TetCube::ComputeTetrahedra3(std::vector<Tetrahedron>::iterator& it, int &tID)
{
	it->Verts[0] = Verts[0];
	it->Verts[1] = Verts[1];
	it->Verts[2] = Verts[3];
	it->Verts[3] = Verts[5];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;

	it->Verts[0] = Verts[0];
	it->Verts[1] = Verts[2];
	it->Verts[2] = Verts[3];
	it->Verts[3] = Verts[6];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;

	it->Verts[0] = Verts[3];
	it->Verts[1] = Verts[5];
	it->Verts[2] = Verts[6];
	it->Verts[3] = Verts[7];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;

	it->Verts[0] = Verts[0];
	it->Verts[1] = Verts[4];
	it->Verts[2] = Verts[5];
	it->Verts[3] = Verts[6];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;

	it->Verts[0] = Verts[0];
	it->Verts[1] = Verts[3];
	it->Verts[2] = Verts[5];
	it->Verts[3] = Verts[6];
	it->set_label(m_label);
	it->ID = tID;
	++tID;
	it++;


}