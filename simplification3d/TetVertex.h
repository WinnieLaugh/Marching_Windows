#ifndef _TET_VERTEX_H_
#define _TET_VERTEX_H_

#include "Cleaverlib/vec3.h"
#include "Cleaverlib/Vertex.h"
#include "SymMat4.h"
#include <cwgUtilities.h>

namespace cwg
{
	class Tetrahedron;
	class TetEdge;
	class TetTriangle;
	class TetrahedralMesh;

	class QMat
	{
	public:
		QMat() { diagval = x = y = z = fval = 0.0; }
		QMat(double _diag, double _x, double _y, double _z, double _fval) : diagval(_diag),
			x(_x), y(_y), z(_z), fval(_fval) {}

		inline void		set2zero() { diagval = x = y = z = fval = 0.0; }

		double			diagval;
		double			x, y, z;
		double			fval;
		

		inline QMat operator+(const QMat& right) { return QMat(diagval+right.diagval, x+right.x, y+right.y, z+right.z, fval+right.fval); }
		inline QMat& operator+=(const QMat& right) { diagval+=right.diagval; x+=right.x; y+=right.y; z+=right.z; fval+=right.fval; return *this; }
		inline QMat operator-(QMat& right) { return QMat(diagval-right.diagval, x-right.x, y-right.y, z-right.z, fval-right.fval); }
		inline QMat operator*(double w) { return QMat(diagval*w, x*w, y*w, z*w, fval*w); }
		inline QMat& operator*=(double w) { diagval*=w; x*=w; y*=w; z*=w; fval*=w; return *this; }
		inline QMat operator/(double w) { return QMat(diagval/w, x/w, y/w, z/w, fval/w); }
	};

	class TetVertex
	{
	public:
		TetVertex() : m_xyz(Cleaver::vec3::zero), m_funcval(0.0), avg_weight(1.0), weightLevel(0), lfs(0), weighted(false), m_border_code(255), ID(ms_nextid++), slideBoundary(false), last_slide(false) {}
		TetVertex(Cleaver::Vertex3D* v) : m_xyz(v->pos()), m_funcval(v->pos().dot(v->pos())), 
			 avg_weight(1.0), weightLevel(0), lfs(0), weighted(false), m_border_code(255), ID(ms_nextid++), slideBoundary(false), last_slide(false) {}
		TetVertex(int id, Cleaver::vec3 xyz, double fval, double w, int bordercode) : ID(id), m_funcval(fval), 
			 avg_weight(1.0), weightLevel(0), lfs(0), weighted(false) , slideBoundary(false), last_slide(false){ m_border_code = bordercode; m_xyz = xyz; }
		~TetVertex() {}

		void								invalidate();
		void								clear();
		inline void							set_full_border_code(unsigned char t) { m_border_code = t; }
		inline unsigned char				get_full_border_code() const { return m_border_code; }	//8 bits, first 2 - border type, rest 6 - which facet 
		inline unsigned char				get_border_code() const { return (0x3F & m_border_code); }	//return the last 6 bits - which facet
		inline unsigned char				border_type() const { return (m_border_code>>6); }	//border type. 0-interior, 1-face, 2-edge, 3-corner
		
		void								assign_label(int lb);
		void								update_labels(TetrahedralMesh* tetmesh);
		int									get_label(int i) const {return m_labels[i];}
		void								clear_label(){m_labels.clear();}

		void								correlate_tet(int tet);
		void								decorrelate_tet(int tet);
		void								correlate_boundary_tri(int tri);
		inline SymMat4						QMatrix() const { return m_Q; }
		inline SymMat4&						QMatrix() { return m_Q; }
		//inline SymMat4						QMatrix_surf() const { return m_Q_surf; }
		//inline SymMat4&						QMatrix_surf() { return m_Q_surf; }

		inline const Cleaver::vec3&			pos() const { return m_xyz; }
		inline Cleaver::vec3&				pos() { return m_xyz; }
		inline Cleaver::vec3&				tempPos() { return temp_xyz; }
		inline Cleaver::vec3&				color() { return m_color; }
		void								setPos(Cleaver::vec3 &p) { m_xyz = p;}
		void								setTempPos(Cleaver::vec3 &p) { temp_xyz = p;}

		inline double						funcval() const { return m_funcval; }
		inline void							update_funcval() { m_funcval = m_xyz.dot(m_xyz); }
		double								weight() const{ return avg_weight; }
		//double							avgWeight(){ return avg_weight; }
		
		inline void							setAvgWeight(double w)  { avg_weight = w; }

		void								ComputeInitialQ();
		void								ComputeInitialQ( TetrahedralMesh* tetmesh );
		inline void							UpdateQ(const SymMat4& anotherQ) { m_Q += anotherQ; }
		inline void							SetQ(const SymMat4& QMat) {m_Q = QMat; }
		inline int							NLabels() const { return m_labels.size(); }

		void								read(std::ifstream& ifile);
		void								write(std::ofstream& ofile) const;

		void								copy_info_from_cleaver_vert(Cleaver::Vertex3D* v);
		inline static void					reset_nextid() { ms_nextid = 0; }

		
		friend std::istream&				operator>>(std::istream& is, TetVertex& obj);
		friend std::ostream&				operator<<(std::ostream& os, const TetVertex& obj);

		int									ID;
		std::vector<int>					Tets;
		std::vector<int>					OrigBoundaryTris;
		std::vector<int>					OuterBoundaryTris;	//triangles on outer boundarys for smoothing
		int									weightLevel;	//set weight
		bool								weighted;	//weighted for once
		
		//=========== sliding boundary - cannot simplify in this slide ==========
		bool								slideBoundary;
		//======================wenhua===========================================
		bool								boundary;
		int									keep;
		bool								last_slide;	//the v is in last slide, going to be saved

		//std::set<int>						parents;
		double								lfs;	//only used on surface vertices - local feature size



		//==========test===========================
		//SymMat4								m_Q_s;	//org qem
		//SymMat4								m_Q_d;	//new qem

		std::vector<int>					m_labels;
	private:
		//QMat								m_Q;
		SymMat4								m_Q;
		//SymMat4								m_Q_surf;

		

		double								m_funcval;
		Cleaver::vec3						m_xyz;
		Cleaver::vec3						temp_xyz;
		Cleaver::vec3						m_color;
		
		//std::vector<double>					m_weight;
		double								avg_weight;
		unsigned char						m_border_code;
		

	private:
		static int							ms_nextid;
	};

	struct TetVertexIDComp
	{
		bool								operator() (const TetVertex* left, const TetVertex* right) { return left->ID < right->ID; }
	};

	std::istream&							operator>>(std::istream& is, TetVertex& obj);
	std::ostream&							operator<<(std::ostream& os, const TetVertex& obj);

	cwg::SymMat4							compute_tri_Q( TetVertex* v0, TetVertex* v1, TetVertex* v2 );
}

#endif