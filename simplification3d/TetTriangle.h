#ifndef _TRIANGLE_H_
#define _TRIANGLE_H_

#include "TetVertex.h"

namespace cwg
{
	class TetTriangle
	{
	public:
		TetTriangle() : ID(ms_next_id++) { m_triple = NULL; Tets[0] = Tets[1] = -1; }
		TetTriangle(int v0, int v1, int v2) : ID(ms_next_id++) { SetVertices(v0, v1, v2); Tets[0] = Tets[1] = -1; }
		~TetTriangle() { SAFE_DELETE(m_triple); }

		inline void							SetVertices(int v0, int v1, int v2) { m_triple = new cwg::ordered_triple(v0,v1,v2); }
		inline int							Vert(int i) const { return (*m_triple)(i); }
		inline static void					reset_nextid() { ms_next_id = 0; }
		
		int									Tets[2];
		int									ID;
		//bool								interSurf;
	private:
		cwg::ordered_triple*				m_triple;
		static int							ms_next_id;
	};
}

#endif