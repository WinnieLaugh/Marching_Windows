#ifndef _TET_EDGE_HEAP_H_
#define _TET_EDGE_HEAP_H_

#include "TetEdge.h"

namespace cwg
{
	template<typename TCOMP>
	struct TetEdgeComp
	{
		TCOMP comp;
		bool operator() (const TetEdge* left, const TetEdge* right)
		{
			return comp(left->cost(), right->cost());
		}
	};

	class TetEdgeHeap
	{
	public:
		TetEdgeHeap();
		~TetEdgeHeap();
		
		void									insert(TetEdge* e);
		void									update(int eid);//
		void									pop();//
		void									remove(int eid);//
		void									clear();

		inline TetEdge*							top_edge() const { return (m_size > 0 ? m_heap[1] : NULL); }
		inline int								size() const { return m_size; }
		inline bool								empty() const { return m_size==0; }

		inline void								resize(int num) { m_heap.resize(num+1, NULL); /*m_hash_pos.resize(num, 0);*/ m_size = num; }
		inline void								set_edge( TetEdge* e, int i ) { m_heap[i+1] = e; /*m_hash_pos[e->id()] = i+1; */ m_hash_pos.push_back(i+1); }
		inline TetEdge*							get_edge( int i ) const { return m_heap[i+1]; }
		void									heap_sort();
		void									random_shuffle(); 

		inline TetEdge*							operator[] (int eid) const { return m_heap[m_hash_pos[eid]]; } 
		//This is only for test! m_hash_pos should be private
		std::vector<int>						m_hash_pos;

	private:
		int										bubble_up(int pos);
		int										bubble_down(int pos);

		std::vector<TetEdge*>					m_heap;
		

		int										m_size;

		TetEdgeComp< std::less<double> >		m_comp_less;
		TetEdgeComp< std::greater<double> >		m_comp_greater;
	};
}

#endif