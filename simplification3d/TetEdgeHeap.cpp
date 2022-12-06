#include "stdafx.h"
#include "TetEdgeHeap.h"

using namespace std;
using namespace cwg;

cwg::TetEdgeHeap::TetEdgeHeap() : m_size(0)
{
	
	m_heap.push_back(NULL);
}

cwg::TetEdgeHeap::~TetEdgeHeap()
{
	m_heap.clear();
	m_size = 0;
	m_hash_pos.clear();
}

void cwg::TetEdgeHeap::clear()
{
	m_heap.clear();
	m_size = 0;
	m_hash_pos.clear();
	m_heap.push_back(NULL);
}

void cwg::TetEdgeHeap::insert( TetEdge* e )
{
	
	//assert(0);
	m_heap.push_back(e);
	if (e->id() > m_hash_pos.size()-1)
		m_hash_pos.resize(e->id()+1);
	m_size++;
	m_hash_pos[e->id()] = m_size;
	m_hash_pos[e->id()] = bubble_up(m_size);
}

int cwg::TetEdgeHeap::bubble_up( int pos )
{
	TetEdge* e = m_heap[pos];
	size_t father_pos = pos >> 1;
	
	for (; pos > 1; pos = father_pos, father_pos >>= 1)
	{
		if ( !m_comp_less(e, m_heap[father_pos]) )
			break;

		m_heap[pos] = m_heap[father_pos];
		m_hash_pos[ m_heap[pos]->id() ] = pos;
	}

	m_heap[pos] = e;
	m_hash_pos[e->id()] = pos;
	return pos;
}

int cwg::TetEdgeHeap::bubble_down( int pos )
{
	TetEdge* e = m_heap[pos];
	size_t child_pos = pos << 1;
	for (; child_pos <= m_size; pos = child_pos, child_pos <<= 1)
	{
		if ( (child_pos != m_size) && m_comp_less(m_heap[child_pos+1], m_heap[child_pos]) )
			child_pos++;

		if ( !m_comp_less(m_heap[child_pos], e) )
			break;

		m_heap[pos] = m_heap[child_pos];
		m_hash_pos[m_heap[pos]->id()] = pos;
	}
	m_heap[pos] = e;
	m_hash_pos[e->id()] = pos;
	return pos;
}

void cwg::TetEdgeHeap::pop()
{
	if (m_size == 0)
		throw exception("void cwg::TetEdgeHeap::pop(): heap size == 0 now, cannot pop!");

	TetEdge* e = m_heap[1];
	m_hash_pos[e->id()] = 0;
	m_heap[1] = m_heap[m_size--];
	m_heap.pop_back();
	if (m_size == 0) return;
	bubble_down(1);
}

void cwg::TetEdgeHeap::heap_sort()
{
	make_heap(m_heap.begin()+1, m_heap.end(), m_comp_greater);
	for (int i=1; i<=m_size; i++)
		m_hash_pos[ m_heap[i]->id() ] = i;
}

void cwg::TetEdgeHeap::update( int eid )
{
	size_t pos = m_hash_pos[eid];
	assert(m_heap[pos]->id() == eid);
	if (pos == 0)
		throw exception("void cwg::TetEdgeHeap::update( int eid ): pos == 0");
	bubble_down(pos);
	bubble_up(pos);
}

void cwg::TetEdgeHeap::random_shuffle()
{
	
	//srand(time(0));
	srand(1);
	std::random_shuffle(m_heap.begin()+1, m_heap.end());
	heap_sort();
}

void cwg::TetEdgeHeap::remove( int eid )
{
	size_t pos = m_hash_pos[eid];
	m_hash_pos[eid] = 0;
	m_heap[pos] = m_heap[m_size];
	m_size--;
	m_heap.pop_back();
	if (m_size == pos-1) 
		return;
	bubble_down(pos);
	bubble_up(pos);
}




