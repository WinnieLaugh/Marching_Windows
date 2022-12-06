#ifndef _CWG_CLEAVER_VERTEX_INDEX_DEF_H_
#define _CWG_CLEAVER_VERTEX_INDEX_DEF_H_

#ifndef _VERTEX_INDEX_
#define _VERTEX_INDEX_

enum cleaver_vertex_index
{
	_O   = -1,// ___      NO MORE VERTICES
	_A   	 ,//
	_B   	 ,//     \__  Lattice Vertices
	_C   	 ,//     /
	_D   	 ,// ___/

	_AB  	 ,//
	_AC  	 ,//
	_AD  	 ,//      \__ Cutpoint Vertices
	_BC  	 ,//      /
	_CD  	 ,//
	_BD  	 ,// ___/

	_ABC 	 ,//
	_ACD 	 ,//     \__  TriplePoint Vertices
	_ABD 	 ,//     /
	_BCD 	 ,// ___/

	_ABCD	  // QuadPoint Vertex
};

#define VERT 0
#define CUT  1
#define TRIP 2
#define QUAD 3

#endif

#endif