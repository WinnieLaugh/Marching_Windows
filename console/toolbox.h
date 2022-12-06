#ifndef _CWG_CONSOLE_TOOLBOX_H_
#define _CWG_CONSOLE_TOOLBOX_H_

#include <cwgUtilities.h>
#include <simplification3d.h>

#include <Visualize.h>

//************************************
// Method:    convert2svdata
// FullName:  convert2svdata
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: int argc
// Parameter: char * * argv
// INPUT: Yating's labelled binary files
// OUTPUT: svdata
//************************************
void				convert2svdata(int argc, char** argv);

//************************************
// Method:    tetmesh_simplify_batch
// FullName:  tetmesh_simplify_batch
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: int argc
// Parameter: char * * argv
// INPUT: svdata
// OUTPUT: simplified plys both in volume and surface
//************************************
void				tetmesh_simplify_batch(int argc, char** argv);

//************************************
// Method:    tetmesh_simplify
// FullName:  tetmesh_simplify
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: int argc
// Parameter: char * * argv
//************************************
void				tetmesh_simplify(int argc, char** argv);

//************************************
// Method:    split_whole_data_into_blocks
// FullName:  split_whole_data_into_blocks
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: int argc
// Parameter: char * * argv
//************************************
void				split_whole_data_into_blocks(int argc, char** argv);
void				save_tsp_blocks(const cwg::TSPVideo& video, 
						int ks, int ke, int is, int ie, int js, int je, 
						const std::string& output_foldername);

//************************************
// Method:    merge_tetm
// FullName:  merge_tetm
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: int argc
// Parameter: char * * argv
//************************************
void				merge_tetm(int argc, char** argv);

void				resize_tetm(int argc, char** argv);


void				improve_tetm(int argc, char** argv);
//*****************************************
#endif