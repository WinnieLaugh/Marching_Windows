// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the SIMPLIFICATION3D_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// SIMPLIFICATION3D_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifndef _SIMPLIFICATION_3D_H_
#define _SIMPLIFICATION_3D_H_

#ifdef SIMPLIFICATION3D_EXPORTS
#define SIMPLIFICATION3D_API __declspec(dllexport)
#else
#define SIMPLIFICATION3D_API __declspec(dllimport)
#endif

#include "TspVideo.h"
#include "TetrahedralMesh.h"

#include <Visualize.h>

namespace cwg
{
	class SIMPLIFICATION3D_API Simplification3D
	{
	public:
		Simplification3D();
		~Simplification3D();

		
		void									LoadVolumeData(const std::string& foldername);

		void									LoadMeshData(const std::string& filename);

		void									GenerateSingleLabelData(int h, int w, int nframes);
		bool									LoadDensity();

		void									CleaverCut();

		//============== sliding initialization ====================
		void									DensityField(const std::string& foldername);
		
		void									SlidingSimplify(Visualize &vis, const std::string& foldername ,double tar_percentage, double edge_var_thr, int slidex, int slidey, int slidez);
		
		void									Simplify(int tar_nverts);

		void									ExportData();
		void									ExportData(const std::string& filename);

		
		void									ImproveQuality(Visualize &vis,int odt_count);
		int										ImproveQuality_window(int odt_count, int starts_x, int starts_y, int starts_z);

		inline int								Get_Orig_Num_of_Verts() const { return m_tetmesh->m_storage_verts.size(); }

		void									ResetID();

		TetrahedralMesh*						m_tetmesh;

		bool									is_double_padded_volume;
		double									cleaver_sparsity;
		int										odt_count;

		Cleaver::vec3							shift_vec;

	private:
		LabelVolume*							m_labelvolume;

		//********************************wenhua test***************************************************
		void									SlidingSimplifyMultiProcess(Visualize &vis, const std::string& foldername, TetrahedralMesh& m_tetmesh, double tar_percentage, const int slides[3], const int sizes[3]);
		void									Init_sliding_wenhua();
		bool									render_mode;	

	};
}


#endif
