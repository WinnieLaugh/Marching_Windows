#ifndef _CWG_TSPVIDEO_H_
#define _CWG_TSPVIDEO_H_

#ifdef SIMPLIFICATION3D_EXPORTS
#define SIMPLIFICATION3D_API __declspec(dllexport)
#else
#define SIMPLIFICATION3D_API __declspec(dllimport)
#endif

#include <cwgUtilities.h>
#include "VideoClass.h"

namespace cwg
{
	class SIMPLIFICATION3D_API TSPVideo
	{
	public:
		TSPVideo(void);
		virtual ~TSPVideo(void);

		virtual void										Load_from_Folder(const std::string& foldername);
		
		inline int											nTsp() const {return m_nTsp;}
		inline int											Width() const {return m_TspVideoLite->Width();}
		inline int											Height() const {return m_TspVideoLite->Height();}
		inline int											nFrames() const {return m_TspVideoLite->nFrames();}
		inline int											operator()(int k, int i, int j) const {return (*m_TspVideoLite)[k].at<int>(i,j);}

		cv::Mat												frame(int k) const { return (*m_TspVideoLite)[k]; }

	protected:
		void												compute_nTsp();
		
		VideoMEM*											m_TspVideoLite;

		int													m_nTsp;
	private:
		
	};

	class SIMPLIFICATION3D_API PaddedTSPVideo : public TSPVideo
	{
	public:
		PaddedTSPVideo(void);
		PaddedTSPVideo(int npaddings);
		~PaddedTSPVideo(void) {}

		void												Load_from_Folder(const std::string& foldername);
		void												GenerateSimpleTSPVideo(int h, int w, int nframes);

		const int											nPaddings;
	};
}

typedef cwg::PaddedTSPVideo LabelVolume;

#endif


