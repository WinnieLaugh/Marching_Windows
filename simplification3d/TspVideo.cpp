#include "StdAfx.h"
#include "TSPVideo.h"

using namespace std;
using namespace cv;

cwg::TSPVideo::TSPVideo(void)
{
	m_nTsp = -1;
	m_TspVideoLite = new VideoMEM();
}

cwg::TSPVideo::~TSPVideo(void)
{
	SAFE_DELETE(m_TspVideoLite);
}

void cwg::TSPVideo::Load_from_Folder( const string& foldername )
{
	m_TspVideoLite->Load_from_Folder(foldername, "svdata");
	compute_nTsp();
}

void cwg::TSPVideo::compute_nTsp()
{
	Mat tmp;
	double minv, maxv, MaxV=0;
	for (int k=0; k<m_TspVideoLite->nFrames(); k++)
	{
		tmp = (*m_TspVideoLite)[k];
		minMaxLoc(tmp, &minv, &maxv);
		if (maxv >= MaxV)
			MaxV = maxv;
	}
	m_nTsp = MaxV + 1;
}

cwg::PaddedTSPVideo::PaddedTSPVideo( void ) : nPaddings(4)
{
	SAFE_DELETE(m_TspVideoLite);
	m_TspVideoLite = (cwg::VideoMEM*)new cwg::PaddedVideoMEM(nPaddings);
}

cwg::PaddedTSPVideo::PaddedTSPVideo( int npaddings ) : nPaddings(npaddings)
{
	SAFE_DELETE(m_TspVideoLite);
	m_TspVideoLite = (cwg::VideoMEM*)new cwg::PaddedVideoMEM(nPaddings);
}

void cwg::PaddedTSPVideo::Load_from_Folder( const string& foldername )
{
	cwg::PaddedVideoMEM* tmpptr = (cwg::PaddedVideoMEM*)m_TspVideoLite;
	if (nPaddings == 4)
		tmpptr->ForceOdd(true);
	else if (nPaddings == 2)
		tmpptr->ForceOdd(false);
	else
		throw exception("void cwg::PaddedTSPVideo::Load_from_Folder( const string& foldername ): nPaddings != 2 or 4");
		
	m_TspVideoLite->Load_from_Folder(foldername, "svdata");
	compute_nTsp();
	for (int k=0; k<m_TspVideoLite->nFrames(); k++)
	{
		for (int i=0; i<m_TspVideoLite->Height(); i++)
		{
			for (int j=0; j<m_TspVideoLite->Width(); j++)
			{
				if ( (*this)(k, i, j) == -1 )
					(*m_TspVideoLite)[k].at<int>(i, j) = m_nTsp;
			}
		}
	}
	m_nTsp += 1;
}

void cwg::PaddedTSPVideo::GenerateSimpleTSPVideo( int h, int w, int nframes )
{
	m_TspVideoLite->generate_single_label(h, w, nframes);
	compute_nTsp();
	for (int k=0; k<m_TspVideoLite->nFrames(); k++)
	{
		for (int i=0; i<m_TspVideoLite->Height(); i++)
		{
			for (int j=0; j<m_TspVideoLite->Width(); j++)
			{
				if ( (*this)(k, i, j) == -1 )
					(*m_TspVideoLite)[k].at<int>(i, j) = m_nTsp;
			}
		}
	}
	m_nTsp += 1;
}

