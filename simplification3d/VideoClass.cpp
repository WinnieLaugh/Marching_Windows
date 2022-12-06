#include "StdAfx.h"
#include "VideoClass.h"
#include <Objbase.h>

using namespace std;
using namespace cv;

cwg::VideoMEM::VideoMEM( int N )
{
	mVideo = vector<Mat>(N);
	m_nframes = N;
}

cwg::VideoMEM::VideoMEM()
{

}

cwg::VideoMEM::~VideoMEM( void )
{

}

int cwg::VideoMEM::State()
{// output: 1: complete; 0: incomplete; -1: empty;
	if (mVideo.size() == 0)
		return -1;

	bool isempty = true;
	bool iscomplete = true;

	for (int k=0; k<mVideo.size(); k++)
	{
		isempty = isempty && mVideo[k].empty();
		iscomplete = iscomplete && !(mVideo[k].empty());
	}

	if (isempty)
		return -1;

	if (iscomplete)
		return 1;

	return 0;
}
cv::Mat cwg::VideoMEM::DataTable()
{
	int numRows = nVoxels();
	int numCols = nChannels();

	int height = Height();
	int width = Width();

	Mat data_table(numRows, numCols, CV_64FC1, Scalar(0.0));

	int offset = 0;
	for (int k=0; k<mVideo.size(); k++)
	{
		Mat tmp = mVideo[k];
		tmp = tmp.reshape(1, height*width);

		Mat tmp_double;
		if (tmp.depth() != CV_64FC1)
			tmp.convertTo(tmp_double, CV_64FC1);
		else
			tmp_double = tmp;

		tmp_double.copyTo(data_table.rowRange(offset, offset+height*width));
		offset += height*width;
	}
	return data_table;
}

void cwg::VideoMEM::SetData( Mat mat_data, int height, int width, int nframes )
{
	assert(mat_data.total() == height*width*nframes);

	if (mat_data.rows == 1 && mat_data.cols != 1)
		mat_data = mat_data.t();

	mVideo = vector<Mat>(nframes);

	int offset = 0;
	Mat tmp;
	for (int k=0; k<nframes; k++)
	{
		tmp = mat_data.rowRange(offset, offset+height*width);
		tmp = tmp.reshape(mat_data.channels(), height);

		mVideo[k] = Mat(height, width, mat_data.type());
		tmp.copyTo(mVideo[k]);

		offset += height*width;
	}
}

void cwg::VideoMEM::Load_from_Folder( const string& foldername, const string& ext )
{
	vector<string> filelist = show_file_list(foldername, ext);
	mVideo = vector<Mat>(filelist.size());

	for (int k=0; k<filelist.size(); k++)
	{
		Mat tmp;
		if (ext == "jpg" || ext == "png" || ext == "bmp" || ext == "tiff")
		{
			tmp = imread(foldername+filelist[k]);
		}
		else
		{
			tmp = loadmat(foldername+filelist[k]);
		}

		mVideo[k] = tmp;
	}

	m_height = mVideo[0].rows;
	m_width = mVideo[0].cols;
	m_nframes = mVideo.size();
	
}

void cwg::VideoMEM::ConvertTo( int rtype, double alpha/*=1.0*/, double beta/*=0.0*/ )
{
	for (int k=0; k<mVideo.size(); k++)
	{
		if (mVideo[k].empty())
			continue;

		mVideo[k].convertTo(mVideo[k], rtype, alpha, beta);
	}
}

void cwg::VideoMEM::Write_to_Folder( const string& foldername, const string& ext, int flag/*=0*/ )
{//flag: 0 ==> images; 1 ==> videos
	if (flag == 0 && (ext == "jpg" || ext == "png" || ext == "tiff" || ext == "bmp"))
	{
		char str[300];
		string filename;
		for (int k=0; k<mVideo.size(); k++)
		{
			if (mVideo[k].empty())
				continue;
			
			sprintf(str, "frame%04d.", k);
			filename = foldername + str + ext;

			imwrite(filename, mVideo[k]);
		}
	}
	else if (flag == 1 && (ext == "mp4" || ext == "avi" || ext == "wmv"))
	{
		throw std::exception("cwg::VideoMEM::Write_to_Folder error, not implemented!");
	}
	else
		throw std::exception("cwg::VideoMEM::Write_to_Folder error");

}

void cwg::VideoMEM::convert_to_double( int type )
{
	if (mVideo[0].type() == CV_32SC1)
		throw std::exception("IntVideo cannot be converted!");
	
	for (int i=0; i<m_nframes; i++)
		mVideo[i].convertTo(mVideo[i], type);
}

cwg::PaddedVideoMEM::PaddedVideoMEM( int nPaddings ) : m_setting_odd(false)
{
	m_nPaddings = nPaddings;
}

cwg::PaddedVideoMEM::~PaddedVideoMEM()
{

}

void cwg::PaddedVideoMEM::Load_from_Folder( const string& foldername, const string& ext )
{
	vector<string> filelist = show_file_list(foldername, ext);
	for (int i=0; i<filelist.size(); i++)
		filelist[i] = foldername + filelist[i];
	
	if (m_setting_odd)
		load_from_folder_odd(filelist);
	else
		load_from_folder_orig(filelist);

	
}

void cwg::PaddedVideoMEM::generate_single_label( int h, int w, int nframes )
{
	mVideo = vector<Mat>(nframes + m_nPaddings*2);

	for (int k=0; k<nframes; k++)
	{
		Mat tmp;

		tmp = Mat::zeros(h, w, CV_32SC1);
		copyMakeBorder(tmp, tmp, m_nPaddings, m_nPaddings, m_nPaddings, m_nPaddings, BORDER_CONSTANT, Scalar::all(-1));

		mVideo[k+m_nPaddings] = tmp;
	}

	m_height  = mVideo[m_nPaddings].rows;
	m_width   = mVideo[m_nPaddings].cols;
	m_nframes = mVideo.size();

	for (int k=0; k<m_nPaddings; k++)
	{
		Mat tmp1, tmp2;
		tmp1 = Mat::ones(m_height, m_width, CV_32SC1)*(-1);
		tmp2 = tmp1.clone();
		mVideo[k] = tmp1;
		mVideo[mVideo.size()-1-k] = tmp2;
	}
}

void cwg::PaddedVideoMEM::generate_quadratic_field( int h, int w, int nframes )
{
	mVideo = vector<Mat>(nframes + m_nPaddings*2);

	for (int k=0; k<nframes; k++)
	{
		Mat tmp;

		tmp = Mat::zeros(h, w, CV_64FC3);
		double intensity = 0;
		for (int i=0; i<h; i++)
		{
			for (int j=0; j<w; j++)
			{
				intensity = i*i+j*j+k*k;
				tmp.at<Vec3d>(i, j) = Vec3d(intensity, intensity, intensity);
			}
		}
		copyMakeBorder(tmp, tmp, m_nPaddings, m_nPaddings, m_nPaddings, m_nPaddings, BORDER_CONSTANT, Scalar::all(0));

		mVideo[k+m_nPaddings] = tmp;
	}

	m_height  = mVideo[m_nPaddings].rows;
	m_width   = mVideo[m_nPaddings].cols;
	m_nframes = mVideo.size();

	for (int k=0; k<m_nPaddings; k++)
	{
		Mat tmp1, tmp2;
		tmp1 = Mat::zeros(m_height, m_width, CV_64FC3);
		tmp2 = tmp1.clone();
		mVideo[k] = tmp1;
		mVideo[mVideo.size()-1-k] = tmp2;
	}
}

void cwg::PaddedVideoMEM::load_from_folder_odd( const std::vector<std::string>& filelist )
{
	int bWidth_odd, bHeight_odd, bFrames_odd;
	bFrames_odd = filelist.size() % 2;

	mVideo = vector<Mat>(filelist.size() + m_nPaddings*2 + 1 - bFrames_odd);
	
	if (filelist.size() == 0)
		return;

	unsigned found = filelist[0].find_last_of(".");
	string ext = filelist[0].substr(found+1);
	
	for (int k=0; k<filelist.size(); k++)
	{
		Mat tmp;
		if (ext == "jpg" || ext == "png" || ext == "bmp" || ext == "tiff")
		{
			tmp = imread(filelist[k]);
			if (k == 0)
			{
				bWidth_odd = tmp.cols % 2;
				bHeight_odd = tmp.rows % 2;
			}
			copyMakeBorder(tmp, tmp, m_nPaddings, m_nPaddings+1-bHeight_odd, 
				m_nPaddings, m_nPaddings+1-bWidth_odd, BORDER_CONSTANT, Scalar::all(0));
		}
		else
		{
			tmp = loadmat(filelist[k]);
			if (k == 0)
			{
				bWidth_odd = tmp.cols % 2;
				bHeight_odd = tmp.rows % 2;
			}
			copyMakeBorder(tmp, tmp, m_nPaddings, m_nPaddings+1-bHeight_odd,
				m_nPaddings, m_nPaddings+1-bWidth_odd, BORDER_CONSTANT, Scalar::all(-1));
		}

		if (bWidth_odd==0)
		{
			int lastcol = tmp.cols - m_nPaddings - 2;
			tmp.col(lastcol).copyTo(tmp.col(lastcol+1));
		}
		if (bHeight_odd==0)
		{
			int lastrow = tmp.rows - m_nPaddings - 2;
			tmp.row(lastrow).copyTo(tmp.row(lastrow+1));
		}

		mVideo[k+m_nPaddings] = tmp;
	}

	m_height = mVideo[m_nPaddings].rows;
	m_width = mVideo[m_nPaddings].cols;
	m_nframes = mVideo.size();

	for (int k=0; k<m_nPaddings; k++)
	{
		Mat tmp1, tmp2;
		if (ext == "jpg" || ext == "png" || ext == "bmp" || ext == "tiff")
			tmp1 = Mat::zeros(m_height, m_width, mVideo[m_nPaddings].type());
		else
			tmp1 = Mat::ones(m_height, m_width, mVideo[m_nPaddings].type())*(-1);
		tmp2 = tmp1.clone();
		mVideo[k] = tmp1;
		mVideo[mVideo.size()-1-k] = tmp2;
	}

	if (bFrames_odd == 0)
		mVideo[m_nPaddings+filelist.size()] = mVideo[m_nPaddings+filelist.size()-1].clone();
}

void cwg::PaddedVideoMEM::load_from_folder_orig( const std::vector<std::string>& filelist )
{
	mVideo = vector<Mat>(filelist.size() + m_nPaddings*2);

	if (filelist.size() == 0)
		return;
	unsigned found = filelist[0].find_last_of(".");
	string ext = filelist[0].substr(found+1);

	for (int k=0; k<filelist.size(); k++)
	{
		Mat tmp;
		if (ext == "jpg" || ext == "png" || ext == "bmp" || ext == "tiff")
		{
			tmp = imread(filelist[k]);
			copyMakeBorder(tmp, tmp, m_nPaddings, m_nPaddings, m_nPaddings, m_nPaddings, BORDER_CONSTANT, Scalar::all(0));
		}
		else
		{
			tmp = loadmat(filelist[k]);
			copyMakeBorder(tmp, tmp, m_nPaddings, m_nPaddings, m_nPaddings, m_nPaddings, BORDER_CONSTANT, Scalar::all(-1));
		}
		mVideo[k+m_nPaddings] = tmp;
	}

	m_height = mVideo[m_nPaddings].rows;
	m_width = mVideo[m_nPaddings].cols;
	m_nframes = mVideo.size();

	for (int k=0; k<m_nPaddings; k++)
	{
		Mat tmp1, tmp2;
		if (ext == "jpg" || ext == "png" || ext == "bmp" || ext == "tiff")
			tmp1 = Mat::zeros(m_height, m_width, mVideo[m_nPaddings].type());
		else
			tmp1 = Mat::ones(m_height, m_width, mVideo[m_nPaddings].type())*(-1);
		tmp2 = tmp1.clone();
		mVideo[k] = tmp1;
		mVideo[mVideo.size()-1-k] = tmp2;
	}
}
