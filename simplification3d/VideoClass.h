#ifndef _VIDEO_CLASS_H_
#define _VIDEO_CLASS_H_

#include <cwgUtilities.h>

namespace cwg
{
	class VideoBase
	{
	public:
		VideoBase() : m_width(-1), m_height(-1), m_nframes(-1) {}
		virtual ~VideoBase() {};

		virtual	inline cv::Mat					operator[](int k) const = 0;
		virtual inline cv::Mat					voxel(int k, int i, int j) const = 0;
		virtual inline cv::Vec3b				operator()(int k, int i, int j) const = 0;

		virtual inline cv::Vec3d				volume_data(int k, int i, int j) const = 0;

		virtual inline int						Height() const {return m_height;}
		virtual inline int						Width() const {return m_width;}
		virtual inline int						nFrames() const {return m_nframes;}
		virtual inline int						nVoxels() const {return m_width*m_height*m_nframes;}

	protected:
		int										m_width;
		int										m_height;
		int										m_nframes;
	};

	class VideoMEM : public VideoBase			
	{
	public:
		VideoMEM(int N);
		VideoMEM();
		virtual ~VideoMEM(void);

		virtual void							Load_from_Folder(const std::string& foldername, const std::string& ext);
		virtual void							Write_to_Folder(const std::string& foldername, const std::string& ext, int flag=0);

		virtual inline void						SetFrame(cv::Mat frame, int k){mVideo[k] = frame; m_height = frame.rows; m_width = frame.cols;}
		
		virtual cv::Mat							operator[](int k) const {return mVideo[k];}
		virtual inline cv::Mat					voxel(int k, int i, int j) const { return mVideo[k].row(i).col(j);}
		virtual inline cv::Vec3b				operator()(int k, int i, int j) const { return (*this)[k].at<cv::Vec3b>(i,j); }
		void									convert_to_double(int type);
		virtual inline cv::Vec3d				volume_data(int k, int i, int j) const { return (*this)[k].at<cv::Vec3d>(i,j); }

		virtual inline int						nChannels() const {return mVideo[0].channels();}
		virtual inline int						Depth()     const {return mVideo[0].depth();}

		virtual void							ConvertTo(int rtype, double alpha=1.0, double beta=0.0);
		virtual int								State(); // Complete | Incomplete | Empty

		virtual cv::Mat							DataTable();
		virtual void							SetData(cv::Mat mat_data, int height, int width, int nframes);
	
		virtual void							generate_single_label(int h, int w, int nframes) {}
		virtual void							generate_quadratic_field(int h, int w, int nframes) {}

	protected:
		std::vector<cv::Mat>					mVideo;
	};

	class PaddedVideoMEM : public VideoMEM
	{
	public:
		PaddedVideoMEM(int nPaddings);
		~PaddedVideoMEM();

		void									Load_from_Folder(const std::string& foldername, const std::string& ext);
		void									generate_single_label(int h, int w, int nframes);
		void									generate_quadratic_field(int h, int w, int nframes);

		inline int								nPaddings() const { return m_nPaddings; }			

		inline void								SetFrame(cv::Mat frame, int k){ throw std::exception("ERROR CALLED!"); }
		inline void								SetData(cv::Mat mat_data, int height, int width, int nframes) { throw std::exception("ERROR CALLED!"); }

		inline void								ForceOdd(bool odd) { m_setting_odd = odd; }

	private:
		void									load_from_folder_odd(const std::vector<std::string>& filelist);
		void									load_from_folder_orig(const std::vector<std::string>& filelist);

	private:
		int										m_nPaddings;
		bool									m_setting_odd;
	};
}

#endif


