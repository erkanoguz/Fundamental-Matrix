
#ifndef FEATURE_MATCHING_H
#define FEATURE_MATCHING_H


#include <opencv2/highgui.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/videoio.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/xfeatures2d.hpp>


using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;


#define NORMAL 0
#define EIGHTPOINT 1


class FeatureMatching
{
public:
	FeatureMatching(Mat& image1, Mat& image2);
	~FeatureMatching() {}
	void run(vector<Point2f>&, vector<Point2f>&, bool state = NORMAL);
	void showResult();
	void readFeatureFromFile(const char*, const char*, vector<Point2f>&, vector<Point2f>&);

private:
	Ptr<SURF> detector;
	vector<KeyPoint> keypoints_1, keypoints_2;
	Mat descriptors_1, descriptors_2;
	Mat image1;
	Mat image2;
	FlannBasedMatcher matcher;
	vector<DMatch> matches;
	vector<Point2f> points1; //(good_matches.size());
	vector<Point2f> points2; //(good_matches.size());
	vector< DMatch > good_matches;
	Mat img_matches;


	/*Functions*/
	vector<Point2f> readFile(const char*);
};




#endif // !FEATURE_MATCHING_H


