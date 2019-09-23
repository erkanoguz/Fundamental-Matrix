#include <opencv2/highgui.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>
#include <opencv2/features2d.hpp>
#include <opencv2/videoio.hpp>
#include <opencv2/opencv.hpp>
#include <vector>
#include <opencv2/xfeatures2d.hpp>
#include <iomanip>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>
#include "FeatureMatching.h"
#include "FundamentalMatrixOperation.h"
#include <fstream>

using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;
using namespace std::chrono;

int computeRank(const Mat1d&);
void showDiffViewImage(Mat , Mat , vector<Point2f>& , vector<Point2f>& );
void printMat(const Mat& );

int main()
{

	Mat image = imread("Data/pic_a.jpg");
	Mat image2 = imread("Data/pic_b.jpg");

	vector<Point2f> pts1;
	vector<Point2f> pts2;

	const char* leftImFileName = "Data/pts2d-pic_a.txt";
	const char* rightImFileName = "Data/pts2d-pic_b.txt";


	FeatureMatching fm(image, image2);
	fm.readFeatureFromFile(leftImFileName, rightImFileName, pts1, pts2);


	Mat testPoint1 = Mat::zeros(1, 3, CV_32F);
	Mat testPoint2 = Mat::zeros(3, 1, CV_32F);

	testPoint1.at<float>(0, 0) = pts1[0].x;
	testPoint1.at<float>(0, 1) = pts1[0].y;
	testPoint1.at<float>(0, 2) = 1;

	testPoint2.at<float>(0, 0) = pts2[0].x;
	testPoint2.at<float>(1, 0) = pts2[0].y;
	testPoint2.at<float>(2, 0) = 1;


	FundamentalMatrixOperation fundamentalMatOp(image, image2);
	Mat fundamental_matrix = fundamentalMatOp.F_Matrix_Normalized_Eight_Point(pts1, pts2);
		
	cout << "Fundamental MAtrix: " << endl;
	printMat(fundamental_matrix);
	fundamentalMatOp.Plot_Epipolar_Lines(pts1, pts2, fundamental_matrix);
		
	// Combined Two Image
	//showDiffViewImage(image, image2, pts1, pts2);


	/************Feature Matcing*************/
	/*
	Mat image = imread("Data/im1.jpg");
	Mat image2 = imread("Data/im2.jpg");

	fm.run(pts1, pts2, EIGHTPOINT);

	testPoint1.at<float>(0, 0) = pts1[0].x;
	testPoint1.at<float>(0, 1) = pts1[0].y;
	testPoint1.at<float>(0, 2) = 1;

	testPoint2.at<float>(0, 0) = pts2[0].x;
	testPoint2.at<float>(1, 0) = pts2[0].y;
	testPoint2.at<float>(2, 0) = 1;

	FundamentalMatrixOperation fundamentalMatOp2(image, image2);
	fundamental_matrix = fundamentalMatOp2.F_Matrix_Normalized_Eight_Point(pts1, pts2);

	cout << "Fundamental MAtrix2: " << endl;
	printMat(fundamental_matrix);
	fundamentalMatOp2.Plot_Epipolar_Lines(pts1, pts2, fundamental_matrix);

	// Combined Two Image
	//showDiffViewImage(image, image2, pts1, pts2);
	*/
	waitKey(0);

}


void showDiffViewImage(Mat image2, Mat image, vector<Point2f>& pts1, vector<Point2f>& pts2)
{
	/*Draw 2 images in single window*/

	int rowsLeft = image.rows;
	int rowsRight = image2.rows;
	int colsLeft = image.cols;
	int colsRight = image2.cols;
	int allCols = colsLeft + colsRight;
	int ch = image.channels();

	Mat combinedIm = Mat::zeros(rowsLeft, allCols, image.type());

	for (int i = 0; i < rowsLeft; i++)
	{
		for (int j = 0; j < colsLeft; j++)
		{
			for (int c = 0; c < ch; c++)
				combinedIm.at<Vec3b>(i, j)[c] = image.at<Vec3b>(i, j)[c];
		}
	}

	for (int i = 0; i < rowsLeft; i++)
	{
		for (int j = colsLeft; j < allCols; j++)
		{
			for (int c = 0; c < ch; c++)
				combinedIm.at<Vec3b>(i, j)[c] = image2.at<Vec3b>(i, j - colsLeft)[c];
		}
	}


	for (int i = 0; i < pts1.size(); i++)
	{
		int x1 = (int)pts1[i].x;
		int y1 = (int)pts1[i].y;
		int x2 = (int)pts2[i].x;
		int y2 = (int)pts2[i].y;

		circle(combinedIm, Point(x1, y1), 3, Scalar(0, 0, 255), 3, 8, 0);
		circle(combinedIm, Point(colsLeft + x2, y2), 3, Scalar(0, 255, 0), 3, 8, 0);


		line(combinedIm, Point(x1, y1), Point(x2 + colsLeft, y2), Scalar(0, 255, 0), 1, 8, 0);

	}
	imshow("Combined Image", combinedIm);
}

int computeRank(const Mat1d& input) 
{
	Mat1d w, u, vT;

	// Compute SVD
	SVD::compute(input, w, u, vT);


	Mat1b nonZeroSingularValues = w > 0.0001;

	int rank = countNonZero(nonZeroSingularValues);

	return rank;	// input matrix rank
}

void printMat(const Mat& input)
{
	for (int i = 0; i < input.rows; i++)
	{
		for (int j = 0; j < input.cols; j++)
		{
			cout << input.at<float>(i, j) << " ";
		}
		cout << endl;
	}
}