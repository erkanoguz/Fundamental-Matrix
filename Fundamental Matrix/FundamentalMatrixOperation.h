#ifndef _FUNDAMENTAL_MATRIX_OPERATION_H
#define _FUNDAMENTAL_MATRIX_OPERATION_H

#include "FeatureMatching.h"




class FundamentalMatrixOperation
{
public:
	FundamentalMatrixOperation(Mat& im1, Mat& im2);
	Mat F_Matrix_Eight_Point(const vector<Point2f>&, const vector<Point2f>&);
	Mat F_Matrix_Normalized_Eight_Point(const vector<Point2f>&, const vector<Point2f>&);
	void Plot_Epipolar_Lines(const vector<Point2f>&, const vector<Point2f>&, const Mat&);

private:
	/*Functions*/
	vector<Point3f> augmentVector(const vector<Point2f>&);
	void normalizePixelCoordinates(const vector<Point2f>&, const vector<Point2f>&, vector<Point2f>&, vector<Point2f>&, Mat&, Mat&);
	void vectorMean(const vector<Point3f>&, float&, float&);
	void calculateDistance(const vector<Point3f>&, Mat&, float, float);
	void normalizePoint(const vector<Point3f>, const Mat&, vector<Point2f>&);
	void vector2FMean(const vector <Point2f>&, float&, float&);
	void printMat(const Mat&);

	Mat image1;
	Mat image2;
};




#endif // !_FUNDAMENTAL_MATRIX_OPERATION_H





