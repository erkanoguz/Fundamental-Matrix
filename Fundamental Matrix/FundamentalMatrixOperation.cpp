#include "FundamentalMatrixOperation.h"




FundamentalMatrixOperation::FundamentalMatrixOperation(Mat& im1, Mat& im2)
{
	image1 = im1; image2 = im2;
}

Mat FundamentalMatrixOperation::F_Matrix_Eight_Point(const vector<Point2f>& points1, const vector<Point2f>& points2)
{
	/*Calculation of fundamental Matrix*/

	/*

	┌       ┐ ┌f11  f12  f13┐ ┌u┐
	|u` v` 1|*|f21  f22  f23|*|v|=0
	└       ┘ └f31  f32  f33┘ └1┘

	┌u'1u1   u'1v1   u'1   v'1u1   v'1v1   v'1   u1   v1   1┐   ┌f11┐
	|u'2u2   u'2v2   u'2   v'2u2   v'2v2   v'2   u2   v2   1|   |f12|
	|u'3u3   u'3v3   u'3   v'3u3   v'3v3   v'3   u3   v3   1|   |f13|
	|u'4u4   u'4v4   u'4   v'4u4   v'4v4   v'4   u4   v4   1|   |f21|
	|u'5u5   u'5v5   u'5   v'5u5   v'5v5   v'5   u5   v5   1| * |f22|=0
	|u'6u6   u'6v6   u'6   v'6u6   v'6v6   v'6   u6   v6   1|   |f23|
	|u'7u7   u'7v7   u'7   v'7u7   v'7v7   v'7   u7   v7   1|   |f31|
	└u'8u8   u'8v8   u'8   v'8u8   v'8v8   v'8   u8   v8   1┘   |f32|
																└f33┘
	Ax = 0

	*/



	// Building A Matrix
	Mat A = Mat::zeros(8, 9, CV_32F);
	double u_prime, v_prime, u, v;

	for (int i = 0; i < 8; i++)
	{
		u_prime = points1[i].x;
		v_prime = points1[i].y;
		u = points2[i].x;
		v = points2[i].y;
		A.at<float>(i, 0) = u_prime * u;
		A.at<float>(i, 1) = u_prime * v;
		A.at<float>(i, 2) = u_prime;
		A.at<float>(i, 3) = v_prime * u;
		A.at<float>(i, 4) = v_prime * v;
		A.at<float>(i, 5) = v_prime;
		A.at<float>(i, 6) = u;
		A.at<float>(i, 7) = v;
		A.at<float>(i, 8) = 1;
	}


	// Compute SVD of A Matrix
	Mat U, VT, S;
	SVD::compute(A.clone(), S, U, VT, 4);
	Mat V = VT.t();		//Transpose of VT
	Mat fundamentalMatrix;

	// The entries of F are the components of the column of V corresponding to the least s.v
	Mat F_vec = V.col(V.cols - 1);

	Mat tempF = Mat::zeros(3, 3, CV_32F);
	tempF.at<float>(0, 0) = F_vec.at<float>(0, 0);
	tempF.at<float>(0, 1) = F_vec.at<float>(1, 0);
	tempF.at<float>(0, 2) = F_vec.at<float>(2, 0);
	tempF.at<float>(1, 0) = F_vec.at<float>(3, 0);
	tempF.at<float>(1, 1) = F_vec.at<float>(4, 0);
	tempF.at<float>(1, 2) = F_vec.at<float>(5, 0);
	tempF.at<float>(2, 0) = F_vec.at<float>(6, 0);
	tempF.at<float>(2, 1) = F_vec.at<float>(7, 0);
	tempF.at<float>(2, 2) = F_vec.at<float>(8, 0);


	// To enforce rank 2 constraint

	Mat F_U, F_VT, F_S;
	SVD::compute(tempF.clone(), F_S, F_U, F_VT);

	F_S.at<float>(F_S.rows - 1, 0) = 0;				// Force rank of F to 2

	Mat F_diag = Mat::zeros(3, 3, CV_32F);

	F_diag.at<float>(0, 0) = F_S.at<float>(0, 0);
	F_diag.at<float>(1, 1) = F_S.at<float>(1, 0);
	F_diag.at<float>(2, 2) = F_S.at<float>(2, 0);


	fundamentalMatrix = F_U * F_diag * F_VT;

	return fundamentalMatrix;
}

Mat FundamentalMatrixOperation::F_Matrix_Normalized_Eight_Point(const vector<Point2f>& points1, const vector<Point2f>& points2)
{
	/*Calculation of fundamental Matrix*/

	/*

	┌       ┐ ┌f11  f12  f13┐ ┌u┐
	|u` v` 1|*|f21  f22  f23|*|v|=0
	└       ┘ └f31  f32  f33┘ └1┘

	┌u'1u1   u'1v1   u'1   v'1u1   v'1v1   v'1   u1   v1   1┐   ┌f11┐
	|u'2u2   u'2v2   u'2   v'2u2   v'2v2   v'2   u2   v2   1|   |f12|
	|u'3u3   u'3v3   u'3   v'3u3   v'3v3   v'3   u3   v3   1|   |f13|
	|u'4u4   u'4v4   u'4   v'4u4   v'4v4   v'4   u4   v4   1|   |f21|
	|u'5u5   u'5v5   u'5   v'5u5   v'5v5   v'5   u5   v5   1| * |f22|=0
	|u'6u6   u'6v6   u'6   v'6u6   v'6v6   v'6   u6   v6   1|   |f23|
	|u'7u7   u'7v7   u'7   v'7u7   v'7v7   v'7   u7   v7   1|   |f31|
	└u'8u8   u'8v8   u'8   v'8u8   v'8v8   v'8   u8   v8   1┘   |f32|
																└f33┘
	Ax = 0

	*/

	vector<Point2f> normalizedP1;
	vector<Point2f> normalizedP2;
	Mat norm1 = Mat::zeros(3, 3, CV_32F);
	Mat norm2 = Mat::zeros(3, 3, CV_32F);

	normalizePixelCoordinates(points1, points2, normalizedP1, normalizedP2, norm1, norm2);


	// Building A Matrix
	Mat A = Mat::zeros(normalizedP1.size(), 9, CV_32F);
	float u_prime, v_prime, u, v;


	// show normalizedP1
	cout << "Normalized P1: " << endl;
	for (int i = 0; i < normalizedP1.size(); i++)
	{
		cout << normalizedP1[i] << endl;
	}


	for (int i = 0; i < normalizedP1.size(); i++)
	{
		u_prime = normalizedP2[i].x;
		v_prime = normalizedP2[i].y;
		u = normalizedP1[i].x;
		v = normalizedP1[i].y;
		
		A.at<float>(i, 0) = u_prime * u;
		A.at<float>(i, 1) = u_prime * v;
		A.at<float>(i, 2) = u_prime;
		A.at<float>(i, 3) = v_prime * u;
		A.at<float>(i, 4) = v_prime * v;
		A.at<float>(i, 5) = v_prime;
		A.at<float>(i, 6) = u;
		A.at<float>(i, 7) = v;
		A.at<float>(i, 8) = 1;

	}



	/* Old A matrix
	for (int i = 0; i < 8; i++)
	{
		u_prime = normalizedP1[i].x;
		v_prime = normalizedP1[i].y;
		u = normalizedP2[i].x;
		v = normalizedP2[i].y;
		A.at<float>(i, 0) = u_prime * u;
		A.at<float>(i, 1) = u_prime * v;
		A.at<float>(i, 2) = u_prime;
		A.at<float>(i, 3) = v_prime * u;
		A.at<float>(i, 4) = v_prime * v;
		A.at<float>(i, 5) = v_prime;
		A.at<float>(i, 6) = u;
		A.at<float>(i, 7) = v;
		A.at<float>(i, 8) = 1;

	}
	*/


	cout << "Mat A: " << endl;
	printMat(A);


	// Compute SVD of A Matrix
	Mat U, VT, S;
	SVD::compute(A.clone(), S, U, VT, 4);
	Mat V = VT.t();		//Transpose of VT
	Mat fundamentalMatrix;

	// The entries of F are the components of the column of V corresponding to the least s.v
	Mat F_vec = V.col(V.cols - 1);

	Mat tempF = Mat::zeros(3, 3, CV_32F);
	tempF.at<float>(0, 0) = F_vec.at<float>(0, 0);
	tempF.at<float>(0, 1) = F_vec.at<float>(1, 0);
	tempF.at<float>(0, 2) = F_vec.at<float>(2, 0);
	tempF.at<float>(1, 0) = F_vec.at<float>(3, 0);
	tempF.at<float>(1, 1) = F_vec.at<float>(4, 0);
	tempF.at<float>(1, 2) = F_vec.at<float>(5, 0);
	tempF.at<float>(2, 0) = F_vec.at<float>(6, 0);
	tempF.at<float>(2, 1) = F_vec.at<float>(7, 0);
	tempF.at<float>(2, 2) = F_vec.at<float>(8, 0);



	// To enforce rank 2 constraint

	Mat F_U, F_VT, F_S;
	SVD::compute(tempF.clone(), F_S, F_U, F_VT);

	F_S.at<float>(F_S.rows - 1, 0) = 0;				// Force rank of F to 2

	Mat F_diag = Mat::zeros(3, 3, CV_32F);

	F_diag.at<float>(0, 0) = F_S.at<float>(0, 0);
	F_diag.at<float>(1, 1) = F_S.at<float>(1, 0);
	F_diag.at<float>(2, 2) = F_S.at<float>(2, 0);


	fundamentalMatrix = F_U * F_diag * F_VT;

	norm2 = norm2.t();
	fundamentalMatrix = norm2 * fundamentalMatrix * norm1;

	return fundamentalMatrix;

}

void FundamentalMatrixOperation::normalizePixelCoordinates(const vector<Point2f>& p1, const vector<Point2f>& p2, vector<Point2f>& normVector1, vector<Point2f>& normVector2, Mat& outputMat1, Mat& outputMat2)
{
	
	// augment input points with 1
	vector<Point3f> augVec1 = augmentVector(p1);
	vector<Point3f> augVec2 = augmentVector(p2);


	// Find the centroid of image

	float xMean1, yMean1, xMean2, yMean2;
	vectorMean(augVec1, xMean1, yMean1);		// Compute mean of left image point 
	vectorMean(augVec2, xMean2, yMean2);		// Compute mean of right image points

	//cout << "mean of pts1 x: " << xMean1 << " mean of pts1 y: " << yMean1 << endl;
	//cout << "mean of pts2 x: " << xMean2 << " mean of pts2 y: " << yMean2 << endl;

	// transform matrix to normalize the point
	//Mat outputMat1 = Mat::zeros(3, 3, CV_32F);
	//Mat outputMat2 = Mat::zeros(3, 3, CV_32F);

	calculateDistance(augVec1, outputMat1, xMean1, yMean1);
	calculateDistance(augVec2, outputMat2, xMean2, yMean2);




	normalizePoint(augVec1, outputMat1, normVector1);
	normalizePoint(augVec2, outputMat2, normVector2);


}

vector<Point3f> FundamentalMatrixOperation::augmentVector(const vector<Point2f>& inputPoints)
{

	// Calculate the number of points
	int numberOfPoints = inputPoints.size();

	// agument the vector with 1
	vector <Point3f> outputVector;
	for (int i = 0; i < numberOfPoints; i++)
	{
		outputVector.push_back(Point3f(inputPoints[i].x, inputPoints[i].y, 1));
	}

	return outputVector;
}

void FundamentalMatrixOperation::vectorMean(const vector <Point3f>& inputVector, float& xMean, float& yMean)
{
	// Initializing sum variable
	float sumX = 0;
	float sumY = 0;

	// Sum both x and y points
	for (int i = 0; i < inputVector.size(); i++)
	{
		sumX += (float)inputVector[i].x;
		sumY += (float)inputVector[i].y;
	}

	// Calculate mean
	xMean = sumX / (float)inputVector.size();
	yMean = sumY / (float)inputVector.size();

}

void FundamentalMatrixOperation::vector2FMean(const vector <Point2f>& inputVector, float& xMean, float& yMean)
{
	// Initializing sum variable
	float sumX = 0;
	float sumY = 0;

	// Sum both x and y points
	for (int i = 0; i < inputVector.size(); i++)
	{
		sumX += (float)inputVector[i].x;
		sumY += (float)inputVector[i].y;
	}

	// Calculate mean
	xMean = sumX / (float)inputVector.size();
	yMean = sumY / (float)inputVector.size();

}



void FundamentalMatrixOperation::calculateDistance(const vector<Point3f>& inputVector, Mat& outputMat, float xMean, float yMean)
{
	float tempVal = 0;
	float tempOut = 0;
	float tempSum = 0;
	float allMean = 0;
	float xSum = 0, ySum = 0;

	Mat temp1 = Mat::zeros(3, 3, CV_32F);
	Mat temp2 = Mat::zeros(3, 3, CV_32F);


	vector<Point2f> centered_a;

	for (int i = 0; i < inputVector.size(); i++)
	{
		centered_a.push_back(Point2f(inputVector[i].x - xMean, inputVector[i].y - yMean)); //* (xMean - inputVector[i].x);
		
		//tempVal = xTemp + yTemp;
		//tempOut = sqrt(tempVal);
		//tempSum += tempOut;
	}
	
	float centeredMeanX, centeredMeanY;
	vector2FMean(centered_a, centeredMeanX, centeredMeanY);

	for (int i = 0; i < centered_a.size(); i++)
	{
		xSum += (centered_a[i].x - centeredMeanX) * (centered_a[i].x - centeredMeanX);
		ySum += (centered_a[i].y - centeredMeanY) * (centered_a[i].y - centeredMeanY);
	}

	float varX = xSum / centered_a.size();
	float varY = ySum / centered_a.size();

	float std_x = sqrt(varX);
	float std_y = sqrt(varY);



	//allMean = tempSum / (float)inputVector.size();	// Calculate all Mean

	// implement transformation matrix for normalization
	temp1.at<float>(0, 0) = (float)1 / std_x;
	temp1.at<float>(0, 1) = 0;
	temp1.at<float>(0, 2) = 0;
	temp1.at<float>(1, 0) = 0;
	temp1.at<float>(1, 1) = 1 / std_y;
	temp1.at<float>(1, 2) = 0;
	temp1.at<float>(2, 0) = 0;
	temp1.at<float>(2, 1) = 0;
	temp1.at<float>(2, 2) = 1;

	temp2.at<float>(0, 0) = 1;
	temp2.at<float>(0, 1) = 0;
	temp2.at<float>(0, 2) = -xMean;
	temp2.at<float>(1, 0) = 0;
	temp2.at<float>(1, 1) = 1;
	temp2.at<float>(1, 2) = -yMean;
	temp2.at<float>(2, 0) = 0;
	temp2.at<float>(2, 1) = 0;
	temp2.at<float>(2, 2) = 1;

	outputMat = temp1 * temp2;


	
}





void FundamentalMatrixOperation::normalizePoint(const vector<Point3f> vec, const Mat& tr, vector<Point2f>& outputVec)
{

	// Convert Point3f to Mat
	Mat inputMat = Mat::zeros(3, vec.size(), CV_32F);

	for (int i = 0; i < vec.size(); i++)
	{
		inputMat.at<float>(0, i) = vec[i].x;
		inputMat.at<float>(1, i) = vec[i].y;
		inputMat.at<float>(2, i) = vec[i].z;
	}

	Mat outputMat = tr * inputMat;

	// Convert Mat to Point2F
	for (int i = 0; i < vec.size(); i++)
	{
		float x = outputMat.at<float>(0, i);
		float y = outputMat.at<float>(1, i);
		outputVec.push_back(Point2f(x, y));
	}

}

void FundamentalMatrixOperation::Plot_Epipolar_Lines(const vector<Point2f>& p1, const vector<Point2f>& p2, const Mat& fundamentalMat)
{

	Mat pp1 = Mat_<float>(1, 3);
	Mat pp2 = Mat_<float>(3, 1);
	Mat ur = Mat::zeros(3, 1, CV_32F);
	Mat ul = Mat::zeros(3, 1, CV_32F);



	float px1, px2, AA, BB, CC;
	float py1, py2;
	float testpx, testpy;

	for (int i = 0; i < p1.size(); i++)
	{
		pp1.at<float>(0, 0) = (float)p1[i].x;
		pp1.at<float>(0, 1) = (float)p1[i].y;
		pp1.at<float>(0, 2) = (float)1;

		pp2.at<float>(0, 0) = (float)p2[i].x;
		pp2.at<float>(1, 0) = (float)p2[i].y;
		pp2.at<float>(2, 0) = (float)1;
		
		//LEFT image lines
		ur = fundamentalMat * (pp1.t());
		
		//Right image lines
		ul = (pp2.t()) * fundamentalMat;

		px1 = 0.0;
		px2 = image2.cols - 1;
		AA = ur.at<float>(0, 0);
		BB = ur.at<float>(1, 0);
		CC = ur.at<float>(2, 0);
		py1 = (-CC - (AA * px1)) / BB;
		py2 = (-CC - (AA * px2)) / BB;

		testpx = pp2.at<float>(0, 0);
		testpy = pp2.at<float>(1, 0);

		line(image2, Point(px1, py1), Point(px2, py2), Scalar(255, 0, 0), 1);
		circle(image2, Point(testpx, testpy), 3, Scalar(0, 255, 0), 1);



		px1 = 0.0;
		px2 = image2.cols - 1;
		AA = ul.at<float>(0, 0);
		BB = ul.at<float>(0, 1);
		CC = ul.at<float>(0, 2);
		py1 = (-CC - (AA * px1)) / BB;
		py2 = (-CC - (AA * px2)) / BB;

		testpx = pp1.at<float>(0, 0);
		testpy = pp1.at<float>(0, 1);

		line(image1, Point(px1, py1), Point(px2, py2), Scalar(255, 0, 0), 1);
		circle(image1, Point(testpx, testpy), 3, Scalar(0, 255, 0), 1);

	}

	imshow("LEFT IMAGE EPIPOLAR LINES", image2);
	imshow("RIGHT IMAGE EPIPOLAR LINES", image1);

}


void FundamentalMatrixOperation::printMat(const Mat& inputMat)
{

	for (int i = 0; i < inputMat.rows; i++)
	{
		for (int j = 0; j < inputMat.cols; j++)
		{
			cout << inputMat.at<float>(i, j) << " ";
		}
		cout << endl;
	}
}