#include "FeatureMatching.h"



FeatureMatching::FeatureMatching(Mat& im1, Mat& im2 )
{

	image1 = im1;
	image2 = im2;

	detector = SURF::create();
	keypoints_1, keypoints_2;
	descriptors_1, descriptors_2;
	//-- Step 2: Matching descriptor vectors using FLANN matcher


}

void FeatureMatching::run(vector<Point2f>& p1, vector<Point2f>& p2, bool state)
{
	//-- Step 1: Detect the keypoints using SURF Detector, compute the descriptors
	int minHessian = 400;
	detector->setHessianThreshold(minHessian);

	detector->detectAndCompute(image1, Mat(), keypoints_1, descriptors_1);
	detector->detectAndCompute(image2, Mat(), keypoints_2, descriptors_2);

	//-- Step 2: Matching descriptor vectors using FLANN matcher
	matcher.match(descriptors_1, descriptors_2, matches);


	double max_dist = 0; double min_dist = 100;

	//-- Quick calculation of max and min distances between keypoints
	for (int i = 0; i < descriptors_1.rows; i++)
	{
		double dist = matches[i].distance;
		if (dist < min_dist) min_dist = dist;
		if (dist > max_dist) max_dist = dist;
	}


	//-- Draw only "good" matches (i.e. whose distance is less than 2*min_dist,
	//-- or a small arbitary value ( 0.02 ) in the event that min_dist is very
	//-- small)
	//-- PS.- radiusMatch can also be used here.
	

	for (int i = 0; i < descriptors_1.rows; i++)
	{
		if (matches[i].distance <= max(2 * min_dist, 0.02))
		{
			good_matches.push_back(matches[i]);
		}
	}

	
	drawMatches(image1, keypoints_1, image2, keypoints_2,
		good_matches, img_matches, Scalar::all(-1), Scalar::all(-1),
		vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);


	float x, y;
	/*Coordinates*/
	cout << "Coordiantes: " << endl;

	for (vector<DMatch>::const_iterator it = good_matches.begin(); it != good_matches.end(); ++it)
	{

		x = keypoints_1[it->queryIdx].pt.x; y = keypoints_1[it->queryIdx].pt.y;
		points1.push_back(Point2f(x, y));
		x = keypoints_2[it->trainIdx].pt.x; y = keypoints_2[it->trainIdx].pt.y;
		points2.push_back(Point2f(x, y));

	}


	if (state)
	{
		for (int i = 0; i < 8; i++)
		{
			p1.push_back(points1[i]);
			p2.push_back(points2[i]);
		}
	}
	else
	{
		p1 = points1;
		p2 = points2;

	}
	

}


void FeatureMatching::showResult()
{
	//-- Show detected matches
	imshow("Good Matches", img_matches);

	
	for (int i = 0; i < (int)good_matches.size(); i++)
	{
		printf("-- Good Match [%d] Keypoint 1: %d  -- Keypoint 2: %d  \n", i, good_matches[i].queryIdx, good_matches[i].trainIdx);
	}

}


/*Read all corespondece from file*/
void FeatureMatching::readFeatureFromFile(const char* fileLeftIm, const char* fileRightIm, vector<Point2f>& p1, vector<Point2f>& p2)
{

	p1 = readFile(fileLeftIm);		// Read left image coordinates
	p2 = readFile(fileRightIm);		// Read right image coordinates

}


vector<Point2f> FeatureMatching::readFile(const char* fileName)
{
	
	// Open File
	fstream imagePoints(fileName);
	vector<Point2f> points;

	float tempX, tempY;

	if (!imagePoints.is_open())
	{
		cout << "File can not opened" << endl;
	}

	while (imagePoints)
	{
		imagePoints >> tempX;
		imagePoints >> tempY;
		points.push_back(Point2f(tempX, tempY));

	}

	points.pop_back();

	return points;
}

