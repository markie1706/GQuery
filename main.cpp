/*
 main.cpp
 QForce
 Purpose: Generates a table of polar coordinates of spatial points around some pre-defined anchor point from random latitude and longitude values. Then query this table of polar coordinates for those within a query radius of a provided latitude and longitude, which must first be converted into polar coordinates relative to the same anchor as the generated spatial points.
 
 Created by Mark Keane on 8/3/16.
 Copyright © 2016 Mark Keane. All rights reserved.
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <ios>
#include <iterator>
#include <array>
#include <sstream>
#include <iomanip>
using namespace std;
using  ns = chrono::nanoseconds;
using get_time = chrono::steady_clock ;
#define PI 3.14159265
#define pi 3.14159265358979323846
#define earthRadiusKm 6371.0
#define maxDistanceFromAnAnchorPoint 9000

// This function converts decimal degrees to radians
double deg2rad(double deg) {
	return (deg * pi / 180);
}

//  This function converts radians to decimal degrees
double rad2deg(double rad) {
	return (rad * 180 / pi);
}


double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
	double lat1r, lon1r, lat2r, lon2r, u, v;
	lat1r = deg2rad(lat1d);
	lon1r = deg2rad(lon1d);
	lat2r = deg2rad(lat2d);
	lon2r = deg2rad(lon2d);
	u = sin((lat2r - lat1r)/2);
	v = sin((lon2r - lon1r)/2);
	return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

inline double round( double val )
{
	if( val < 0 ) return ceil(val - 0.5);
	return floor(val + 0.5);
}

double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

const int rows = 500000;
const int columns = 2;

static double latsAndLongsOfDataPoints [rows][columns];

static int arrayOfPoints [rows][columns];
static double arrayOfPointsAngles [rows][columns];
//Accuracy
static double arrayOfPointsLatLongs [rows][columns];



int main(int argc, const char * argv[]) {
	//Gets random number generation initiated correctly:
	srand ( time(NULL) );

	const int databaseAnchorRows = 5;
	const int databaseAnchorsColumns = 2;
	//databaseAnchorsAsIntegers have latitude and longitude values for anchor points in predefined geographic places. These are the values that new and query geo points are compared too, with the new or query points being converted into polar coordinates based on the anchor point they are closest too. Ints are used to make figuring which anchor point a new or query point is closest too faster.
	int databaseAnchorsAsIntegers [databaseAnchorRows][databaseAnchorsColumns] = {
		{42, -120}, //West north america
		{41, -97},	//Central north america
		{41, -78},	//East north america
		{47, 4},	//Western europe
		{51, 30}	//Eastern europe
	};
	
	//----------------------------------------------------------------------------------------------------------------------------
	//-----		-1) Pre Setup. We need to create a number of random latitude and longitude pairs to be converted to polar coordinates around the nearest anchor point to them, which are then queried. The number of random points to create and query is defined in the variable "rows"------
	//----------------------------------------------------------------------------------------------------------------------------
	double _minLatitude = 34.00000;
	double _maxLatitude = 48.00000;
	double _minLongitude = -125.00000;
	double _maxLongitude = -115.00000;
	for (int v = 0 ; v < rows ; v++) {
		latsAndLongsOfDataPoints[v][0] = round(fRand(_minLatitude, _maxLatitude) * 10000.0) / 10000.0;
		latsAndLongsOfDataPoints[v][1] = round(fRand(_minLongitude, _maxLongitude) * 10000.0) / 10000.0;
	}
	
	//Convert the randomly generated lats and longs to polar coordinates based on the nearest anchor to them
	for (int v = 0 ; v < rows ; v++) {
		int radiusOfPointToAdd = distanceEarth((double)databaseAnchorsAsIntegers[0][0], (double)databaseAnchorsAsIntegers[0][1], latsAndLongsOfDataPoints[v][0], latsAndLongsOfDataPoints[v][1]);
		
		//s1 and s2 are the two sides of the triangle (third side is the radiusOfTheQueryPointToAnchor)
		int s1 = distanceEarth((double)databaseAnchorsAsIntegers[0][0], (double)databaseAnchorsAsIntegers[0][1], (double)databaseAnchorsAsIntegers[0][0], latsAndLongsOfDataPoints[v][1]);
		int s2 = distanceEarth(latsAndLongsOfDataPoints[v][0], latsAndLongsOfDataPoints[v][1], (double)databaseAnchorsAsIntegers[0][0], latsAndLongsOfDataPoints[v][1]);
		//Because distanceEarth gives an absolute value response, negative signs are added here if needed.
		if (latsAndLongsOfDataPoints[v][1] < databaseAnchorsAsIntegers[0][1]) {
			s1 = s1 - (s1 * 2);
		}
		if (latsAndLongsOfDataPoints[v][1] < databaseAnchorsAsIntegers[0][0]) {
			s2 = s2 - (s2 * 2);
		}
		if (s1 == 0) {
			s1 = 1;
		}
		double thetaOfPointToAdd = (atan(s2/s1) * 180 / PI);
		thetaOfPointToAdd = round(thetaOfPointToAdd * 1000.0) / 1000.0;
		if (databaseAnchorsAsIntegers[0][0] < latsAndLongsOfDataPoints[v][0] && latsAndLongsOfDataPoints[v][1] < databaseAnchorsAsIntegers[0][1]) {
			thetaOfPointToAdd = 180.0 - ((thetaOfPointToAdd + 90.0) - 180.0);
		} else if (databaseAnchorsAsIntegers[0][0] > latsAndLongsOfDataPoints[v][0] && latsAndLongsOfDataPoints[v][1] < databaseAnchorsAsIntegers[0][1]) {
			thetaOfPointToAdd = 180.0 - ((thetaOfPointToAdd + 90.0) - 180.0);
		} else if(databaseAnchorsAsIntegers[0][0] > latsAndLongsOfDataPoints[v][0] && latsAndLongsOfDataPoints[v][1] > databaseAnchorsAsIntegers[0][1]) {
			thetaOfPointToAdd = 180.0 - ((thetaOfPointToAdd + 270.0) - 180.0);
		}
		
		arrayOfPoints[v][0] = radiusOfPointToAdd;
		arrayOfPointsAngles[v][0] = thetaOfPointToAdd;
		
		//Accuracy
		arrayOfPointsLatLongs[v][0] = latsAndLongsOfDataPoints[v][0];
		arrayOfPointsLatLongs[v][1] = latsAndLongsOfDataPoints[v][1];
	}
	
	
	//Convert arrayOfPoints and arrayOfPointsAngles into ascending order
	for (int i = 0; i < rows; ++i)
	{
		for (int j = i + 1; j < rows; ++j)
		{
			if (arrayOfPoints[i][0] > arrayOfPoints[j][0])
			{
				int a =  arrayOfPoints[i][0];
				arrayOfPoints[i][0] = arrayOfPoints[j][0];
				arrayOfPoints[j][0] = a;
				
				double b =  arrayOfPointsAngles[i][0];
				arrayOfPointsAngles[i][0] = arrayOfPointsAngles[j][0];
				arrayOfPointsAngles[j][0] = b;
			}
		}
	}
	
	//----------------------------------------------------------------------------------------------------------------------------
	//-----		0) Setup. Above creates random sets of spatial points in polar coordinates. Below sets up anchor points, timers and inputs to the algo.		------
	//----------------------------------------------------------------------------------------------------------------------------
	//Timers to measure how fast the algo goes.
	auto begin = chrono::high_resolution_clock::now();
	
	//INPUTS TO ALGORITHM, a lat and long value and number of kilometers around it you want to query.  qr is query radius. algorithm should return the randomly generated points in step -1 within qr km of queryPointLatitude and queryPointLongitude
	int queryPointLatitude = 37;
	int queryPointLongitude = -119;
	double queryPointLatDouble = 37.317752;
	double queryPointLongDouble = -118.904297;
	int qr = 25;
	
	//----------------------------------------------------------------------------------------------------------------------------
	//-----		1) Figure out which anchor point and associated table of spatial points is closest to the query point		------
	//----------------------------------------------------------------------------------------------------------------------------
	int smallestNetDistance = -1;
	int indexOfDBAnchorIsClosestTo = -1;
	//for loop figures out the index of which databaseAnchorsAsIntegers is the closest to the query point parameter provided.
	for (int i = 0 ; i < databaseAnchorRows; i++) {
		int latitudeDistance = abs(databaseAnchorsAsIntegers[i][0] - queryPointLatitude);
		int longitudeDistance = abs(databaseAnchorsAsIntegers[i][1] - queryPointLongitude);
		if(smallestNetDistance == -1) {
			smallestNetDistance = longitudeDistance + latitudeDistance;
			indexOfDBAnchorIsClosestTo = i;
		} else {
			if (longitudeDistance + latitudeDistance < smallestNetDistance) {
				smallestNetDistance = longitudeDistance + latitudeDistance;
				indexOfDBAnchorIsClosestTo = i;
			}
		}
	}
	
	//----------------------------------------------------------------------------------------------------------------------------
	//-----		2) Convert the query point to polar coordinates around the closest anchor point		------------------------------
	//----------------------------------------------------------------------------------------------------------------------------
	//Radius is found simply by calculating the distance of the query point lat long and the closest anchor point.
	int radiusOfQueryPointToAnchor = distanceEarth((double)databaseAnchorsAsIntegers[indexOfDBAnchorIsClosestTo][0], (double)databaseAnchorsAsIntegers[indexOfDBAnchorIsClosestTo][1], queryPointLatDouble, queryPointLongDouble);
	//s1 and s2 are the two sides of the triangle (third side is the radiusOfTheQueryPointToAnchor)
	int s1 = distanceEarth((double)databaseAnchorsAsIntegers[indexOfDBAnchorIsClosestTo][0], (double)databaseAnchorsAsIntegers[indexOfDBAnchorIsClosestTo][1], (double)databaseAnchorsAsIntegers[indexOfDBAnchorIsClosestTo][0], queryPointLongDouble);
	int s2 = distanceEarth(queryPointLatDouble, queryPointLongDouble, (double)databaseAnchorsAsIntegers[indexOfDBAnchorIsClosestTo][0], queryPointLongDouble);
	//Because distanceEarth gives an absolute value response, negative signs are added here if needed.
	if (queryPointLongitude < databaseAnchorsAsIntegers[indexOfDBAnchorIsClosestTo][1]) {
		s1 = s1 - (s1 * 2);
	}
	if (queryPointLongitude < databaseAnchorsAsIntegers[indexOfDBAnchorIsClosestTo][0]) {
		s2 = s2 - (s2 * 2);
	}
	//The angle of the polar coordinate is calculated below based on simple trig, before being converted into the right frame so as to represent the angle formed from the postive y axis going clockwise to the query point.
	double thetaOfQueryPointToAnchor = (atan(s2/s1) * 180 / PI);
	cout << "Angle of query point before it has been added to and rounded: " << thetaOfQueryPointToAnchor << "\n";
	//thetaOfQueryPointToAnchor is stored to 3 significant digits as this means angles can be represented accurately up to 10,000 miles away without loss of accuracy. Proven on paper via non right triangle math, law of cosines. Know two sides and one angle.
	thetaOfQueryPointToAnchor = round(thetaOfQueryPointToAnchor * 1000.0) / 1000.0;
	
	if (databaseAnchorsAsIntegers[indexOfDBAnchorIsClosestTo][0] < queryPointLatitude && queryPointLongitude < databaseAnchorsAsIntegers[indexOfDBAnchorIsClosestTo][1]) {
		cout << "Was in second quadrant.\n";
		thetaOfQueryPointToAnchor = 180.0 - ((thetaOfQueryPointToAnchor + 90.0) - 180.0);
	} else if (databaseAnchorsAsIntegers[indexOfDBAnchorIsClosestTo][0] > queryPointLatitude && queryPointLongitude < databaseAnchorsAsIntegers[indexOfDBAnchorIsClosestTo][1]) {
		cout << "Was in third quadrant.\n";
		thetaOfQueryPointToAnchor = 180.0 - ((thetaOfQueryPointToAnchor + 90.0) - 180.0);
	} else if(databaseAnchorsAsIntegers[indexOfDBAnchorIsClosestTo][0] > queryPointLatitude && queryPointLongitude > databaseAnchorsAsIntegers[indexOfDBAnchorIsClosestTo][1]) {
		cout << "Was in fourth quadrant.\n";
		thetaOfQueryPointToAnchor = 180.0 - ((thetaOfQueryPointToAnchor + 270.0) - 180.0);
	}
	
	cout << "Radius of query point: " << radiusOfQueryPointToAnchor << "\nAngle of query point: " << thetaOfQueryPointToAnchor << "\n";
	
	//----------------------------------------------------------------------------------------------------------------------------
	//-----		3) Query the database of spatial points in polar coordinates closest to the query point, searching for points roughly within the query radius from the query point		------------------------------
	//----------------------------------------------------------------------------------------------------------------------------
	int r = radiusOfQueryPointToAnchor;
	double t = thetaOfQueryPointToAnchor;
	int minR = r - qr;
	int maxr = r + qr;
	double minT = t - (atan(double(qr)/double(r)) * 180 / PI);
	double maxT = t + (atan(double(qr)/double(r)) * 180 / PI);
	
	vector<int> filteredRadii;
	vector<double> filteredAngles;
	vector<int> finalRadii;
	vector<double> finalAngles;
	vector<int> completeRadii;
	vector<double> completeAngles;
	
	//Accuracy
	vector<double> filteredLats;
	vector<double> filteredLongs;
	vector<double> finalLats;
	vector<double> finalLongs;
	vector<double> completeLats;
	vector<double> completeLongs;
	
	//ALTERNATIVE FIRST FILTER, assumes stored radius and angle values are sorted by radius.
	int minSavedRadius = arrayOfPoints[0][0];
	int midSavedRadius = arrayOfPoints[rows/2][0];
	int maxSavedRadius = arrayOfPoints[rows-1][0];
	double slope;
	int startIndex;
	if (r >= minSavedRadius && r <= midSavedRadius) {
		double slope = (double) (midSavedRadius - minSavedRadius) / ((double)rows/2.0);
		int startIndex = ((r-minSavedRadius) / slope);
	} else {
		double slope = (double) (maxSavedRadius - midSavedRadius) / ((double)rows/2.0);
		int startIndex = (rows/2) + ((r-midSavedRadius) / slope);
	}
	
	bool exitedRangeUpper = false;
	bool exitedRangeLower = false;
	
	int lowerBound = startIndex - 1;
	int upperBound = startIndex + 1;
	
	if (arrayOfPoints[startIndex][0] <= maxr && arrayOfPoints[startIndex][0] >= minR) {
		filteredRadii.push_back(arrayOfPoints[startIndex][0]);
		filteredAngles.push_back(arrayOfPointsAngles[startIndex][0]);
		//Accuracy
		filteredLats.push_back(arrayOfPointsLatLongs[startIndex][0]);
		filteredLongs.push_back(arrayOfPointsLatLongs[startIndex][1]);
		
		do {
			if (arrayOfPoints[lowerBound][0] <= maxr && arrayOfPoints[lowerBound][0] >= minR) {
				filteredRadii.push_back(arrayOfPoints[lowerBound][0]);
				filteredAngles.push_back(arrayOfPointsAngles[lowerBound][0]);
				//Accuracy
				filteredLats.push_back(arrayOfPointsLatLongs[lowerBound][0]);
				filteredLongs.push_back(arrayOfPointsLatLongs[lowerBound][1]);
			} else {
				exitedRangeLower = true;
			}
			
			if (arrayOfPoints[upperBound][0] <= maxr && arrayOfPoints[upperBound][0] >= minR) {
				filteredRadii.push_back(arrayOfPoints[upperBound][0]);
				filteredAngles.push_back(arrayOfPointsAngles[upperBound][0]);
				//Accuracy
				filteredLats.push_back(arrayOfPointsLatLongs[upperBound][0]);
				filteredLongs.push_back(arrayOfPointsLatLongs[upperBound][1]);
			} else {
				exitedRangeUpper = true;
			}
			
			upperBound++;
			lowerBound--;
			
		} while(exitedRangeUpper == false || exitedRangeLower == false);
	} else {
		if(arrayOfPoints[startIndex][0] > maxr) {
			do {
				startIndex = startIndex - 10;
			} while (arrayOfPoints[startIndex][0] > maxr);
			
			int lowerBound = startIndex - 1;
			int upperBound = startIndex + 1;
			
			do {
				if (arrayOfPoints[lowerBound][0] <= maxr && arrayOfPoints[lowerBound][0] >= minR) {
					filteredRadii.push_back(arrayOfPoints[lowerBound][0]);
					filteredAngles.push_back(arrayOfPointsAngles[lowerBound][0]);
					//Accuracy
					filteredLats.push_back(arrayOfPointsLatLongs[lowerBound][0]);
					filteredLongs.push_back(arrayOfPointsLatLongs[lowerBound][1]);
				} else {
					exitedRangeLower = true;
				}
				
				if (arrayOfPoints[upperBound][0] <= maxr && arrayOfPoints[upperBound][0] >= minR) {
					filteredRadii.push_back(arrayOfPoints[upperBound][0]);
					filteredAngles.push_back(arrayOfPointsAngles[upperBound][0]);
					//Accuracy
					filteredLats.push_back(arrayOfPointsLatLongs[upperBound][0]);
					filteredLongs.push_back(arrayOfPointsLatLongs[upperBound][1]);
				} else {
					exitedRangeUpper = true;
				}
				
				upperBound++;
				lowerBound--;
				
			} while(exitedRangeUpper == false || exitedRangeLower == false);
			
		} else {
			do {
				startIndex = startIndex + 10;
			} while (arrayOfPoints[startIndex][0] < minR);
			
			int lowerBound = startIndex - 1;
			int upperBound = startIndex + 1;
			
			do {
				if (arrayOfPoints[lowerBound][0] <= maxr && arrayOfPoints[lowerBound][0] >= minR) {
					filteredRadii.push_back(arrayOfPoints[lowerBound][0]);
					filteredAngles.push_back(arrayOfPointsAngles[lowerBound][0]);
					//Accuracy
					filteredLats.push_back(arrayOfPointsLatLongs[lowerBound][0]);
					filteredLongs.push_back(arrayOfPointsLatLongs[lowerBound][1]);
				} else {
					exitedRangeLower = true;
				}
				
				if (arrayOfPoints[upperBound][0] <= maxr && arrayOfPoints[upperBound][0] >= minR) {
					filteredRadii.push_back(arrayOfPoints[upperBound][0]);
					filteredAngles.push_back(arrayOfPointsAngles[upperBound][0]);
					//Accuracy
					filteredLats.push_back(arrayOfPointsLatLongs[upperBound][0]);
					filteredLongs.push_back(arrayOfPointsLatLongs[upperBound][1]);
				} else {
					exitedRangeUpper = true;
				}
				
				upperBound++;
				lowerBound--;
				
			} while(exitedRangeUpper == false || exitedRangeLower == false);
		}
	}
	
	for(int h = 0 ; h < filteredAngles.size() ; h++) {
		if(filteredAngles[h] <= maxT && filteredAngles[h] >= minT) {
			finalRadii.push_back(filteredRadii[h]);
			finalAngles.push_back(filteredAngles[h]);
			//Accuracy
			finalLats.push_back(filteredLats[h]);
			finalLongs.push_back(filteredLongs[h]);
		}
	}
	for ( int x = 0 ; x < finalRadii.size() ; x++) {
		if (finalAngles[x] < thetaOfQueryPointToAnchor) {
			double afterTan = tan(((finalAngles[x] - minT) * (PI / 180)));
			double part1 = sqrt(pow(qr,2) - pow((r*afterTan), 2)) ;
			double part2 = qr + minR;
			double highY = part1 + part2;//sqrt(pow(qr,2) - pow((r*tan(finalAngles[x] - minT)) * 180 / PI, 2)) + qr + minR;
			double lowY =  -part1 + part2;//-sqrt(pow(qr,2) - pow((r*tan(finalAngles[x] - minT)) * 180 / PI, 2)) + qr + minR;
			if (finalRadii[x] <= highY && finalRadii[x] >= lowY) {
				completeRadii.push_back(finalRadii[x]);
				completeAngles.push_back(finalAngles[x]);
				//Accuracy
				completeLats.push_back(finalLats[x]);
				completeLongs.push_back(finalLongs[x]);
			}
		} else {
			double afterTan = tan(((finalAngles[x] - thetaOfQueryPointToAnchor) * (PI / 180)));
			double part1 = sqrt(pow(qr,2) - pow((r*afterTan), 2)) ;
			double part2 = qr + minR;
			double highY = part1 + part2;//sqrt(pow(qr,2) - pow((r*tan(finalAngles[x] - minT)) * 180 / PI, 2)) + qr + minR;
			double lowY =  -part1 + part2;//-sqrt(pow(qr,2) - pow((r*tan(finalAngles[x] - minT)) * 180 / PI, 2)) + qr + minR;
			if (finalRadii[x] <= highY && finalRadii[x] >= lowY) {
				completeRadii.push_back(finalRadii[x]);
				completeAngles.push_back(finalAngles[x]);
				//Accuracy
				completeLats.push_back(finalLats[x]);
				completeLongs.push_back(finalLongs[x]);
			}
		}
	}
	auto end = chrono::high_resolution_clock::now();
	cout << fixed << "Scanned " << rows << " spatial points in " << chrono::duration_cast<chrono::nanoseconds>(end-begin).count() / 1000000000.0 << " seconds, finding " << completeRadii.size() <<  " points in the query radius." <<  "\n";
	
	return 0;
}


