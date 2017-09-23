//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
// Copyright (c) 2015 Mahyar Khayatkhoei
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"

using namespace Util;

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType) {
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType) {
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint) {
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints) {
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window) {
#ifdef ENABLE_GUI

	//================DELETE THIS PART AND THEN START CODING===================
	static bool flag = false;
	if (!flag) {
		std::cerr << "ERROR>>>>Member function drawCurve is not implemented!" << std::endl;
		flag = true;
	}
	//=========================================================================

	// Robustness: make sure there is at least two control point: start and end points
	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	// Note that you must draw the whole curve at each frame, that means connecting line segments between each two points on the curve

	return;
#endif
}
//Comp function sort curvePoint by distance
bool disComp(const CurvePoint &P1, const CurvePoint &P2) {
	float d1 = sqrt(pow(P1.position.x, 2) + pow(P1.position.y, 2) + pow(P1.position.z, 2));
	float d2 = sqrt(pow(P2.position.x, 2) + pow(P2.position.y, 2) + pow(P2.position.z, 2));
	return (d1<d2);
}
// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints() {

	std::sort(controlPoints.begin(), controlPoints.end(), disComp);
	return;
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
// Note that this function should return false if the end of the curve is reached, or no next point can be found
bool Curve::calculatePoint(Point& outputPoint, float time) {
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;

	// Define temporary parameters for calculation
	unsigned int nextPoint;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	// Note that nextPoint is an integer containing the index of the next control point
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve given the next control point (nextPoint)
	if (type == hermiteCurve) {
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve) {
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Roboustness
bool Curve::checkRobust() {
	if (type == hermiteCurve && controlPoints.size() >= 2) {

		return true;
	}

	else if (type == catmullCurve && controlPoints.size() >= 3) {

		return true;
	}

	return false;
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time) {
	for (int i = 0; i < controlPoints.size(); i++) {
		if (controlPoints[i].time > time) {
			nextPoint = i;
			return true;
		}
	}
	return false;
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;

	int sPoint = nextPoint-1;		//starting point
	int ePoint = nextPoint;			//end point

	intervalTime = controlPoints[ePoint].time - controlPoints[sPoint].time;
	normalTime = (time - controlPoints[sPoint].time) / intervalTime;

	// Calculate position at t = time on Hermite curve
	float newPositionX = (controlPoints[sPoint].position.x * (2*normalTime*normalTime*normalTime-3*normalTime*normalTime+1))
						+(controlPoints[sPoint].tangent.x * (normalTime*normalTime*normalTime-2*normalTime*normalTime+normalTime))
						+(controlPoints[ePoint].position.x * (-2*normalTime*normalTime*normalTime+3*normalTime*normalTime))
						+(controlPoints[ePoint].tangent.x * (normalTime*normalTime*normalTime-normalTime*normalTime));

	float newPositionY = (controlPoints[sPoint].position.y * (2*normalTime*normalTime*normalTime-3*normalTime*normalTime+1))
						+(controlPoints[sPoint].tangent.y * (normalTime*normalTime*normalTime-2*normalTime*normalTime+normalTime))
						+(controlPoints[ePoint].position.y * (-2*normalTime*normalTime*normalTime+3*normalTime*normalTime))
						+(controlPoints[ePoint].tangent.y * (normalTime*normalTime*normalTime-normalTime*normalTime));

	float newPositionZ = (controlPoints[sPoint].position.z * (2*normalTime*normalTime*normalTime-3*normalTime*normalTime+1))
						+(controlPoints[sPoint].tangent.z * (normalTime*normalTime*normalTime-2*normalTime*normalTime+normalTime))
						+(controlPoints[ePoint].position.z * (-2*normalTime*normalTime*normalTime+3*normalTime*normalTime))
						+(controlPoints[ePoint].tangent.z * (normalTime*normalTime*normalTime-normalTime*normalTime));
	// Return result
	newPosition = Point(newPositionX,newPositionY,newPositionZ);
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time) {
	Point newPosition;

	//================DELETE THIS PART AND THEN START CODING===================
	static bool flag = false;
	if (!flag) {
		std::cerr << "ERROR>>>>Member function useCatmullCurve is not implemented!" << std::endl;
		flag = true;
	}
	//=========================================================================

	// Calculate position at t = time on Catmull-Rom curve

	// Return result
	return newPosition;
}