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
	if (!checkRobust()) {
		return;
	}
	Point startingPoint;
	Point finalPoint;
	int size = controlPoints.size();
	float lastPointTime = controlPoints[size-1].time;
	float i;
	for(i=0; i < lastPointTime; i+=window){
		calculatePoint(startingPoint,i);
		if(i+window < lastPointTime){
			calculatePoint(finalPoint,i+window);
		}
		else{
			calculatePoint(finalPoint,lastPointTime);
		}
		DrawLib::drawLine(startingPoint, finalPoint, curveColor, curveThickness);
	}

	// Robustness: make sure there is at least two control point: start and end points
	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	// Note that you must draw the whole curve at each frame, that means connecting line segments between each two points on the curve

	return;
#endif
}
//Comp function sort curvePoint by distance
bool disComp(const CurvePoint &P1, const CurvePoint &P2) {
	
	return (P1.time<P2.time);
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
	newPosition.x = (controlPoints[sPoint].position.x * (2*normalTime*normalTime*normalTime-3*normalTime*normalTime+1))
						+(controlPoints[sPoint].tangent.x * (normalTime*normalTime*normalTime-2*normalTime*normalTime+normalTime) * intervalTime)
						+(controlPoints[ePoint].position.x * (-2*normalTime*normalTime*normalTime+3*normalTime*normalTime))
						+(controlPoints[ePoint].tangent.x * (normalTime*normalTime*normalTime-normalTime*normalTime) * intervalTime);

	newPosition.y = (controlPoints[sPoint].position.y * (2*normalTime*normalTime*normalTime-3*normalTime*normalTime+1))
						+(controlPoints[sPoint].tangent.y * (normalTime*normalTime*normalTime-2*normalTime*normalTime+normalTime)* intervalTime)
						+(controlPoints[ePoint].position.y * (-2*normalTime*normalTime*normalTime+3*normalTime*normalTime))
						+(controlPoints[ePoint].tangent.y * (normalTime*normalTime*normalTime-normalTime*normalTime)* intervalTime);

	newPosition.z = (controlPoints[sPoint].position.z * (2*normalTime*normalTime*normalTime-3*normalTime*normalTime+1))
						+(controlPoints[sPoint].tangent.z * (normalTime*normalTime*normalTime-2*normalTime*normalTime+normalTime)* intervalTime)
						+(controlPoints[ePoint].position.z * (-2*normalTime*normalTime*normalTime+3*normalTime*normalTime))
						+(controlPoints[ePoint].tangent.z * (normalTime*normalTime*normalTime-normalTime*normalTime)* intervalTime);
	// Return result
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time) {
	Point newPosition;
	//1.calculate Tangent line (M0, M1)
	Point p0=controlPoints[nextPoint-1].position;
	Point p1=controlPoints[nextPoint].position;


	Point m0;
	if (nextPoint-1 == 0) {//if p0 is first point
		Point current = controlPoints[nextPoint - 1].position;
		Point after = controlPoints[nextPoint - 1 + 1].position;

		m0.x = after.x - current.x;
		m0.y = after.y - current.y;
		m0.z = after.z - current.z;
	}
	else if (nextPoint - 1 == controlPoints.size() - 1) {//if p0 is last point
		Point current = controlPoints[nextPoint - 1].position;
		Point before = controlPoints[nextPoint - 1 - 1].position;

		m0.x = current.x - before.x;
		m0.y = current.y - before.y;
		m0.z = current.z - before.z;
	}
	else {//if p0 is point in the middle
		Point before = controlPoints[nextPoint - 1 - 1].position;
		Point after = controlPoints[nextPoint - 1 + 1].position;
		Point current = controlPoints[nextPoint - 1].position;

		m0.x = 0.5f*(after.x - before.x);
		m0.y = 0.5f*(after.y - before.y);
		m0.z = 0.5f*(after.z - before.z);
	}
	Point m1;
	if (nextPoint == 0) {//if p1 is first point
		Point current = controlPoints[nextPoint].position;
		Point after = controlPoints[nextPoint + 1].position;

		m1.x = after.x - current.x;
		m1.y = after.y - current.y;
		m1.z = after.z - current.z;
	}
	else if (nextPoint == controlPoints.size() - 1) {//if p1 is last point
		Point current = controlPoints[nextPoint].position;
		Point before = controlPoints[nextPoint - 1].position;

		m1.x = current.x - before.x;
		m1.y = current.y - before.y;
		m1.z = current.z - before.z;
	}
	else {//if p1 is point in the middle
		Point before = controlPoints[nextPoint - 1].position;
		Point after = controlPoints[nextPoint + 1].position;
		Point current = controlPoints[nextPoint].position;

		m1.x = 0.5f*(after.x - before.x);
		m1.y = 0.5f*(after.y - before.y);
		m1.z = 0.5f*(after.z - before.z);
	}


	//2.calculate equation H(t)
	float TimeInterval=controlPoints[nextPoint].time-controlPoints[nextPoint-1].time;
	float t=(time-controlPoints[nextPoint-1].time)/TimeInterval;

	newPosition.x=(2.0f * t * t * t - 3.0f * t * t + 1.0f) * p0.x 
					+ (t * t * t - 2.0f * t * t + t) * m0.x 
					+ (-2.0f * t * t * t + 3.0f * t * t) * p1.x 
					+ (t * t * t - t * t) * m1.x;
	newPosition.y=(2.0f * t * t * t - 3.0f * t * t + 1.0f) * p0.y 
					+ (t * t * t - 2.0f * t * t + t) * m0.y 
					+ (-2.0f * t * t * t + 3.0f * t * t) * p1.y 
					+ (t * t * t - t * t) * m1.y;
	newPosition.z=(2.0f * t * t * t - 3.0f * t * t + 1.0f) * p0.z
					+ (t * t * t - 2.0f * t * t + t) * m0.z 
					+ (-2.0f * t * t * t + 3.0f * t * t) * p1.z 
					+ (t * t * t - t * t) * m1.z;
	// Return result
	return newPosition;
}