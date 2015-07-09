#ifndef _MathOperator_hh_
#define _MathOperator_hh_
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <iostream>
namespace TTbarAnalysis 
{
	class MathOperator 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			MathOperator ();
			virtual ~MathOperator ();
		//
		//	Methods
		//
			static float getModule(const std::vector< int > & v);
			static float getModule(const std::vector< float > & v);
			static float getModule(const double * vector1);
			static float getModule(const float * vector1);
			
			static std::vector< float > getAngles(std::vector< float > & direction);
			static float getAngle(const double * vector1, const double * vector2);
			static float getAngle(const float * vector1, const float * vector2);

			static bool approximatelyEqual(const double * start1, const double * end);
			static float getDistance(const double * start, const double * end);
			static float getDistance(const float * start, const float * end);
			static std::vector< float > * vectorProduct(const std::vector< float > & v1,const std::vector< float > & v2);
			
			static float getDistanceTo(const std::vector< int > & vectorPoint1,const std::vector< float > & vector1, const std::vector< int > * pointOfLine );
			static float getDistanceTo(const double * vectorPoint1, std::vector< float > & vector1,const double * pointOfLine );
			
			static float getDistanceBtw(const double * vectorPoint1, std::vector< float > & vectorDirection1,const double * vectorPoint2, std::vector< float > & vectorDirection2 );
			
			static std::vector< float > getDirection(std::vector< int > & vectorPoint1, std::vector< int > & vectorPoint2);
			static std::vector< float > getDirection(const double * vectorPoint1, const double * vectorPoint2);
			static std::vector< float > getDirection(const double * vectorPoint);
			
			static float getTrace(const std::vector< float > & matrix, int dimension);
			
			static std::vector< std::vector< int > * > * GetMagicNumbers();
			static std::vector< int > * getPoint(int x, int y, int z);
			static double * getPtOnVector(const double * momentum, const float * target);
			static double getMissingPt(std::vector< const double * > & vectors, const float * target);
			static double * toDoubleArray(const float * target, int size);
		private:
		//
		//	Data
		//
			/* data */
		//
		//	Private methods
		//
			static double * castIntToDouble(int * array);
	};
}
#endif
