#include "MathOperator.hh"
using std::vector;
namespace TTbarAnalysis 
{
	bool MathOperator::approximatelyEqual(const double * start1, const double * start2, double precision)
	{
		float distance = getDistance(start1, start2);
		return distance < precision;
	}
	bool MathOperator::approximatelyEqual(const double * start1, const float * start2, double precision)
	{
		//float distanceCut = 5.0;
		double end[3];
		for (int i = 0; i < 3; i++) 
		{
			end[i] = start2[i];
		}
		float distance = getDistance(start1, end);
		return distance < precision;
	}
	
	float MathOperator::getDistance(const double * start, const double * end)
	{
		vector<float> vec;
		for (int i = 0; i < 3; i++) 
		{
			vec.push_back(start[i] - end[i]);
		}
		float sqr = 0.0;
		for (int i = 0; i < 3; i++) 
		{
		        sqr += vec[i]*vec[i];
		}
		return sqrt(sqr);
	}
	float MathOperator::getDistance(const float * start, const float * end)
	{
		vector<float> vec;
		for (int i = 0; i < 3; i++) 
		{
			vec.push_back(start[i] - end[i]);
		}
		float sqr = 0.0;
		for (int i = 0; i < 3; i++) 
		{
		        sqr += vec[i]*vec[i];
		}
		return sqrt(sqr);
	}
	
	
	float MathOperator::getModule(const vector< int > & v)
	{
	       float module = 0.0;
	       for (unsigned int i = 0; i < v.size(); i++)
	       {
	               module += v[i]*v[i];
	       }
	       module = sqrt(module);
	       return module;
	}
	
	float MathOperator::getModule(const vector< float > & v)
	{
		float module = 0.0;
		for (unsigned int i = 0; i < v.size(); i++)
		{
		        module += v[i]*v[i];
		}
		module = sqrt(module);
		return module;
	}
	float MathOperator::getModule(const double * vector1)
	{
		float module = 0.0;
		for (int i = 0; i < 3; i++) 
		{
			module += vector1[i]*vector1[i];
		}
		return sqrt(module);
	}
	float MathOperator::getAngle(const double * vector1, const double * vector2)
	{
		float module1 = getModule(vector1);
		float module2 = getModule(vector2);
		if (module1 < 0.00000000001 || module2 < 0.000000000001) 
		{	
			return 10e10;
		}
		double product = 0.0;
		for (int i = 0; i < 3; i++) 
		{
			product += vector1[i]*vector2[i];
		}
		return acos((float)product/module1/module2);
	}
	vector< float > MathOperator::getAngles(vector< float > & direction)
	{
		vector< float > result;
		float epsilon = 0.00001;
		float semi = 1.5708;
		float pi = 2*semi;
		float phi = 0.0;
		if (direction[0] > 0.0 && direction[1] > 0.0 - epsilon) 
		{
			phi = atan(direction[1] / direction[0]); //(direction[0] < epsilon && direction[0] > 0.0 - epsilon)?
		}
		if (direction[0] < 0.0 && direction[1] > 0.0) 
		{
			phi = semi - atan(direction[1] / direction[0]) ;
		}
		if (direction[0] < 0.0 && direction[1] < 0.0 + epsilon) 
		{
			phi =  atan(direction[1] / direction[0]) + pi;
		}
		if (direction[0] > 0.0 && direction[1] < 0.0 - epsilon) 
		{
			phi = semi - atan(direction[1] / direction[0]) + pi;
		}
		if (direction[1] > 0.0 && direction[0] < 0.0 + epsilon && direction[0] > 0.0 -  epsilon) 
		{
			phi = semi;
		}
		if (direction[1] < 0.0 && direction[0] < 0.0 + epsilon && direction[0] > 0.0 -  epsilon) 
		{
			phi = pi + semi;
		}
		float teta = acos(direction[2]);
		result.push_back(phi);
		result.push_back(teta);
		return result;
	}
	vector< float > MathOperator::getDirection(vector< int > & vectorPoint1, vector< int > & vectorPoint2)
	{
		/*double * arr1 = MathOperator::castIntToDouble(&vectorPoint1[0]);
		double * arr2 = MathOperator::castIntToDouble(&vectorPoint2[0]);
		vector< float > result = MathOperator::getDirection(arr1,arr2);
		delete [] arr1;
		delete [] arr2;
		return result;*/
		vector< float > vector1;
		for (int i = 0; i < 3; i++)
		{
		        vector1.push_back((float)(vectorPoint1[i] - vectorPoint2[i]));
		}
		float module = getModule(vector1);
		for (int i = 0; i < 3; i++)
		{
		        vector1[i] = vector1[i]/module;
		}
		return vector1;
	}
	
	vector< float > MathOperator::getDirection(const double * vectorPoint1, const double * vectorPoint2)
	{
		vector< float > vector1;
		for (int i = 0; i < 3; i++)
		{
		        vector1.push_back(vectorPoint1[i] - vectorPoint2[i]);
		}
		float module = getModule(vector1);
		for (int i = 0; i < 3; i++)
		{
		        vector1[i] = vector1[i]/module;
		}
		return vector1;
	}
	vector< float > MathOperator::getDirection(const double * vectorPoint)
	{
		float module = getModule(vectorPoint);
		vector< float > vector1;
		for (int i = 0; i < 3; i++) 
		{
			vector1.push_back(vectorPoint[i] / module);
		}
		return vector1;
	}

	vector< float > * MathOperator::vectorProduct(const vector< float > & v1, const vector< float > & v2)
	{
		vector<float> * result = new vector<float>();
		result->push_back(v1[1]*v2[2]-v1[2]*v2[1]);
		result->push_back(v1[2]*v2[0]-v1[0]*v2[2]);
		result->push_back(v1[0]*v2[1]-v1[1]*v2[0]);
		return result;
	}

	float MathOperator::getDistanceTo( const vector< int > & vectorPoint1, const vector< float > & vector1,const vector< int > * point )
	{
		float result = 0.0;
		//vector< float > vector1 = getDirection(vectorPoint1,vectorPoint2);
		vector< float > vector2;
		for (int i = 0; i < 3; i++)
		{
		        vector2.push_back((float)(vectorPoint1[i] - point->at(i)));
		}
		vector< float > * product = vectorProduct(vector1,vector2);
		result = getModule(*product)/getModule(vector1);
		delete product;
		return result;
	}

	float MathOperator::getDistanceTo(const double * vectorPoint1, vector< float > & vector1,const double * point )
	{
		float result = 0.0;
		vector< float > vector2;
		for (int i = 0; i < 3; i++)
		{
		        vector2.push_back(vectorPoint1[i] - point[i]);
		}
		vector< float > * product = vectorProduct(vector1,vector2);
		result = getModule(*product)/getModule(vector1);
		delete product;
		return result;
	}

	vector< vector< int > * > * MathOperator::GetMagicNumbers()
	{
		vector< vector< int > * > * result = new vector< vector< int > * >();
		for (int x = -2; x < 3; x++)
		{
		        for (int y = -2; y < 3; y++)
			{
			        int z = 2;
				result->push_back(getPoint(x,y,z));
			}
			for (int y = -2; y < 3; y+=4)
			{
			        int z = 1;
				result->push_back(getPoint(x,y,z));
			}
			int y = -2;
			int z = 0;
			result->push_back(getPoint(x,y,z));

		}
		for (int x = -2; x < 3; x+=4)
		{
		        for (int y = -1; y < 2; y++)
			{
				int z = 1;
				result->push_back(getPoint(x,y,z));

			}
			int y = -1;
			int z = 0;
			result->push_back(getPoint(x,y,z));
		}
		result->push_back(getPoint(2,0,0));
		/*vector<int> additional = *getPoint(3,4,6);
		for (int i = 0; i < 3; i++) 
		{
			for (int j = -1; j < 2; j+=2) 
			{
				result->push_back(getPoint(additional[i],j,0));
				result->push_back(getPoint(additional[i],0,j));
				result->push_back(getPoint(0,additional[i],j));
				result->push_back(getPoint(j,additional[i],0));
				result->push_back(getPoint(j,0,additional[i]));
				result->push_back(getPoint(0,j,additional[i]));
			}
		}*/
		int z = 3;
		for (int y = -2; y < 3; y++) 
		{
			for (int x = -2; x < 3; x++) 
			{
				if (x == 0 && y == 0) 
				{
					continue;
				}
				result->push_back(getPoint(x,y,z));
			}
		}
		z = 4;
		for (int y = -1; y < 2; y++) 
		{
			for (int x = -1; x < 2; x++) 
			{
				if (x == 0 && y == 0) 
				{
					continue;
				}
				result->push_back(getPoint(x,y,z));
			}
		}
		return result;
	}
	vector< int > * MathOperator::getPoint(int x, int y, int z)
	{
		vector< int > * point = new vector< int >();
		point->push_back(x);
		point->push_back(y);
		point->push_back(z);
		return point;
	}
	float MathOperator::getPt(const double * momentum)
	{
		return sqrt(momentum[0] * momentum[0] + momentum[1] * momentum[1]);
	}
	float MathOperator::getRapidity(const double * momentum)
	{
		float p = getModule(momentum);
		float pL = momentum[2];
		float result = -1.0;
		if (p > pL) 
		{
			result = 0.5 * log((p+pL)/(p-pL));
		}
		return result;
	}
	double * MathOperator::getPtOnVector(const double * momentum, const float * target)
	{
		double * converted = toDoubleArray(target, 3);
		/*for (int i = 0; i < 3; i++) 
		{
			std::cout << i << ": " << converted[i];
		}
		std::cout << '\n';*/
		vector< float > direction = getDirection(converted);
		double * pt = new double[3];
		double product = 0.0;
		for (int i = 0; i < 3; i++) 
		{
			product += momentum[i] * direction[i];
		}
		for (int i = 0; i < 3; i++) 
		{
			pt[i] = momentum[i] - direction[i] * product;
		}
		//std::cout << "Pl: " << product << " |Pt|: " << getModule(pt) << '\n';
		return pt;
	}
	double MathOperator::getMissingPt(vector< const double * > & vectors, const float * target)
	{
		double sum[3];
		vector< double * > pts;
		for (unsigned int i = 0; i < vectors.size(); i++) 
		{
			pts.push_back(getPtOnVector(vectors[i], target));
		}
		for (int i = 0; i < 3; i++) 
		{
			for (unsigned int j = 0; j < vectors.size(); j++) 
			{
				sum[i] += pts[j][i]; 
			}
		}
		return getModule(sum);
	}
	double * MathOperator::toDoubleArray(const float * target, int size)
	{
		double * array = new double[3]();
		for (int i = 0; i < size; i++) 
		{
			array[i] = target[i];
		}
		return array;
	}
}
