#ifndef ChargeInducer_h
#define ChargeInducer_h

#include <string>
#include <map>

#include <random>

class SimDigitalGeomCellId ;
struct AsicKey ;

class ChargeInducer
{
	public :
		ChargeInducer() ;
		virtual ~ChargeInducer() ;
		virtual float getCharge(SimDigitalGeomCellId* cellID) = 0 ;

		void setSeed(unsigned int value) ;

	protected :
		std::mt19937 generator ;
} ;

class UniformPolya : public ChargeInducer
{
	public :
		UniformPolya(float _qbar , float _theta) ;
		~UniformPolya() ;

		virtual float getCharge(SimDigitalGeomCellId* cellID) ;


	protected :
		std::gamma_distribution<float> gammadist ;

} ;

class AsicPolya : public UniformPolya
{
	public :
		AsicPolya(float _qbar , float _theta , std::string fileName) ;
		~AsicPolya() ;

		virtual float getCharge(SimDigitalGeomCellId* cellID) ;


	protected :
		void readFile(std::string fileName) ;

		std::map<AsicKey , std::gamma_distribution<float> > polyaMap ;
} ;

#endif //ChargeInducer_h
