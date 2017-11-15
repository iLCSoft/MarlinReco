#ifndef ChargeSpreader_h
#define ChargeSpreader_h

#include <marlin/Global.h>

#include <map>

struct AsicKey ;
class SimDigitalGeomCellId ;

struct ChargeSpreaderParameters
{
		float cellSize = 10.f ;
		float range = 30.f ;
		float padSeparation = 0.f ;

		//erf
		std::vector<float> erfWidth = {2} ;
		std::vector<float> erfWeigth = {1} ;

		//exact
		float d = 1.f ;
} ;


class ChargeSpreader
{
	public :
		ChargeSpreader() ;
		virtual ~ChargeSpreader() ;

		virtual void setParameters(ChargeSpreaderParameters param) { parameters = param ; }

		virtual void init() = 0 ;

		typedef std::pair<int,int> I_J_Coordinates ;
		virtual void addCharge( float charge , float posI , float posJ , SimDigitalGeomCellId* ) ;
		void newHit(float cellSize_) { chargeMap.clear() ; parameters.cellSize = cellSize_ ; }

		const std::map<I_J_Coordinates,float>& getChargeMap() const { return chargeMap ; }

	protected :
		virtual float computeIntegral(float x1 , float x2 , float y1 , float y2) const = 0 ;

		std::map<I_J_Coordinates,float> chargeMap ;
		ChargeSpreaderParameters parameters ;

		float normalisation = 0.f ;
} ;


class GaussianSpreader : public ChargeSpreader
{
	public :
		GaussianSpreader() ;
		virtual ~GaussianSpreader() ;
		virtual void init() ;

	protected :
		virtual float computeIntegral(float x1 , float x2 , float y1 , float y2) const ;
} ;

class ExactSpreader : public ChargeSpreader
{
	public :
		ExactSpreader() ;
		virtual ~ExactSpreader() ;
		virtual void init() ;

	protected :
		float computeIntegral(float x1 , float x2 , float y1 , float y2) const ;
} ;

class ExactSpreaderPerAsic : public ExactSpreader
{
	public :
		ExactSpreaderPerAsic(std::string fileName) ;
		virtual ~ExactSpreaderPerAsic() ;

		virtual void setParameters(ChargeSpreaderParameters param) { parameters = param ; dGlobal = parameters.d ; }

		virtual void addCharge(float charge, float posI, float posJ , SimDigitalGeomCellId* cellID) ;

	protected :

		float dGlobal = 1.0f ;
		void readFile(std::string fileName) ;


		std::map<AsicKey,float> dMap ;
} ;


#endif //ChargeSpreader_h
