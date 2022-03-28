// ////////////////////////////////////////////////////////////////////// //
//                                                                        //
// SiStripDigi - Marlin Processor - digitizing data from Si strip sensors //
//                                                                        //
// ////////////////////////////////////////////////////////////////////// //

/**
\addtogroup SiStripDigi SiStripDigi
@{

 *\mainpage notitle
 *
 *<br><br><h1>Digitization of Si Microstrip Detectors - Marlin Processors' Documentation</h1>
 *
 *<h4 style="text-align: center"> version 01 </h4>
 *
 *<h3>Description of SiStripDigi processor:</h3> SiStripDigi (a Marlin processor)
 * represents a detailed strip detector digitizer that's mainly intended for detailed
 * simulation studies; but in principle, can be used for full physics studies. The
 * geometry information required for data processing is being accessed on the fly via
 * a geometry interface: SiStripGeom that simply utilizes the information obtained
 * from given GEAR *.xml file. As regards the solid-state physics, the processor
 * naturally implements the following processes for e-h pairs created by a
 * high-energy particle: drift in electric field (simplified just as electric
 * field of a p-n junction - oriented in negative direction of X axis of local
 * reference system, i.e. electrons drift in positive direction X, holes vice
 * versa), diffusion due to multiple collisions and Lorentz shift in magnetic
 * field; furthermore, it also takes into account electronics processes:
 * mutual micro-strip crosstalks (dependent on AC/DC coupling), electronics noise, the
 * type of current collected by strips and set if there are floating string in between
 * real read-out strips. As DSSDs are assumed to be digitized, both electrons and
 * holes are collected. Electrons drift to strips oriented in R-Phi direction, holes
 * to strips oriented in Z direction. A user can naturally influence these effects
 * when changing processor parameters: sensor depletion voltage (set by default to 60 V),
 * sensor bias voltage (set by default to 150 V), temperature of a silicon wafer (set by
 * default to 300 K), interstrip capacitance (set by default to 6 pF), strip-to-backplane
 * capacitance (set by default to 0 pF), coupling capacitance (set by default to 120 pF)
 * and sigma of CMS-like noise (set by default to 1200 electrons). Finally, to achieve
 * higher effectivity one can influence the simulation precision using these parameters:
 * space precision (set by default to 20 um), relative Lorentz angle precision (set by
 * default to 1%) and relative drift time precision (set by default to 1%). As the space
 * precision might be different from the one set in Geant 4, the Landau fluctuations
 * of deposited energy are simulated here too, using internal fluctuator SiEnergyFluct,
 * which exactly follows a Geant4 class called G4UniversalFluctuation. If the particle is
 * regarded as low energy (below MIP beta*gamma factor), the energy is distributed
 * uniformly using the info from Geant 4, otherwise internal fluctuation is used. As an
 * output one obtains a standard LCIO collection of TrackerPulses, so to obtain hits,
 * one still must perform clustering, for instance SiStripClus processor can do it ...
 *
 *<h3>Description of SiStripClus processor:</h3>SiStripClus (a Marlin processor)
 * provides a cluster finding algorithm based on COG (centre of gravity) method (cluster
 * size is lower than 3) or on head-tail analog method (cluster-size is equal or higher than
 * 3) and thus transforms electric pulses (TrackerPulses - obtained from SiStripDigi processor)
 * into real 2D hits (RPhi + Z) (resp. 3D hits, where the 3rd component is calculated
 * as a center of a sensor). First, a seed strip in Z direction is found (seed strip is a
 * strip with signal higher than roughly 5 * noise level - set by user), than the processor
 * looks for neighbouring strips (with signals higher than roughly 2 * noise level - set by
 * user) and thus finds a 1D cluster in Z. The same is performed in R-Phi direction (pitch in
 * R-Phi might differ in forward region based on current Z position) and final 2D hit, resp.
 * 3D hit is calculated (including "blind" hits, so-called ghosts). All the computation is
 * performed in local reference system (using information from the geometry interface:
 * SiStripGeom) and final 3D hits are transformed back into global ref. system. In
 * order to obtain quasi-realistic hit covariance matrices, the resolution can be simulated
 * using this processor (use #define ROOT_OUTPUT) and then saved as input Marlin resolution
 * parameters ...
 *
 *
 * If you have any questions or comments, please write an email to the author:
 * <a class="email" href="mailto:drasal@ipnp.troja.mff.cuni.cz"> Zbynek Drasal</a>,
 * Charles University, Prague.
 *
 * CLHEP   Namespace of Class Library in High Energy Physics
 * gear    Namespace of GEometry Api for Reconstruction
 * lcio    Namespace of Linear Collider InputOutput
 * marlin  Namespace of Modular Analysis and Reconstruction tool for LINear collider
 * std     Standard library namespace
 * sistrip SiStripDigi processor namespace
 */

#ifndef SISTRIPDIGI_H
#define SISTRIPDIGI_H 1

// Define ROOT output if needed
//#define ROOT_OUTPUT_LAND

#include <vector>
#include <queue>
#include <map>
#include<string>

// Include CLHEP header files
#include <CLHEP/Random/RandGauss.h>
#include <CLHEP/Vector/ThreeVector.h>

// Include Digi header files
#include "DigiCluster.h"
#include "SiEnergyFluct.h"
#include "Signal.h"
#include "SimTrackerDigiHit.h"
#include "SiStripGeom.h"

// Include LCIO header files
#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <IMPL/SimTrackerHitImpl.h>

// Include Marlin
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <streamlog/streamlog.h>

// Include ROOT classes
#ifdef ROOT_OUTPUT_LAND
#include <TFile.h>
#include <TTree.h>
#endif

class RombIntSolver;

namespace sistrip {

// FIXME: Use the StripType enumerate
#define STRIPFRONT RPhi   // Const denoting strips in sensors 1,2 (front)
#define STRIPREAR    Z    // Const denoting strips in sensors 3,4 (rear)
//FIXME: TO DELETE: Solo para compilar
#define STRIPRPHI RPhi   // Const denoting strips in sensors 1,2 (front)
#define STRIPZ    Z    // Const denoting strips in sensors 3,4 (rear)

#define ROUNDEPS    5.   // Avoid rounding errors - space precision for round. errors in um
#define SEED    12586    // Random generator initialization

// Typedefs
typedef const std::vector< std::string >           ConstStringVec;
typedef       std::vector< std::string >           StringVec;
typedef       std::queue < std::string >           StringQue;
typedef       std::vector< DigiCluster *>          DigiClusterVec;
typedef       std::vector< EVENT::LCCollection *>  LCCollectionVec;
typedef const std::vector< EVENT::SimTrackerHit *> ConstSimTrackerHitVec;
typedef       std::vector< EVENT::SimTrackerHit *> SimTrackerHitVec;
typedef       std::vector< SimTrackerDigiHit *>    SimTrackerDigiHitVec;
typedef       std::map< int, Signal *>             StripChargeMap;   // strip ID, signal
typedef       std::map< int, StripChargeMap *>     SensorStripMap;   // sensor ID, strip map
typedef       std::map<SimTrackerHit *, float>     SimHitMap;        // hit, weight

/**
\addtogroup SiStripDigiProcessor SiStripDigiProcessor
\ingroup SiStripDigi
@{
Marlin processor intended for detailed digitization of silicon strip sensors.
//!
//! @author Z. Drasal, Charles University, Prague
//!
*/
class SiStripDigi : public marlin::Processor 
{
	public:
		//!Method that returns a new instance of this processor
		virtual Processor * newProcessor() { return new SiStripDigi(); }
		
		//!Constructor - set processor description and register processor parameters
		SiStripDigi();
		
		//!Method called at the beginning of data processing - used for initialization
		virtual void init();
		
		//!Method called for each run - used for run header processing
		virtual void processRunHeader(LCRunHeader * run);
		
		//!Method called for each event - used for event data processing
		virtual void processEvent(LCEvent * event);
		
		//!Method called after each event - used for data checking
		virtual void check(LCEvent * event);
		
		//!Method called after all data processing
		virtual void end();
		
	protected:
		
		// MAIN DIGI METHOD
		
		//!Method digitizing given SimTrackerDigiHit - takes into account all relevant
		//!physical processes: landau fluctuations, drift, diffusion, Lorentz shift 
		//!in magnetic field (input parameter: a pointer to digitized hit, output 
                //!parameter: a sensor map of strips with total integrated charge and time 
		//!when particle crossed the sensor)
		void digitize(const SimTrackerDigiHit * hit, SensorStripMap & sensorMap);
		
		
		// OTHER METHODS
		
		//!Method that calculates e-h clusters, where clusters are worked out from
		//!the hits (simulated in Geant 4 and transformed into local reference
		//!system) using predefined space precision. As the precision might be
		//!different from a precision in Geant 4, the Landau fluctuations are
		//!performed here too, using internal Landau fluctuator. Required parameters:
		//!a pointer to SimTrackerDigiHit and two vectors of pointers to electron,
		//!resp. hole clusters (pairs) )
		void calcClusters(const SimTrackerDigiHit * hit, DigiClusterVec & hClusterVec);
		//DigiClusterVec & eClusterVec, DigiClusterVec & hClusterVec);
		
		//!Method that calculates crosstalk effect, i.e. each total charge is 
		//!redistributed according the following relation: 
		//!          Q_neigh = Q_centr * C_inter/(C_inter + C_back + C_coupl),
		//!where neigh denotes neighbouring, centr central, inter interstrip,
		//!back strip2backplane and coupl coupling.
		void calcCrossTalk(SensorStripMap & sensorMap);
		
		//!Method generating random noise using Gaussian distribution and add this 
		//!effect to the final results.
		void genNoise(SensorStripMap & sensorMap);
		
		//!Method transforming given SimTrackerHit into local ref. system of each 
		//!sensor, resp. wafer, where the center is positioned such as x, y and z 
		//!coordinates are always positive and x is in direction of thickness.
		void transformSimHit(SimTrackerDigiHit * simDigiHit);
		
		// GET METHODS
		
		//!Get method - returns ELECTRON diffusivity (parameters: electron mobility
		//!in cm^2/s, respectively position in cm). The diffusivity is calculated
		//!according to the Einstein relation:
		//!
		//!  &nbsp;&nbsp; D = k*temp/e * mobil(E(x),temp) in cm^2/s
		//!
		//!Where k represents Boltzmann constant and e an elementary charge.
		double getElecDiffusivity(double pos);
		
		//!Get method - returns HOLE diffusivity (parameters: hole mobility in
		//!cm^2/s, respectively position in cm). The diffusivity is calculated
		//!according to the Einstein relation:
		//!
		//!  &nbsp;&nbsp; D = k*temp/e * mobil(E(x),temp) in cm^2/s
		//!
		//!Where k represents Boltzmann constant and e an elementary charge.
		double getHoleDiffusivity(double pos);
		
		//!Get method - returns total ELECTRON drift time (parameters: electron 
		//!position in cm). The total time is obtained as a numerical solution of
		//!following equation of motion:
		//!
		//!  &nbsp;&nbsp; dx/dt = v = mobil(E(x), temp) * E(x)
		//!  &nbsp;&nbsp; Int(x0,d){1/(mobility(E(x), temp) * E(x)) * dx} = t
		//!  &nbsp;&nbsp; Int(x0,d){1/v(E(x), temp) * dx} = t
		//!
		//!Where mobil represents electron mobility, d sensor thickness, which
		//!corresponds to sensor back side (DSSDs - back strip side), i.e. the
		//!one with high voltage set; and v to electron actual velocity.
		double getElecDriftTime(double pos);
		
		//!Get method - returns total HOLE drift time (parameters: hole position
		//!in cm). The total time is obtained as a numerical solution of following
		//!equation of motion:
		//!
		//!  &nbsp;&nbsp; dx/dt = v = mobil(E(x), temp) * E(x)
		//!  &nbsp;&nbsp; Int(x0,0){1/(mobility(E(x), temp) * E(x))*dx} = t
		//!  &nbsp;&nbsp; Int(x0,0){1/v(E(x), temp) * dx} = t
		//!
		//!Where mobil represents hole mobility, 0 corresponds to sensor strip side
		//!(with high voltage set to 0) and v to hole actual velocity.
		double getHoleDriftTime(double pos);
		
		//!Get method - returns electric intensity (parameters: position X in cm).
		//!The intensity is calculated using analytically expressible relation
		//!between electron, resp. hole, position and electric field, derived for
		//!simple abrupt p-n junction with depletion voltage Vd and bias voltage Vb.
		//!
		//! Relation:
		//!
		//!  &nbsp;&nbsp; E(x) = -(Vb+Vd/d - 2*x/d^2*Vd) in V/cm
		//!
		//!Where d represents sensor thickness.
		double getEField(double pos);
		
		//!Get method - returns ELECTRON mobility (parameters: electric intenzity in
		//!V/cm, respectively position X in cm). For the region where the mobility
		//!is dependent on the electric field, applied along <111> direction,
		//!following parametrization can be used:
		//!
		//!  &nbsp;&nbsp; mobil = vm/Ec * 1/(1 + (E(x)/Ec)^beta)^(1/beta) in cm^2/V.s
		//!
		//!  &nbsp;&nbsp; vm = 1.53*10^9*temp^(-0.87) cm/s
		//!  &nbsp;&nbsp; Ec = 1.01*temp^(1.55) V/cm
		//!  &nbsp;&nbsp; beta = 2.57*10^(-2)*temp(0.66)
		//!
		//!For more details about parametrization see: A review of some charge
		//!transport properties of silicon by C.Jacobini et. al. (Solid-State
		//!Electronics, 1977, Vol. 20, p. 77-89)
		double getElecMobility(double pos);
		
		//!Get method - returns HOLE mobility (parameters: electric intenzity in
		//!V/cm, respectively position X in cm). For the region where the mobility
		//!is dependent on the electric field, applied along <111> direction,
		//!following parametrization can be used:
		//!
		//!  &nbsp;&nbsp; mobil = vm/Ec * 1/(1 + (E(x)/Ec)^beta)^(1/beta) in cm^2/V.s
		//!
		//!  &nbsp;&nbsp; vm = 1.62*10^8*temp^(-0.52) cm/s
		//!  &nbsp;&nbsp; Ec = 1.24*temp^(1.68) V/cm
		//!  &nbsp;&nbsp; beta = 0.46*temp(0.17)
		//!
		//!For more details about parametrization see: A review of some charge
		//!transport properties of silicon by C.Jacobini et. al. (Solid-State
		//!Electronics, 1977, Vol. 20, p. 77-89)
		double getHoleMobility(double pos);
		
		//!Method used for precise calculation of tan(Lorentz angle) - an angle of
		//!ELECTRON deflection due to magnetic field. (parameters: position z in
		//!cm). The Lorentz angle calculation utilizes so-called Romberg integration
		//!method, quite fast and reasonably precise method. The precision of the
		//!calculation is set via variable: _epsAngle by a user.
		//!
		//! Relation:
		//!
		//!  &nbsp;&nbsp; tan(thetaL) =  Int(x,d){mobil(E(z))*r*B*dz}/Int(z,d){dz}
		//!  &nbsp;&nbsp; shiftX      = -Int(x,d){mobil(E(z))*r*By*dz}
		//!  &nbsp;&nbsp; shiftY      =  Int(x,d){mobil(E(z))*r*Bx*dz}
		//!
		//!Where r = 1.13 + 0.0008*(temp-273) represents Hall scattering factor
		//!for electrons and d sensor thickness, corresponding to sensor back side.
		CLHEP::Hep3Vector getElecLorentzShift(double posZ);
		
		//!Method used for precise calculation of tan(Lorentz angle) - an angle of
		//!HOLES deflection due to magnetic field. (parameters: position z in
		//!cm). The Lorentz angle calculation utilizes so-called Romberg integration
		//!method, quite fast and reasonably precise method. The precision of
		//!the calculation is set via variable: _epsAngle by a user.
		//!
		//! Relation:
		//!
		//!  &nbsp;&nbsp; tan(thetaL) = Int(z,0){mobil(E(z))*r*B*dz}/Int(z,0){dz}
		//!
		//!
		//!Where r = 0.72 - 0.0005*(temp-273) represents Hall scattering factor
		//!for holes and 0 corresponds to sensor strip side.
		CLHEP::Hep3Vector getHoleLorentzShift(double pos);

                //!Get method - returns ELECTRON actual velocity, calculated as follows:
                //!
                //!  &nbsp;&nbsp; dz/dt = v = mobil(E(x), temp) * E(x)
                //!
                //!Where E corresponds to elec. intenzity a mobil to electron mobility.
                //!(parameters: electron position in cm)
		double getElecVelocity(double pos);
                
                //!Get method - returns actual ELECTRON inverse velocity, calculated as 
		//!follows:
                //!
                //!  &nbsp;&nbsp; dt/dx = v^-1 = 1/(mobil(E(x), temp) * E(x))
                //!
                //!Where E corresponds to elec. intenzity a mobil to electron mobility.
                //!(parameters: electron position in cm)
		double getElecInvVelocity(double pos);
                
                //!Get method - returns HOLE actual velocity, calculated as follows:
                //!
                //!  &nbsp;&nbsp; dx/dt = v = mobil(E(x), temp) * E(x)
                //!
                //!Where E corresponds to elec. intenzity a mobil to hole mobility.
                //!(parameters: hole position in cm)
                double getHoleVelocity(double pos);
                
                //!Get method - returns actual HOLE inverse velocity, calculated as follows:
                //!
                //!  &nbsp;&nbsp; dt/dx = v^-1 = 1/(mobil(E(x), temp) * E(x))
                //!
                //!Where E corresponds to elec. intenzity a mobil to hole mobility.
                //!(parameters: hole position in cm)
                double getHoleInvVelocity(double pos);
                
                
                // PRINT METHODS
                
                //!Method printing cluster info
                void printClusterInfo( const DigiCluster & cluster) const;
                
                //!Method printing clusters info (parameter: info - extra information)
                void printClustersInfo( std::string info, 
				const DigiClusterVec & clusterVec) const;
                
                //!Method printing hit info (parameter: info - extra information)
                void printHitInfo( std::string info, const SimTrackerDigiHit * hit) const;
                
                //!Method printing processor parameters
                void printProcessorParams() const;
                
                //!Method printing info about signals at each strip
                void printStripsInfo( std::string info, 
				const SensorStripMap & sensorMap) const;
                
                
                // VARIABLES
		
		// Collection names
		std::string _inColName;         //!< LCIO input collection name
		std::string _outColName;        //!< LCIO output collection name
		std::string _relColNamePlsToSim;//!< LCIO relation collection name -
		                                //!< TrkPulse (Digit)  <-> SimTrkHit

		// Digitization parameters - set by users   
		std::string _subdetector;
		
		float _Vdepl;                   //!< Sensor depletion voltage in volts
		float _Vbias;                   //!< Sensor bias voltage in volts
		float _temp;                    //!< Sensor temperature in Kelvins
		
		bool  _electronicEffects;       //!< Add crosstalk + noise?
		float _capInterStrip;           //!< Mutual interstrip capacitance
		float _capBackPlane;            //!< Strip-to-backplane capacitance
		float _capCoupl;                //!< AC coupling - capacitance
		float _elNoise;                 //!< CMS-like (common mode subtracted) noise added to the signal
		bool  _floatStripsRPhi;         //!< Is every even strip floating in R-Phi?
		bool  _floatStripsZ;            //!< Is every even strip floating in Z?
		
		float _epsSpace;                //!< Absolute digi precision in space in um
		float _epsAngle;                //!< Relative digi precision in Lorentz angle
		float _epsTime;                 //!< Relative digi precision in Drift time
		
		// Digitization parameters for Landau fluctuations - set by users
		bool   _landauFluct;            //!< Define if internal Landau fluctuations (instead of Geant4) used
		float  _landauBetaGammaCut;     //!< Below this beta*gamma factor internal Landau fluctuations not used
		double _prodThreshOnDeltaRays;  //!< Production threshold cut on delta electrons in keV (for Landau fluct.) - use the same as in Geant4 (80keV ~ 0.05 mm)

		double _pitchFront;            //!< Pitch in the front sensors (in the middle)
		double _pitchRear;             //!< Pitch in the rear sensors (in the middle)

		// Digitization parameters for background studies - set by users
		bool   _integrationWindow;      //!< Use integration window?
		double _startIntegration;       //!< Start time of integration of the sensors in ns (everything before this value will not be digitized)
		double _stopIntegration;        //!< Stop time of integration of the sensors in ns(everything after this value will not be digitized)
		
		// Digitization parameters - obtained from Gear xml file
		int       _currentCellID;       //!< Actual layer+ladder+sensor ID encoded into one number
		short int _currentLayerID;      //!< Actual layer ID
		short int _currentLadderID;     //!< Actual ladder ID
		short int _currentSensorID;     //!< Actual sensor ID
		float     _sensorThick;         //!< Actual sensor - Si wafer thickness in system of units
		
		// Magnetic field - obtained from Gear xml file and transformed into 
		//                  local system
		CLHEP::Hep3Vector _magField;    //!< Magnetic field in T in detector reference system
		
		// Geometry parameters
		SiStripGeom * _geometry;        //!< All geometry information from Gear xml file
		
		// Random generator
		CLHEP::RandGauss * _genGauss;   //!< Random number generator - Gaussian distribution
		
		// Simulator of Landau fluctuations in Si
		SiEnergyFluct * _fluctuate;
		
		// Root output
#ifdef ROOT_OUTPUT_LAND		
		TFile * _rootFile;
		TTree * _tree;

		std::map<int,double> _variables; //!< Variable to be stored
#endif
		
	private:
		
		double _timeCPU;                //!< CPU time
		int    _nRun;                   //!< Run number
		int    _nEvent;                 //!< Event number
		
}; // Class
/** @} */

} // Namespace

/** @} */
#endif // SISTRDIGI_H
