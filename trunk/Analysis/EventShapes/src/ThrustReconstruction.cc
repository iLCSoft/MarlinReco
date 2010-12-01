#include "ThrustReconstruction.h"
#include <iostream>
#include <fstream>
#include <vector>

// #include <CLHEP/Vector/ThreeVector.h>
// #include <CLHEP/Random/RanluxEngine.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCIO.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCTOOLS.h> 


using namespace lcio ;
using namespace marlin ;
using namespace std ;


ThrustReconstruction aThrustReconstruction ;


ThrustReconstruction::ThrustReconstruction() 
   : Processor("ThrustReconstruction") {
  
  // modify processor description
  _description = "Calculates thrust axis and thrust value of event using different algorithms" ;
  
  // register steering parameters: 
  // name, description, class-variable, default value

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			   "inputCollectionName" ,
			   "Name of collection of reconstructed particles used for thrust reconstruction"  ,
			   _inputCollectionName ,
			   std::string("SelectedReconstructedParticle") ) ;

  registerProcessorParameter( "typeOfThrustFinder" ,
      "Type of thrust reconstruction algorithm to be used:\n#\t1 : Tasso algorithm\n#\t2 : JetSet algorithm"  ,
      _typeOfThrustFinder , 2 ) ;

}

void ThrustReconstruction::init() { 
  _min = 2; _max = 0;
  // usually a good idea to
  printParameters() ;

  // config ranlux 
  filename = "Ranlux.coonf";
  ifstream rndcfgfile( filename.c_str() );
  if (!rndcfgfile)
    {
      long int ss=1234;
      myrnd.setSeeds(&ss,4);
      myrnd.showStatus();
    }
  else
    {
      rndcfgfile.close();
      myrnd.restoreStatus(filename.c_str());
      myrnd.showStatus();
    } // if file not exist

}

void ThrustReconstruction::processRunHeader( LCRunHeader* run) { 
  // run->parameters().setValue("thrust",12300321);
} 

void ThrustReconstruction::processEvent( LCEvent * evt ) { 

  // get pointer to collection vec of input particles  
  _inParVec = evt->getCollection(_inputCollectionName) ;

  if (_inParVec->getTypeName()!=LCIO::RECONSTRUCTEDPARTICLE) 
    {
      std::stringstream errorMsg;
      errorMsg << "\nProcessor: ThrustReconstruction \n" <<
        "Collection is of wrong type (" << _inParVec->getTypeName() <<
        "). Processor requires collection tpye " << LCIO::RECONSTRUCTEDPARTICLE 
        << "\n" ;
      throw Exception(errorMsg.str());
    }

  // Clear Vector of Hep3Vectors, to hold only momenta of this event.
  if (!_partMom.empty()) _partMom.clear();

  for (int n=0;n<_inParVec->getNumberOfElements() ;n++) 
    {
      ReconstructedParticle* aPart = dynamic_cast<ReconstructedParticle*>( _inParVec->getElementAt(n) );
      if ( aPart == NULL )
	throw Exception( std::string("Particle in ReconstructedParticle collection is not ReconstructedParticle") );
      
      const double* partMom = aPart->getMomentum();
      _partMom.push_back( Hep3Vector(partMom[0], partMom[1], partMom[2]) );
    } // for n

  // Reset the Class variables for Output
  _principleThrustValue = -1;
  _majorThrustValue     = -1;
  _minorThrustValue     = -1;
  _principleThrustAxis.set(0,0,0);
  _majorThrustAxis.set(0,0,0);
  _minorThrustAxis.set(0,0,0);

  // Switch to the desired type of thrust finder
  if (_typeOfThrustFinder == 1) 
    {
      TassoThrust();
    } 
  else if (_partMom.size()<=1)
    {
      TassoThrust();
    }
  else if (_typeOfThrustFinder == 2) 
    {
      JetsetThrust();
    }
  // ###write
  //    evt->parameters().setValue("thrust",_principleThrustValue);

  FloatVec thrax;

  thrax.clear();
  thrax.push_back(_principleThrustAxis.x());
  thrax.push_back(_principleThrustAxis.y());
  thrax.push_back(_principleThrustAxis.z());

  _inParVec->parameters().setValue("principleThrustValue",_principleThrustValue);
  _inParVec->parameters().setValues("principleThrustAxis",thrax);

  if (_typeOfThrustFinder == 2) 
    {
      thrax.clear();
      thrax.push_back(_majorThrustAxis.x());
      thrax.push_back(_majorThrustAxis.y());
      thrax.push_back(_majorThrustAxis.z());

      _inParVec->parameters().setValue("majorThrustValue",_majorThrustValue);
      _inParVec->parameters().setValues("majorThrustAxis",thrax);

      thrax.clear();
      thrax.push_back(_minorThrustAxis.x());
      thrax.push_back(_minorThrustAxis.y());
      thrax.push_back(_minorThrustAxis.z());

      _inParVec->parameters().setValue("minorThrustValue",_minorThrustValue);
      _inParVec->parameters().setValues("minorThrustAxis",thrax);

      float Oblateness;
      Oblateness = _majorThrustValue - _minorThrustValue;
      _inParVec->parameters().setValue("Oblateness",Oblateness);
      if ( (_majorThrustValue < 0) || (_minorThrustValue < 0) )
	{
	  _inParVec->parameters().setValue("Oblateness",-1);
	}
    }

  // Ausgabe der Werte:
    cout << " thrust: " << _principleThrustValue << 
      " TV: " << _principleThrustAxis << endl; 
        cout << "  major: " << _majorThrustValue << 
        " TV: " << _majorThrustAxis << endl; 
        cout << "  minor: " << _minorThrustValue << 
        " TV: " << _minorThrustAxis << endl; 
  if (_principleThrustValue >= _max) _max = _principleThrustValue;
  if (_principleThrustValue <= _min) _min = _principleThrustValue;
}

void ThrustReconstruction::end(){ 
  
  //  cout << "max, min: " << _max << " " << _min << endl;

  myrnd.saveStatus( filename.c_str() );
}

int ThrustReconstruction::TassoThrust(){
  int ThrustError = 0; 
  Hep3Vector tvec;

  // No particle in Event: Error
  if (_inParVec->getNumberOfElements()<=0) 
    {
      ThrustError = -1;
      _principleThrustValue = 0;
      _principleThrustAxis.set(0,0,0);
    }
  // only one Particle in Event: Thrust direction = direction of particle
  else if (_inParVec->getNumberOfElements()==1) 
    {
      _principleThrustValue = 1; 
      _principleThrustAxis = _partMom[0];
    }
  else 
    {
      Hep3Vector ptm, ptot, pt;
      std::vector<Hep3Vector> pc; 
      float sp,u, pp, tmax, t;

      sp = 0;
      for (int i=0;i < _inParVec->getNumberOfElements();i++)
	{
	  pp = _partMom[i].mag();
	  sp += pp;
	  ptot += _partMom[i];  
	} // for i
	// ###
      for (int m = 0; m <= 2; m++ )
	ptot[m] *= 0.5; 
      tmax = 0;
      for (int k = 1; k < _inParVec->getNumberOfElements(); k++)
	{
	  for (int j = 0; j <= k-1;j++)
	    {
              // cross product
	      tvec = _partMom[j].cross(_partMom[k]);
	      pt = -1 * ptot;
	      for (int l = 0; l < _inParVec->getNumberOfElements(); l++)
		{
		  if (l==k) continue;
		  if (l==j) continue;
		  u = _partMom[l] * tvec;
		  if (u<0) continue;
		  pt += _partMom[l];
	        } // for l

	      while(!pc.empty())
		{
		  pc.pop_back();
		}
		// note: the order is important!!!
	      pc.push_back(pt);
	      pc.push_back(pt + _partMom[k]);
	      pc.push_back(pt + _partMom[j]);
	      pc.push_back(pc[2] + _partMom[k]);
	      for (int m = 0; m <= 3; m++ )
		{
		  t = pc[m].mag2();
		  if (t <= tmax) continue;
		  tmax = t;
		  ptm = pc[m];
		} // for m 
	    } // for j
	} // for k
      _principleThrustValue = 2 * sqrt(tmax) / sp;
      tvec = ptm;
    } // end else 

  // Normalization of thrust vector
  double ax = 0; 
  ax = tvec.mag();
  if (ax != 0) 
    {
      ax = 1/ax;
      _principleThrustAxis = ax * tvec; 
    } 
  else 
    { 
      ThrustError = -1; 
      _principleThrustValue = -1;
      _principleThrustAxis.set(0,0,0);
    }
  return ThrustError; 
}

int ThrustReconstruction::JetsetThrust(){
  const int nwork=11,iFastMax = 4,iGood=2;  
  const float dConv=0.0001; // 0.0001
  int sgn;
  double theta=0,phi=0;
  double thp,thps,tds,tmax,dOblateness;
  vector<Hep3Vector> TAxes(3),Fast(iFastMax+1),Workv(nwork);
  vector<double> Workf(nwork),dThrust(3);
  Hep3Vector tdi,tpr,mytest;

  tmax = 0;
  for ( unsigned int i=0; i < _partMom.size(); i++)
    tmax += _partMom[i].mag();

  // pass = 0: find thrust axis
  // pass = 1: find major axis
  for ( int pass=0; pass <= 1; pass++ ) 
    {
      if ( pass == 1 )
	{
	  phi   = TAxes[0].phi();
	  theta = TAxes[0].theta();
	  for ( unsigned  int i = 0;i < _partMom.size(); i++)
	    {
	      _partMom[i].rotateZ(-phi);
	      _partMom[i].rotateY(-theta);
	    }
	  TAxes[0].set(0,0,1);
	} // if pass == 1

      // Find the ifast highest momentum particles and
      // put the highest in Fast[0], next in Fast[1],....Fast[iFast-1].
      // Fast[iFast] is just a workspace.

      for ( unsigned  int i = 0; i < Fast.size(); i++ )
	Fast[i].set(0,0,0);

      for ( unsigned int i = 0; i < _partMom.size(); i++ )
	{
	  for ( int ifast = iFastMax -1; ifast >= 0 ; ifast-- ) 
	    {
	      if (_partMom[i].mag2() > Fast[ifast].mag2() ) 
		{
		  Fast[ifast + 1] = Fast[ifast]; 
		  if (ifast == 0) Fast[ifast] = _partMom[i]; 
		} 
	      else 
		{
		  Fast[ifast + 1] = _partMom[i]; 
		  break;
		} // if p>p_fast
	    } // for ifast 
	} // for i 

      // Find axis with highest thrust (case 0)/ highest major (case 1).

      for ( unsigned int iw = 0; iw < Workv.size(); iw++ ) 
	{
	  Workf[iw] = 0.;
	}
      int p = (int) min( iFastMax,_partMom.size() ) - 1 ;
      int nc = 1 << p;
      for ( int n = 0; n < nc; n++ )
	{
	  tdi.set(0,0,0);
	  for (int i = 0; i < min(iFastMax,nc) ; i++)
	    {
	      if ( (1 << (i+1)) * ( (n + (1<<i)) / (1<<(i+1)) ) >= n+1) //i+1 
		{ sgn = -1;} else {sgn = 1;}
	      tdi += sgn*Fast[i];
	      if (pass==1) tdi.setZ(0);
	    } // for i 
	  tds = tdi.mag2(); 
	  for ( int iw = (int) min(n,9); iw >= 0; iw-- )
	    {
	      if (tds > Workf[iw])
		{
		  Workf[iw+1] = Workf[iw]; 
		  Workv[iw+1] = Workv[iw]; 
		  if (iw == 0) 
		    { Workv[iw] = tdi; Workf[iw] = tds;} 
  		}  
	      else // if tds 
		{
		  Workv[iw+1] = tdi;
		  Workf[iw+1] = tds;
  		} // if tds 
	    } // for iw
	} // for n 

      // Iterate direction of axis until stable maximum.

      dThrust[pass] = 0;
      int nagree = 0;
      for ( int iw = 0; iw < min(nc,10) && nagree < iGood; iw++ )
	{
	  thp = 0;
	  thps = -99999.;
	  while ( thp > thps + dConv )
	    {
	      thps = thp;
	      if ( thp <= 1E-10 )
		{ tdi = Workv[iw]; } else { tdi=tpr; }
	      tpr.set(0,0,0);
	      for ( unsigned int i = 0; i < _partMom.size(); i++ )
		{
		  sgn = (int) sign(1,tdi.dot(_partMom[i]));
		  tpr += sgn*_partMom[i];
		  if (pass == 1) { tpr.setZ(0); } // ###
		} // for i 
	      thp = tpr.mag()/tmax;
	    } // while 
	  // Save good axis. Try new initial axis until enough
	  // tries agree.
	  if ( thp < dThrust[pass] - dConv ) continue;
	  if ( thp > dThrust[pass] + dConv )
	    {
	      nagree = 0;
	      // 	      if (myrnd.flat() > 0.49999)
	      // 		{sgn = 1;} else {sgn=-1;}
	      sgn = 1; 
	      TAxes[pass] = sgn*tpr/(tmax*thp);
	      dThrust[pass] = thp;
	    } // if thp
	  nagree++;
	} // for iw (2)
    } // for pass ...

  // Find minor axis and value by orthogonality.
  if (myrnd.flat() > 0.49999)
    {sgn = 1;} else {sgn=-1;}
  TAxes[2].set( -sgn*TAxes[1].y(), sgn*TAxes[1].x(), 0);
  thp = 0.;
  for ( unsigned int i = 0; i < _partMom.size(); i++ )
    {
      thp += fabs(TAxes[2].dot(_partMom[i]) );
    } // for i 
  dThrust[2] = thp/tmax;

  // Rotate back to original coordinate system.
  for ( unsigned int i = 0;i < TAxes.size(); i++)
    {
      TAxes[i].rotateY(theta); 
      TAxes[i].rotateZ(phi); 
    }
  dOblateness = dThrust[1] - dThrust[2];

  _principleThrustValue = dThrust[0];
  _majorThrustValue     = dThrust[1];
  _minorThrustValue     = dThrust[2];
  _principleThrustAxis  =   TAxes[0];
  _majorThrustAxis      =   TAxes[1];
  _minorThrustAxis      =   TAxes[2];


  return  0;
}

//______________________________________________________________
// helper function to get sign of b
double ThrustReconstruction::sign(double a, double b) 
{
  if ( b < 0 ) 
    { return -fabs(a); } else { return fabs(a); }
}
//______________________________________________________________
double ThrustReconstruction::min(double a, double b) 
{
  if ( a < b ) 
    { return a; } else { return b; }
}
//______________________________________________________________
