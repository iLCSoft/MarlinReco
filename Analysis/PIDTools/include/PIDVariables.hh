/*
 * PIDVariables.hh
 *
 * PIDVariable
 * Helper class that contains and manages a variable
 * used by PIDTools processors.
 *
 * PIDVariables
 * Helper class that creates and manages a map of PIDVariable
 * objects, Update()s it from a reconstructed particle etc.
 * This class takes care and standardises the calculation of
 * variables used for PID.
 *
 *  Created on: Dec 23, 2015
 *      Author: S. Lukic
 */

#ifndef PIDVARIABLES_H_
#define PIDVARIABLES_H_ 1

#include "TVector3.h"
#include "TString.h"
#include <LCIOSTLTypes.h>
#include "EVENT/Cluster.h"
#include "EVENT/Track.h"
#include "EVENT/ReconstructedParticle.h"

#include "PIDParticles.hh"

class TRandom3;

using EVENT::FloatVec;


/*****************************************************************
 *
 * Individual variables. Base class first, then specific
 *
 ****************************************************************/

class PIDVariable_base {
public:
  PIDVariable_base(const char*name, const char *description, const char *unit) :
    _value(0.), _name(name), _description(description), _unit(unit)
    {};
  virtual ~PIDVariable_base() {};

  // Binary masks for interpreting the return code
  // of the Update() function
  static const short MASK_EmptyClusters   = 1 ;
  static const short MASK_EmptyTracks     = 1 << 1;
  static const short MASK_EmptyShapes     = 1 << 2;
  static const short MASK_ZerodEdx        = 1 << 3;
  static const short MASK_InvalidMomentum = 1 << 4;


  float Value() const { return _value; };
  const char *Name() const { return _name; };
  const float *Address() const { return &_value; };

  const char *AxisTitle() const { if(_unit[0] == '\0') return _description;
                                  else return Form("%s (%s)", _description, _unit); };
  const char *Description() const { return _description; };
  const char *Unit() const { return _unit; } ;

  virtual int Update(EVENT::ReconstructedParticle*);
  virtual int Update(const EVENT::ClusterVec cluvec, const EVENT::TrackVec trax, const TVector3 p3) = 0;

  virtual void SetOutOfRange() = 0;

  static double BetheBloch(const PIDParticles::PIDParticle_base* hypothesis, const float p);
  // Some PIDVariables have occasional discrete values in VALID events.
  // MVA does not like this, so we will add tiny smearing in such cases
  // If varRand is not created by the user, no smearing will be added
  static TRandom3 *varRand;

protected:
  float _value;
  const char *_name;
  const char *_description;
  const char *_unit;
};



class PID_CaloTotal : public PIDVariable_base
{
public:
  PID_CaloTotal();
  virtual int Update(const EVENT::ClusterVec, const EVENT::TrackVec, const TVector3 p3);
  virtual void SetOutOfRange() { _value = -1.; }

  static const float pCut; // Minimum permissible tracker momentum
};

class PID_CaloEFrac: public PIDVariable_base
{
public:
  PID_CaloEFrac();
  virtual int Update(const EVENT::ClusterVec, const EVENT::TrackVec, const TVector3 p3);
  virtual void SetOutOfRange() { _value = -1.; }

  static const float caloCut;
};

class PID_CaloMuSys : public PIDVariable_base
{
public:
  PID_CaloMuSys();
  virtual int Update(const EVENT::ClusterVec, const EVENT::TrackVec, const TVector3 p3);
  virtual void SetOutOfRange() { _value = -1.; }
  static const float muSysCut;
};

class PID_CluShapeChi2: public PIDVariable_base
{
public:
  PID_CluShapeChi2();
  virtual int Update(const EVENT::ClusterVec, const EVENT::TrackVec, const TVector3 p3);
  virtual void SetOutOfRange() { _value = -1.; }
};

class PID_CluShapeLDiscr: public PIDVariable_base
{
public:
  PID_CluShapeLDiscr();
  virtual int Update(const EVENT::ClusterVec, const EVENT::TrackVec, const TVector3 p3);
  virtual void SetOutOfRange() { _value = -FLT_MAX; }
};

class PID_CluShapeTDiscr: public PIDVariable_base
{
public:
  PID_CluShapeTDiscr();
  virtual int Update(const EVENT::ClusterVec, const EVENT::TrackVec, const TVector3 p3);
  virtual void SetOutOfRange() { _value = -1.; }
};

class PID_CluShapeXl20: public PIDVariable_base
{
public:
  PID_CluShapeXl20();
  virtual int Update(const EVENT::ClusterVec, const EVENT::TrackVec, const TVector3 p3);
  virtual void SetOutOfRange() { _value = -1.; }
};


class PID_dEdxChi2: public PIDVariable_base
{
public:
  PID_dEdxChi2(const PID_dEdxChi2&);
  PID_dEdxChi2(const PIDParticles::PIDParticle_base* hypothesis, const float dEdx_MIP=1.35e-7);
  ~PID_dEdxChi2();
  virtual int Update(const EVENT::ClusterVec, const EVENT::TrackVec, const TVector3 p3);
  virtual void SetOutOfRange() { _value = -1.e5; }

//private:
  const PIDParticles::PIDParticle_base* _hypothesis;
  const double _dEdx_MIP;
};

class PID_dEdxLogChi2: public PIDVariable_base
{
public:
  PID_dEdxLogChi2(const PID_dEdxLogChi2&);
  PID_dEdxLogChi2(const PIDParticles::PIDParticle_base* hypothesis, const float dEdx_MIP=1.35e-7);
  ~PID_dEdxLogChi2();
  virtual int Update(const EVENT::ClusterVec, const EVENT::TrackVec, const TVector3 p3);
  virtual void SetOutOfRange() { _value = -1.e3; }

//private:
  const PIDParticles::PIDParticle_base* _hypothesis;
  const double _dEdx_MIP;
};



/*****************************************************************
 *
 * Sets of variables. Base class first, then variations
 *
 ****************************************************************/


class PIDVariables_base {

public:
  PIDVariables_base();
  PIDVariables_base(EVENT::ReconstructedParticle*);
  virtual ~PIDVariables_base() = 0;

  typedef std::vector<PIDVariable_base*> VarVec;
  typedef PIDParticles::ParticleMap ParticleMap;

  const VarVec* GetVariables() const { return &_varVec; };
  float GetP() const { return _p; }

  virtual int Update(EVENT::ReconstructedParticle*);
  virtual int Update(const EVENT::ClusterVec, const EVENT::TrackVec, const TVector3 p);
  virtual void SetOutOfRange();

  virtual void ClearVars();

protected:
  VarVec _varVec;
  float _p;
  virtual void Populate() = 0;
};


/***  PIDVariables for the LikelihoodPIDProcessor ***/

class PIDVariables_LLPID : public PIDVariables_base
{
public:
  PIDVariables_LLPID();
  PIDVariables_LLPID(EVENT::ReconstructedParticle*);
  virtual ~PIDVariables_LLPID();

protected:
  virtual void Populate();
};


/***  PIDVariables for the MvaPidProcessor ***/

class PIDVariables_MvaPid : public PIDVariables_base
{
public:
  PIDVariables_MvaPid();
  PIDVariables_MvaPid(EVENT::ReconstructedParticle*);
  virtual ~PIDVariables_MvaPid();

  virtual int Update(EVENT::ReconstructedParticle*);
  virtual void SetOutOfRange();

  FloatVec* GetMvaVariables() { return &_mvaVars; }

protected:
  virtual void Populate();

  // Copy of all variables for the TMVA::Reader and TMVA::Factory
  // As they do not accept const pointers
  FloatVec _mvaVars;
  void RefreshMvaVars();

};



#endif // PIDVARIABLES_H_
