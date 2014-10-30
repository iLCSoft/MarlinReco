#ifndef _ScintillatorPpdDigi_h_
#define _ScintillatorPpdDigi_h_

class ScintillatorPpdDigi {

 public:

  ScintillatorPpdDigi();
  ~ScintillatorPpdDigi(){}

  // expected # photoelectrons / MIP
  void setPEperMIP(float x)           {_pe_per_mip=x;}

  // calibration factor from input hit energy to MIPs
  void setCalibMIP(float x)           {_calib_mip=x;}

  // #pixels of PPD
  void setNPix(int x)                 {_npix=x;}

  // random miscalibration of total #pixels (as a fraction of pixel number: 0.05 = 5% miscalibration)
  void setRandomMisCalibNPix(float x) {_misCalibNpix=x;}

  // spread in pixel capacitance (as a fraction: 0.05 = 5% spread)
  void setPixSpread(float x)          {_pixSpread=x;}

  // electronics noise (in MIP units)
  void setElecNoise(float x)          {_elecNoise=x;}

  float getDigitisedEnergy( float energy );

  void printParameters();

 private:

  float _pe_per_mip;
  float _calib_mip;
  float _npix;
  float _misCalibNpix;
  float _pixSpread;
  float _elecNoise;

};

#endif
