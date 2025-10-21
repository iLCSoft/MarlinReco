#include "Use_TrueJet.h"

// ----- include for verbosity dependend logging ---------

#ifdef MARLIN_USE_AIDA
#include <AIDA/ICloud1D.h>
#include <AIDA/IHistogramFactory.h>
#include <marlin/AIDAProcessor.h>
// #include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA

#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace lcio;
using namespace marlin;
using namespace std;

Use_TrueJet aUse_TrueJet;

Use_TrueJet::Use_TrueJet() : Processor("Use_TrueJet") {

  // modify processor description
  _description = "Use_TrueJet does whatever it does ...";

  // register steering parameters: name, description, class-variable, default value

  // Inputs: MC-particles, Reco-particles, the link between the two

  registerInputCollection(LCIO::MCPARTICLE, "MCParticleCollection", "Name of the MCParticle collection",
                          _MCParticleColllectionName, std::string("MCParticlesSkimmed"));

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "RecoParticleCollection",
                          "Name of the ReconstructedParticles input collection", _recoParticleCollectionName,
                          std::string("PandoraPFOs"));

  registerInputCollection(LCIO::LCRELATION, "RecoMCTruthLink", "Name of the RecoMCTruthLink input collection",
                          _recoMCTruthLink, std::string("RecoMCTruthLink"));

  // Inputs: True jets (as a recoparticle, will be the sum of the _reconstructed particles_
  // created by the true particles in each true jet, in the RecoMCTruthLink sense.
  // link jet-to-reco particles, link jet-to-MC-particles.

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "TrueJets", "Name of the TrueJetCollection input collection",
                          _trueJetCollectionName, std::string("TrueJets"));

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "FinalColourNeutrals",
                          "Name of the FinalColourNeutralCollection input collection",
                          _finalColourNeutralCollectionName, std::string("FinalColourNeutrals"));

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "InitialColourNeutrals",
                          "Name of the InitialColourNeutralCollection input collection",
                          _initialColourNeutralCollectionName, std::string("InitialColourNeutrals"));

  registerInputCollection(LCIO::LCRELATION, "TrueJetPFOLink", "Name of the TrueJetPFOLink input collection",
                          _trueJetPFOLink, std::string("TrueJetPFOLink"));

  registerInputCollection(LCIO::LCRELATION, "TrueJetMCParticleLink",
                          "Name of the TrueJetMCParticleLink input collection", _trueJetMCParticleLink,
                          std::string("TrueJetMCParticleLink"));

  registerInputCollection(LCIO::LCRELATION, "FinalElementonLink rueJetMCParticleLink",
                          "Name of the  FinalElementonLink input collection", _finalElementonLink,
                          std::string("FinalElementonLink"));

  registerInputCollection(LCIO::LCRELATION, "InitialElementonLink",
                          "Name of the  InitialElementonLink input collection", _initialElementonLink,
                          std::string("InitialElementonLink"));

  registerInputCollection(LCIO::LCRELATION, "FinalColourNeutralLink",
                          "Name of the  FinalColourNeutralLink input collection", _finalColourNeutralLink,
                          std::string("FinalColourNeutralLink"));

  registerInputCollection(LCIO::LCRELATION, "InitialColourNeutralLink",
                          "Name of the  InitialColourNeutralLink input collection", _initialColourNeutralLink,
                          std::string("InitialColourNeutralLink"));
}

void Use_TrueJet::init() {

  streamlog_out(DEBUG6) << "   init called  " << std::endl;

  // usually a good idea to
  printParameters();

  _nRun = 0;
  _nEvt = 0;
}

void Use_TrueJet::processRunHeader(LCRunHeader*) { _nRun++; }

void Use_TrueJet::processEvent(LCEvent* evt) {

  streamlog_out(MESSAGE) << std::endl;
  streamlog_out(MESSAGE) << " ==================================================== " << std::endl;
  streamlog_out(MESSAGE) << " processing event: " << evt->getEventNumber() << "   in run:  " << evt->getRunNumber()
                         << std::endl;
  streamlog_out(MESSAGE) << " ==================================================== " << std::endl;
  streamlog_out(MESSAGE) << std::endl;

  LCCollection* rmclcol = NULL;
  try {
    rmclcol = evt->getCollection(_recoMCTruthLink);
  } catch (lcio::DataNotAvailableException& e) {
    rmclcol = NULL;
  }
  // tj is a pointer to a Trujet_Parser, with the data of this processor object:
  TrueJet_Parser* tj = this;
  // this method gets all the collections needed + initialises a few convienent variables.
  tj->getall(evt);
  if (tjcol == NULL) {
    return;
  }

  streamlog_out(DEBUG5) << " Number of jets is " << tj->njets() << endl;
  if (tj->njets() == 0) {
    return;
  }
  int njets = tj->njets();
  IntVec sibl[50]; // This will be used to store siblings of each jet
  double total_Etrue = 0;
  double total_Eisr = 0;
  double total_Eowl = 0;
  int leptons = 0;
  bool cluster_in_event = false;
  bool top_in_event = false;
  bool gluonsplit_in_event = false;
  for (int ijet = 0; ijet < njets; ijet++) { // start jet loop

    // enery, moemtum, type of the jets. Choose between true, seen or true-of-seen values
    // type: 1 string, 2 lepton, 3 cluster, 4 ISR, 5 overlay.
    if (abs(type_jet(ijet) % 100) == 3) {
      cluster_in_event = true;
    }
    if (abs(type_jet(ijet)) > 100) {
      gluonsplit_in_event = true;
    }
    if (abs((initial_elementon(ijet) != NULL ? initial_elementon(ijet)->getPDG() : 0)) == 6) {
      top_in_event = true;
    }
    streamlog_out(DEBUG5) << endl;
    streamlog_out(DEBUG5) << " jet " << ijet << " has type " << type_jet(ijet) << endl;
    streamlog_out(DEBUG5) << " Seen :         " << Eseen(ijet);
    for (int jj = 0; jj < 3; jj++) {
      streamlog_out(DEBUG5) << " " << pseen(ijet)[jj];
    };
    streamlog_out(DEBUG5) << " " << Mseen(ijet) << endl;
    streamlog_out(DEBUG4) << "               ";
    for (int jj = 0; jj < 4; jj++) {
      streamlog_out(DEBUG4) << " " << p4seen(ijet)[jj];
    };
    streamlog_out(DEBUG4) << endl;
    streamlog_out(DEBUG5) << " True :         " << Etrue(ijet);
    for (int jj = 0; jj < 3; jj++) {
      streamlog_out(DEBUG5) << " " << ptrue(ijet)[jj];
    };
    streamlog_out(DEBUG5) << " " << Mtrue(ijet) << endl;
    streamlog_out(DEBUG4) << "               ";
    for (int jj = 0; jj < 4; jj++) {
      streamlog_out(DEBUG4) << " " << p4true(ijet)[jj];
    };
    streamlog_out(DEBUG4) << endl;
    if (rmclcol != NULL) {
      streamlog_out(DEBUG5) << " True of seen : " << Etrueseen(ijet);
      for (int jj = 0; jj < 3; jj++) {
        streamlog_out(DEBUG5) << " " << ptrueseen(ijet)[jj];
      };
      streamlog_out(DEBUG5) << " " << Mtrueseen(ijet) << endl;
      streamlog_out(DEBUG4) << "               ";
      for (int jj = 0; jj < 4; jj++) {
        streamlog_out(DEBUG4) << " " << p4trueseen(ijet)[jj];
      };
      streamlog_out(DEBUG4) << endl;
    }
    streamlog_out(DEBUG5) << " elementon :    " << Equark(ijet);
    for (int jj = 0; jj < 3; jj++) {
      streamlog_out(DEBUG5) << " " << pquark(ijet)[jj];
    };
    streamlog_out(DEBUG5) << " " << Mquark(ijet) << endl;
    streamlog_out(DEBUG4) << "               ";
    for (int jj = 0; jj < 4; jj++) {
      streamlog_out(DEBUG4) << " " << p4quark(ijet)[jj];
    };
    streamlog_out(DEBUG4) << endl;
    streamlog_out(DEBUG5) << endl;

    // make some global sums:
    if (abs(type_jet(ijet)) % 100 < 4 || abs(type_jet(ijet)) % 100 == 6) {
      total_Etrue += Etrue(ijet);
      if (abs(type_jet(ijet)) % 100 == 2) {
        leptons++;
      }
    } else if (abs(type_jet(ijet)) % 100 == 4) {
      total_Eisr += Etrue(ijet);
    } else if (abs(type_jet(ijet)) % 100 == 5) {
      if (rmclcol != NULL) {
        total_Eowl += Etrueseen(ijet);
      }
    }

    // get reco particles of the jet, as a ReconstructedParticleVec.
    streamlog_out(DEBUG5) << " Number of PFOs used : " << seen_partics(ijet).size() << endl;
    if (seen_partics(ijet).size() > 0) {
      streamlog_out(DEBUG4) << " list of PFOs with energy and jet number " << endl;
      ReconstructedParticleVec recos = seen_partics(ijet); // note syntax.
      double recoE = 0;
      for (unsigned kk = 0; kk < recos.size(); kk++) {
        streamlog_out(DEBUG4) << recos[kk] << " " << recos[kk]->getEnergy() << " " << recojet(recos[kk]) << endl;
        recoE += recos[kk]->getEnergy();
      }
      streamlog_out(DEBUG5) << " Total energy: " << recoE << endl;
    }

    // get true particles of the jet, as a MCParticleVec.
    streamlog_out(DEBUG5) << " Number of MCPs used : " << true_partics(ijet).size() << endl;
    if (true_partics(ijet).size() > 0) {
      streamlog_out(DEBUG4) << " list of MCPs with energy and jet number " << endl;
      MCParticleVec mcps = true_partics(ijet); // note syntax.
      for (unsigned kk = 0; kk < mcps.size(); kk++) {
        streamlog_out(DEBUG4) << mcps[kk] << " " << mcps[kk]->getEnergy() << " " << mcpjet(mcps[kk]) << endl;
      }
    }

    // study the siblings, ie. the jet(s) that forms a colour neutral with the current one, either final (ie. at the end
    // of the parton-shower), or initial (ie. at the beginning). The siblings are given as an IntVec of indicies into
    // the jets ReconstructedParticleVec

    // make a copy. final_siblings clears the vector at each call. Note syntax. aibl is IntVec sibl[50].

    sibl[ijet] = final_siblings(ijet);

    // print from sibl
    streamlog_out(DEBUG3) << " Number final Siblings to this jet :" << sibl[ijet].size() << endl;
    if (sibl[ijet].size() > 0) {
      streamlog_out(DEBUG3) << " Siblings to this jet are/is :";
      for (unsigned jj = 0; jj < sibl[ijet].size(); jj++) {
        streamlog_out(DEBUG3) << " " << sibl[ijet][jj];
      };
      streamlog_out(DEBUG3) << endl;
    }
    streamlog_out(DEBUG5) << " Number final siblings to this jet :" << final_siblings(ijet).size() << endl;
    if (final_siblings(ijet).size() > 0) {
      streamlog_out(DEBUG5) << " Siblings to this jet are/is :";
      for (unsigned jj = 0; jj < final_siblings(ijet).size(); ++jj) {
        streamlog_out(DEBUG5) << " " << final_siblings(ijet)[jj];
      };
      streamlog_out(DEBUG5) << endl;
    }
    streamlog_out(DEBUG5) << " Number initial siblings to this jet :" << initial_siblings(ijet).size() << endl;
    if (initial_siblings(ijet).size() > 0) {
      streamlog_out(DEBUG5) << " Siblings to this jet are/is :";
      for (unsigned jj = 0; jj < initial_siblings(ijet).size(); ++jj) {
        streamlog_out(DEBUG5) << " " << initial_siblings(ijet)[jj];
      };
      streamlog_out(DEBUG5) << endl;
    }

    // study colour neutrals the jet comes from. Get the list (as index in the finalcns/initialcns
    // ReconstructedParticleVec:s)
    streamlog_out(DEBUG5) << " final and initial colour-neutral number of the jet " << final_cn(ijet) << " "
                          << initial_cn(ijet) << endl;

    // study initial/final particle (quark/lepton/ISR photon) the jet comes from. as the corresponding MCParticle* .
    // Here we just print the PDG code. Note the check for a NULL pointer (overlay jet !)
    streamlog_out(DEBUG5) << " flavour of final and initial quark, lepton or photon  "
                          << (final_elementon(ijet) != NULL ? final_elementon(ijet)->getPDG() : 0) << " "
                          << (initial_elementon(ijet) != NULL ? initial_elementon(ijet)->getPDG() : 0) << endl;

  } // end of jet loop

  // (Just to check that  sibl really kept the values)
  for (int ijet = 0; ijet <= njets; ijet++) {
    streamlog_out(DEBUG3) << ijet << endl;
    streamlog_out(DEBUG3) << " Number Siblings to this jet :" << sibl[ijet].size() << endl;
    if (sibl[ijet].size() > 0) {
      streamlog_out(DEBUG3) << " Siblings to this jet are/is :";
      for (unsigned jj = 0; jj < sibl[ijet].size(); jj++) {
        streamlog_out(DEBUG3) << " " << sibl[ijet][jj];
      };
      streamlog_out(DEBUG3) << endl;
    }
  }

  double total_Ecn = 0;
  int _nicn = nicn();
  for (int iicn = 0; iicn < _nicn; iicn++) {
    total_Ecn += E_icn(iicn);
  }
  //
  streamlog_out(DEBUG5) << " Number of final colour-neutrals is " << nfcn() << endl;
  int _nfcn = nfcn();
  for (int ifcn = 0; ifcn < _nfcn; ifcn++) { // loop over final colour neutrals

    // get type (same meaning as for jets) and E/P and M.
    streamlog_out(DEBUG5) << endl;
    streamlog_out(DEBUG5) << " fcn " << ifcn << " has type " << type_fcn_parent(ifcn) << endl;
    streamlog_out(DEBUG5) << " True :         " << E_fcn(ifcn);
    for (int jj = 0; jj < 3; jj++) {
      streamlog_out(DEBUG5) << " " << p_fcn(ifcn)[jj];
    };
    streamlog_out(DEBUG5) << " " << M_fcn(ifcn) << endl;
    streamlog_out(DEBUG4) << "               ";
    for (int jj = 0; jj < 4; jj++) {
      streamlog_out(DEBUG4) << " " << p4_fcn(ifcn)[jj];
    };
    streamlog_out(DEBUG4) << endl;

    // the other way wrt the jet-loop: get list of jets comming from this colour-neutral. As an InTVec
    //  with index into the jets ReconstructedParticleVec.
    streamlog_out(DEBUG5) << " Number of jets from this cn:" << jets_of_final_cn(ifcn).size() << endl;
    if (jets_of_final_cn(ifcn).size() > 0) {
      if (jets_of_final_cn(ifcn).size() > 1) {
        streamlog_out(DEBUG5) << " jets are :";
      } else {
        streamlog_out(DEBUG5) << " jet is :  ";
      }
      for (unsigned jj = 0; jj < jets_of_final_cn(ifcn).size(); ++jj) {
        streamlog_out(DEBUG5) << " " << jets_of_final_cn(ifcn)[jj];
      };
      streamlog_out(DEBUG5) << endl;
    }
    IntVec fcn_jets = jets_of_final_cn(ifcn); // make a copy of the list
    {
      // when comparing final cn:es with jets, *dont* include FSRs: they enter into the jet, but not the cn,
      // since they were radiated earlier in the parton-shower.
      _COUNT_FSR = 0; // this tells to ignore FSR
      double p4_fr_jets[4] = {0, 0, 0, 0};
      for (unsigned kk = 0; kk < fcn_jets.size(); kk++) {
        for (int jj = 0; jj < 4; jj++) {
          p4_fr_jets[jj] += p4true(fcn_jets[kk])[jj]; // 4-moentum of the cn, using the values from the true
                                                      // particles contributing to the jets of the colour neutral
        }
      }
      _COUNT_FSR = 1; // set it back
      // consitancy check: Is Eand p calculated in different ways (sum of true particles, or from the colour-neutral
      // itself) the same ?
      double M = sqrt(p4_fr_jets[0] * p4_fr_jets[0] - p4_fr_jets[1] * p4_fr_jets[1] - p4_fr_jets[2] * p4_fr_jets[2] -
                      p4_fr_jets[3] * p4_fr_jets[3]);
      streamlog_out(DEBUG5) << " true p4 and mass of these:      ";
      for (int jj = 0; jj < 4; jj++) {
        streamlog_out(DEBUG5) << " " << p4_fr_jets[jj];
      };
      streamlog_out(DEBUG5) << " " << M << endl;
      // Normally, the mass should be identical. However, due to the B-field, decyas of longlived charged
      // particles have their momentum *direction* changed, so check energy as well
      float margin = 0.007;
      // If this is a tau cn, be more forgiving. The theory is that there is an issue with the
      // the crossing-angle boost that becomes too big in aa/BB events which effects tau:s most
      if (abs(pdg_fcn_parent(ifcn)) == 15 || abs(pdg_fcn_parent(ifcn)) == 16) {
        margin = 0.02;
      }
      if (abs(p4_fr_jets[0] - E_fcn(ifcn)) / E_fcn(ifcn) > margin && p4_fr_jets[0] > 0.5) {
        if (((!cluster_in_event) && (!top_in_event)) ||
            (abs((total_Ecn - total_Eisr) - total_Etrue) / total_Etrue > 0.001)) {
          streamlog_out(WARNING) << " event/run " << evt->getEventNumber() << "/" << evt->getRunNumber()
                                 << ", final cn " << ifcn << ": E of cn:    " << E_fcn(ifcn)
                                 << " E of jets :    " << p4_fr_jets[0] << endl;
          streamlog_out(WARNING) << " event/run " << evt->getEventNumber() << "/" << evt->getRunNumber()
                                 << ", final cn " << ifcn << ": mass of cn: " << M_fcn(ifcn) << " mass of jets : " << M
                                 << endl;
          streamlog_out(WARNING) << " event/run " << " cluster_in_event : " << cluster_in_event
                                 << " gluonsplit_in_event : " << gluonsplit_in_event
                                 << "  top_in_event : " << top_in_event << endl;
          streamlog_out(WARNING) << " event/run " << " Totals : all non-isr initials = " << total_Ecn - total_Eisr
                                 << " all final non-ISR " << total_Etrue << endl;
        }
      }
    }
    {
      // now do it again, this time with the seen values.
      double p4_fr_jets[4] = {0, 0, 0, 0};
      for (unsigned kk = 0; kk < fcn_jets.size(); kk++) {
        for (int jj = 0; jj < 4; jj++) {
          p4_fr_jets[jj] += p4seen(fcn_jets[kk])[jj];
        }
      }
      streamlog_out(DEBUG5) << " seen p4 and mass of these:      ";
      for (int jj = 0; jj < 4; jj++) {
        streamlog_out(DEBUG5) << " " << p4_fr_jets[jj];
      };
      streamlog_out(DEBUG5) << " "
                            << sqrt(p4_fr_jets[0] * p4_fr_jets[0] - p4_fr_jets[1] * p4_fr_jets[1] -
                                    p4_fr_jets[2] * p4_fr_jets[2] - p4_fr_jets[3] * p4_fr_jets[3])
                            << endl;
    }
    {
      if (rmclcol != NULL) {
        // and again, with the true values of seen particles.
        double p4_fr_jets[4] = {0, 0, 0, 0};
        for (unsigned kk = 0; kk < fcn_jets.size(); kk++) {
          for (int jj = 0; jj < 4; jj++) {
            p4_fr_jets[jj] += p4trueseen(fcn_jets[kk])[jj];
          }
        }
        streamlog_out(DEBUG5) << " true-seen p4 and mass of these: ";
        for (int jj = 0; jj < 4; jj++) {
          streamlog_out(DEBUG5) << " " << p4_fr_jets[jj];
        };
        streamlog_out(DEBUG5) << " "
                              << sqrt(p4_fr_jets[0] * p4_fr_jets[0] - p4_fr_jets[1] * p4_fr_jets[1] -
                                      p4_fr_jets[2] * p4_fr_jets[2] - p4_fr_jets[3] * p4_fr_jets[3])
                              << endl;
      }
    }

    streamlog_out(DEBUG5) << " Number of elementons making this cn:" << elementons_final_cn(ifcn).size();
    streamlog_out(DEBUG5) << " Pdg of the parent " << pdg_fcn_parent(ifcn);
    streamlog_out(DEBUG5) << " Pdgs of all elementons making this cn: ";
    for (unsigned ipid = 0; ipid < pdg_fcn_comps(ifcn).size(); ipid++) {
      streamlog_out(DEBUG5) << pdg_fcn_comps(ifcn)[ipid] << " ";
    }
    streamlog_out(DEBUG5) << endl;

  } // end loop over final colour neutrals

  streamlog_out(DEBUG5) << " Number of initial  colour-neutrals is " << nicn() << endl;
  for (int iicn = 0; iicn < _nicn; iicn++) { // loop over initial colour neutrals

    // E/P/type of initial colour neutrals.
    streamlog_out(DEBUG5) << endl;
    streamlog_out(DEBUG5) << " icn " << iicn << " has type " << type_icn_parent(iicn) << endl;
    streamlog_out(DEBUG5) << " True :         " << E_icn(iicn);
    for (int jj = 0; jj < 3; jj++) {
      streamlog_out(DEBUG5) << " " << p_icn(iicn)[jj];
    };
    streamlog_out(DEBUG5) << " " << M_icn(iicn) << endl;
    streamlog_out(DEBUG4) << "               ";
    for (int jj = 0; jj < 4; jj++) {
      streamlog_out(DEBUG4) << " " << p4_icn(iicn)[jj];
    };
    streamlog_out(DEBUG4) << endl;

    // get and list the jets, once again as an IntVec with inicies into the ets ReconstructedParticleVec.
    streamlog_out(DEBUG5) << " Number of jets from this cn:" << jets_of_initial_cn(iicn).size() << endl;
    if (jets_of_initial_cn(iicn).size() > 0) {
      if (jets_of_initial_cn(iicn).size() > 1) {
        streamlog_out(DEBUG5) << " jets are :";
      } else {
        streamlog_out(DEBUG5) << " jet is :  ";
      }
      for (unsigned jj = 0; jj < jets_of_initial_cn(iicn).size(); ++jj) {
        streamlog_out(DEBUG5) << " " << jets_of_initial_cn(iicn)[jj];
      };
      streamlog_out(DEBUG5) << endl;
    }

    // sum up energies of all true particles contributing to each true jet from this colour neutral.
    IntVec icn_jets = jets_of_initial_cn(iicn);
    double p4_fr_jets[4] = {0, 0, 0, 0};
    for (unsigned kk = 0; kk < icn_jets.size(); kk++) {
      for (int jj = 0; jj < 4; jj++) {
        p4_fr_jets[jj] += p4true(icn_jets[kk])[jj];
      }
    }

    // Consiency check: is the sum-of-true and the E/P of the colour neutral the same ?
    double M = sqrt(p4_fr_jets[0] * p4_fr_jets[0] - p4_fr_jets[1] * p4_fr_jets[1] - p4_fr_jets[2] * p4_fr_jets[2] -
                    p4_fr_jets[3] * p4_fr_jets[3]);
    streamlog_out(DEBUG5) << " true p4 and mass of these: ";
    for (int jj = 0; jj < 4; jj++) {
      streamlog_out(DEBUG5) << " " << p4_fr_jets[jj];
    };
    streamlog_out(DEBUG5) << " " << M << endl;
    // Normally, the mass should be identical. However, due to the B-field, decyas of longlived
    // charged particles have their momentum *direction* changed, so check energy as well
    float margin = 0.007;
    // If this is a tau cn, be more forgiving. The theory is that there is an issue with the
    // the crossing-angle boost that becomes too big in aa/BB events which effects tau:s most
    if (abs(pdg_icn_comps(iicn)[0]) == 15 || abs(pdg_icn_comps(iicn)[0]) == 16) {
      margin = 0.02;
    }
    if (abs(p4_fr_jets[0] - E_icn(iicn)) / E_icn(iicn) > margin) {
      if (((!cluster_in_event) && (!top_in_event)) ||
          (abs((total_Ecn - total_Eisr) - total_Etrue) / total_Etrue > 0.001)) {
        streamlog_out(WARNING) << " event/run " << evt->getEventNumber() << "/" << evt->getRunNumber()
                               << " , initial cn " << iicn << " E of cn:    " << E_icn(iicn)
                               << " E of jets :    " << p4_fr_jets[0] << endl;
        streamlog_out(WARNING) << " event/run " << evt->getEventNumber() << "/" << evt->getRunNumber()
                               << " , initial cn " << iicn << " mass of cn: " << M_icn(iicn) << " mass of jets : " << M
                               << endl;
        streamlog_out(WARNING) << " event/run " << " cluster_in_event : " << cluster_in_event
                               << " gluonsplit_in_event : " << gluonsplit_in_event
                               << "  top_in_event : " << top_in_event << endl;
        streamlog_out(WARNING) << " event/run " << " Totals : all non-isr initials = " << total_Ecn - total_Eisr
                               << " all final non-ISR " << total_Etrue << endl;
      }
    }

    // get the list (MCParticleVec) of the elementons making the colour-neutral, and print their PDG.
    streamlog_out(DEBUG5) << " Number of elementons making this cn:" << elementons_initial_cn(iicn).size() << endl;
    streamlog_out(DEBUG5) << " Pdg of the parent " << pdg_icn_parent(iicn) << endl;
    streamlog_out(DEBUG5) << " Pdgs/jet/parent jet/type  of elementons making this cn: " << endl;
    for (unsigned jj = 0; jj < pdg_icn_comps(iicn).size(); jj++) {
      if (type_icn_comps(iicn)[jj] < 100) {
        streamlog_out(DEBUG5) << setw(10) << pdg_icn_comps(iicn)[jj] << setw(4) << jets_of_initial_cn(iicn)[jj]
                              << setw(4) << 0 << setw(4) << type_icn_comps(iicn)[jj] % 100 << endl;
      } else {
        streamlog_out(DEBUG5) << setw(10) << pdg_icn_comps(iicn)[jj] << setw(4) << jets_of_initial_cn(iicn)[jj]
                              << setw(4) << type_icn_comps(iicn)[jj] / 100 - 1 << setw(4)
                              << type_icn_comps(iicn)[jj] % 100 << endl;
      }
    };
    streamlog_out(DEBUG5) << endl;
  }
  // consitency check: Did the sum of MCParticles really add up to the initial energy ?!
  streamlog_out(DEBUG5) << " Total energies : " << " E(physics event)= " << total_Ecn
                        << " E(physics event) - E(isr)= " << total_Ecn - total_Eisr
                        << " E(sum non-isr true parts)= " << total_Etrue << " E(isr)= " << total_Eisr
                        << " E(overlay,seen)= " << total_Eowl << endl;
  if (abs(total_Ecn - total_Eisr - total_Etrue) / total_Etrue > 0.009) {
    streamlog_out(WARNING) << " event/run " << evt->getEventNumber() << "/" << evt->getRunNumber()
                           << ": Mis-match in energy between initial elementons and sum of true particles: "
                           << " E(physics event)= " << total_Ecn - total_Eisr << " E(sum true parts)= " << total_Etrue
                           << endl;
  }

  // Now the other way round: Get the PFOs or MCPs and check which jet each one belongs to:
  LCCollection* pfocol = NULL;
  try {
    pfocol = evt->getCollection(_recoParticleCollectionName);
  } catch (lcio::DataNotAvailableException& e) {
    // streamlog_out(WARNING) <<  _recoParticleCollectionName   << " collection not available" << std::endl;
    pfocol = NULL;
  }

  if (pfocol != NULL) {
    int nPFO = pfocol->getNumberOfElements();
    streamlog_out(DEBUG4) << " list of pfos and their true jets  " << endl;
    // Loop the PFOs : the complicated way showing use of the navigators
    for (int j = 0; j < nPFO; j++) {
      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(pfocol->getElementAt(j));
      LCObjectVec jetvec = reltjreco->getRelatedFromObjects(pfo);
      // jet the PFO belongs to, as a ReconstructedParticle*. jetindex() gives the index in the jets
      // ReconstructedParticleVec, which is the language
      // of all the printing above.
      streamlog_out(DEBUG3) << pfo << " ";
      for (unsigned kk = 0; kk < jetvec.size(); kk++) {
        streamlog_out(DEBUG3) << " " << jetindex(dynamic_cast<ReconstructedParticle*>(jetvec[kk]));
      }
      streamlog_out(DEBUG3) << endl;
    }
    // Same thing, the easy way. In addition, shows the initial and final cn:s of each PFO
    for (int j = 0; j < nPFO; j++) {
      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(pfocol->getElementAt(j));
      streamlog_out(DEBUG4) << pfo << " " << recojet(pfo) << " " << recoicn(pfo) << " " << recofcn(pfo) << std::endl;
    }
  }
  LCCollection* mcpcol = NULL;
  try {
    mcpcol = evt->getCollection(_MCParticleColllectionName);
  } catch (lcio::DataNotAvailableException& e) {
    mcpcol = NULL;
  }

  if (mcpcol != NULL) {
    int nMCP = mcpcol->getNumberOfElements();
    streamlog_out(DEBUG4) << " list of mcps and their true jets  " << endl;
    // Loop the MCPs : the complicated way showing use of the navigators
    for (int j = 0; j < nMCP; j++) {
      MCParticle* mcp = dynamic_cast<MCParticle*>(mcpcol->getElementAt(j));
      LCObjectVec jetvec = reltjmcp->getRelatedFromObjects(mcp);
      // jet the MCP belongs to, as a ReconstructedParticle*. jetindex() gives the index in the jets
      // ReconstructedParticleVec, which is the language
      // of all the printing above.
      streamlog_out(DEBUG3) << mcp << " ";
      for (unsigned kk = 0; kk < jetvec.size(); kk++) {
        streamlog_out(DEBUG3) << " " << jetindex(dynamic_cast<ReconstructedParticle*>(jetvec[kk]));
      }
      streamlog_out(DEBUG3) << endl;
    }
    // Same thing, the easy way. In addition, shows the initial and final cn:s of each MCP
    for (int j = 0; j < nMCP; j++) {
      MCParticle* mcp = dynamic_cast<MCParticle*>(mcpcol->getElementAt(j));
      if (mcpjet(mcp) > -1000) { // -1000 for mcp not included in any jet - the case for e.g. the incomming particles
        streamlog_out(DEBUG4) << mcp << " " << mcpjet(mcp) << " " << mcpicn(mcp) << " " << mcpfcn(mcp) << std::endl;
      }
    }
  }

  if (tj) {
    delall(); // clean up
  }
  _nEvt++;
}

void Use_TrueJet::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void Use_TrueJet::end() {

  //   std::cout << "TrueJet::end()  " << name()
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;
}
