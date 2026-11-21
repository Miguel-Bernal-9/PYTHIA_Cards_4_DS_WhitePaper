#include "Pythia8/Pythia.h"

// Preferably use HepMC3, but alternatively HepMC2.
#ifndef HEPMC2
#include "Pythia8Plugins/HepMC3.h"
#else
#include "Pythia8Plugins/HepMC2.h"
#endif

#include <fstream>

using namespace Pythia8;

// =====================================================================
// CARD SUMMARY (read this first)
// ---------------------------------------------------------------------
//  Output name       : Nc3Nf8_sFOn_pp_3pi_decay_lam_10_probvec_75_v8315
//                      (.hepmc and _HVparameters.txt in outDir)
//  Pythia version    : 8.315. It needs the 8.313+ settings (setLambda,
//                      setabsigma, bLund, sigmaLund, rFact): with 8.307 init()
//                      fails. In 8.315 the ratio m(rho_v)/m(rho0) rescales the
//                      HV stop mass AND other HV string parameters (mJoin,
//                      tolerances); for 8.313 use the 8.313 cards
//  Process           : pp, 13 TeV -> Z' (m = 1 TeV, width 20 GeV, 800-1200 GeV)
//                      -> qv qvbar. The Z' is forced to decay to dark quarks only
//  Dark sector       : SU(3) (Nc = 3), Nf = 8, separateFlav = on
//                      Lambda = 10 GeV (setLambda = on), alphaOrder = 1,
//                      probVector = 0.75
//                      m(qv) = 10.2 GeV, m(pi_v) = 6 GeV, m(rho_v) = 25 GeV
//                      (the same mass for every flavour)
//  Hadronization     : setabsigma = 2, aLund = 0.3, bLund = 0.0077 (= 0.8/10.2^2),
//                      sigmaLund = 5.1 (= 0.5*10.2), rFact = 1: the same values as
//                      the 8.307 defaults (aLund, bmqv2, sigmamqv, rFactqv)
//  Dark shower cutoff: pTminFSR is raised internally to 1.1*Lambda = 11 GeV
//                      (the warning 'Hidden Valley pTmin too low' is expected)
//  HV stop mass      : 32.24 GeV, the 8.315 DEFAULT: StringFragmentation:stopMass
//                      (1.0, not changed) * m(rho_v)/m(rho0) = 25/0.77549.
//                      This is NOT the 8.307 value (9 GeV): fewer and harder
//                      dark pions than in 8.307. For the 8.307 behaviour use
//                      the _stopmass9 card
//  Dark hadron decays: ALL OFF. Every dark pion and dark rho is stable and is in
//                      the final state (status 1). The '3pi_decay' in the name
//                      only refers to the original Snowmass benchmark: here NO
//                      dark hadron decays, so there are no SM particles from the
//                      dark sector and r_inv is NOT defined for these samples.
//                      Use them only for generator-level dark-hadron studies
//  PDG codes written : separateFlav = on -> 4900ijs with 8 >= i >= j >= 1
//                      (s = 1 pi_v, s = 3 rho_v): 36 pion + 36 rho labels
//  Random seed       : Random:seed = 0 -> taken from the clock, so every run gives
//                      different events. Use a fixed seed to reproduce a sample
//  Events            : 50000 (about 10 GB of HepMC3 ASCII)
//  WARNING           : running the card again OVERWRITES the .hepmc and the
//                      _HVparameters.txt files with the same outName
//  Validation        : 1000 events, seed 12345, Pythia 8.315 (final state, per event):
//                      <N(pi_v)> = 3.21, <N(rho_v)> = 8.93,
//                      <pT(pi_v)> = 32.7 GeV, <pT(rho_v)> = 67.6 GeV
//                      (statistical uncertainty about 1-2% on each number)
// ---------------------------------------------------------------------
// Tags used in the comments of the code:
//
// [SFon] changes with respect to the separateFlav = off card. With
// separateFlav = on every dark quark (4900101 ... 4900108) and every dark
// meson 4900ij1 / 4900ij3 (8 >= i >= j >= 1) has its own particle data, so all
// of them must be set; otherwise Pythia uses its defaults (qv: 5 GeV, all
// mesons: 10 GeV). The dark Z' BR (0.4) is split as 0.4/8 per flavour.
//
// [8.313 fix] NOT needed in 8.315: the 8.313 bug (the lightest vector meson
// was searched among all 4900ij3 even with separateFlav = off) was fixed in
// 8.314. In 8.315 only 4900113 is used with separateFlav = off, and all
// 4900ij3 with i, j <= Nf with separateFlav = on (all set to 25 GeV here).
// The printout 'lightest vector used by Pythia' must say 25 GeV.
//
// [stable] all dark hadrons are stable (mayDecay = off for every species).
// A counter in the event loop checks that no dark meson decayed and prints
// '[stable] ... (all stable, OK)' at the end of the log.
//
// Note on r_inv (only for the ORIGINAL benchmark, with decays on):
// r_inv ~ 1 - 1/Nf * Br(pi_0 -> visible), because only the Nf diagonal pions
// can decay to the SM. It does not apply to this card (no decays).
// =====================================================================

// --- Output files ---

const string outDir  = "/scratch/maab1g24/Datasets_IRIDIS_SVJs/HEPMC_Files/DS_TF_Files/";
const string outName = "Nc3Nf8_sFOn_pp_3pi_decay_lam_10_probvec_75_v8315";

// --- [SFon] Number of dark flavours and dark sector masses (used in the loops below) ---

const int    NFLAV = 8;
const double MASS_QV  = 10.2;   // dark quarks
const double MASS_PIV  = 6.0;    // dark pions (spin 0)
const double MASS_RHOV = 25.0;   // dark rhos  (spin 1)

// --- [SFon] PDG codes of the dark mesons 4900ijs, i >= j (s = 1 pion, 3 rho) ---

int darkMeson(int i, int j, int s) { return 4900000 + 100 * i + 10 * j + s; }

//==========================================================================

// --- Print all the HiddenValley settings and the dark particle data ---
// --- (values used in the generation, read after pythia.init()) ---

void printHVParameters(Pythia& pythia, ostream& os) {

  os << "\n =========== Hidden Valley settings used in the generation ===========\n";
  os << " Pythia version: " << PYTHIA_VERSION << "\n\n";

  for (auto& f : pythia.settings.getFlagMap("HiddenValley:"))
    os << "  " << setw(34) << left << f.second.name << " = "
       << (f.second.valNow ? "on" : "off")
       << "   (default " << (f.second.valDefault ? "on" : "off") << ")\n";
  for (auto& m : pythia.settings.getModeMap("HiddenValley:"))
    os << "  " << setw(34) << left << m.second.name << " = " << m.second.valNow
       << "   (default " << m.second.valDefault << ")\n";
  for (auto& p : pythia.settings.getParmMap("HiddenValley:"))
    os << "  " << setw(34) << left << p.second.name << " = " << p.second.valNow
       << "   (default " << p.second.valDefault << ")\n";
  for (auto& w : pythia.settings.getWordMap("HiddenValley:"))
    os << "  " << setw(34) << left << w.second.name << " = " << w.second.valNow << "\n";

  // --- pTminFSR actually used in the dark shower: for alphaOrder = 1 Pythia ---
  // --- raises it internally to 1.1*Lambda if it is lower (see the warning ---
  // --- "Hidden Valley pTmin too low"), but the setting above is not updated ---

  double lambdaHV = pythia.settings.parm("HiddenValley:Lambda");
  double pTminHV  = pythia.settings.parm("HiddenValley:pTminFSR");
  if (pythia.settings.mode("HiddenValley:alphaOrder") == 1)
    os << "\n  Effective dark shower cutoff = max(pTminFSR, 1.1*Lambda) = "
       << max(pTminHV, 1.1 * lambdaHV) << " GeV\n";

  // --- Lightest vector meson and HV string stop mass actually used (8.313+): ---
  // --- same scan as HiddenValleyFragmentation::init(): 4900113 plus 4900ij3 ---
  // --- for i = 2..nFlav, j <= i; in 8.313 even with separateFlav = off (bug) ---
  int nFlavHV = pythia.settings.mode("HiddenValley:nFlav");
  double mVec = pythia.particleData.m0(4900113);
  if (pythia.settings.flag("HiddenValley:separateFlav") || PYTHIA_VERSION_INTEGER < 8314)
    for (int i = 2; i <= nFlavHV; ++i)
      for (int j = 1; j <= i; ++j)
        mVec = min(mVec, pythia.particleData.m0(4900003 + 100 * i + 10 * j));
  double kHV = mVec / pythia.particleData.m0(113);
  os << "  lightest vector used by Pythia  = " << mVec << " GeV\n"
     << "  k = m(rho_v)/m(rho0)            = " << kHV << "\n"
     << "  StringFragmentation:stopMass   = "
     << pythia.settings.parm("StringFragmentation:stopMass") << "\n"
     << "  HV stop mass = stopMass * k    = "
     << pythia.settings.parm("StringFragmentation:stopMass") * kHV
     << " GeV   (8.307: 1.5 * m(pi_v) = " << 1.5 * pythia.particleData.m0(4900111) << " GeV)\n";

  // --- [SFon] check that all masses are the ones requested ---

  double piMin = 1e9, piMax = 0., rhoMin = 1e9, rhoMax = 0., qMin = 1e9, qMax = 0.;
  for (int i = 1; i <= NFLAV; ++i) {
    qMin = min(qMin, pythia.particleData.m0(4900100 + i));
    qMax = max(qMax, pythia.particleData.m0(4900100 + i));
    for (int j = 1; j <= i; ++j) {
      piMin  = min(piMin,  pythia.particleData.m0(darkMeson(i, j, 1)));
      piMax  = max(piMax,  pythia.particleData.m0(darkMeson(i, j, 1)));
      rhoMin = min(rhoMin, pythia.particleData.m0(darkMeson(i, j, 3)));
      rhoMax = max(rhoMax, pythia.particleData.m0(darkMeson(i, j, 3)));
    }
  }
  int nSpecies = NFLAV * (NFLAV + 1) / 2;
  os << "  dark quark masses (" << NFLAV << ")         = " << qMin   << " - " << qMax   << " GeV\n"
     << "  dark pion  masses (" << nSpecies << ")        = " << piMin  << " - " << piMax  << " GeV\n"
     << "  dark rho   masses (" << nSpecies << ")        = " << rhoMin << " - " << rhoMax << " GeV\n";

  // --- Particle data: Z', dark quarks, all dark mesons, invisible state ---

  os << "\n =========== Dark sector particle data ===========\n";
  os << "  " << setw(10) << left << "id" << setw(12) << "name" << setw(12) << "m0"
     << setw(12) << "mWidth" << setw(12) << "mMin" << setw(12) << "mMax"
     << "mayDecay\n";
  vector<int> ids = {4900023};
  for (int i = 1; i <= NFLAV; ++i) ids.push_back(4900100 + i);
  for (int s : {1, 3})
    for (int i = 1; i <= NFLAV; ++i)
      for (int j = 1; j <= i; ++j) ids.push_back(darkMeson(i, j, s));
  ids.push_back(53);
  for (int id : ids) {
    ParticleDataEntryPtr pd = pythia.particleData.findParticle(id);
    if (!pd) continue;
    os << "  " << setw(10) << left << id << setw(12) << pd->name()
       << setw(12) << pd->m0() << setw(12) << pd->mWidth()
       << setw(12) << pd->mMin() << setw(12) << pd->mMax()
       << (pd->mayDecay() ? "on" : "off") << "\n";
    for (int i = 0; i < pd->sizeChannels(); ++i) {
      DecayChannel& ch = pd->channel(i);
      os << "      channel " << i << ": onMode = " << ch.onMode()
         << ", bRatio = " << ch.bRatio() << ", meMode = " << ch.meMode()
         << ", products =";
      for (int j = 0; j < ch.multiplicity(); ++j) os << " " << ch.product(j);
      os << "\n";
    }
  }
  os << " ====================================================================\n\n";
}

//==========================================================================

int main() {

  // Interface for conversion from Pythia8::Event to HepMC event.
  // Specify file where HepMC events will be stored.
  Pythia8ToHepMC toHepMC(outDir + outName + ".hepmc");

  // --- Stop if the HepMC file could not be created ---
  {
    ifstream test(outDir + outName + ".hepmc");
    if (!test.good()) {
      cout << " ERROR: could not create " << outDir + outName + ".hepmc" << endl;
      return 1;
    }
  }

  // Generator.
  Pythia pythia;

  // --- 1) Settings that will be used in a main program ---

  pythia.readString("Main:numberOfEvents = 50000");     // number of events to generate
  pythia.readString("Main:timesAllowErrors = 3");       // abort run after this many flawed events

  // --- Random seed = 0 gives random seed based on time ---

  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");

  // --- 2) Settings related to output in init(), next() and stat() ---

  pythia.readString("Init:showChangedSettings = on");      // list changed settings
  pythia.readString("Init:showAllSettings = on");          // list all settings
  pythia.readString("Init:showChangedParticleData = on");  // list changed particle data
  pythia.readString("Init:showAllParticleData = on");      // list all particle data
  pythia.readString("Next:numberCount = 1000");            // print message every n events
  pythia.readString("Next:numberShowLHA = 1");             // print LHA information n times
  pythia.readString("Next:numberShowInfo = 1");            // print event information n times
  pythia.readString("Next:numberShowProcess = 1");         // print process record n times
  pythia.readString("Next:numberShowEvent = 1");           // print event record n times
  pythia.readString("Stat:showPartonLevel = on");          // additional statistics on MPI

  // --- 3) Beam parameter settings. Values below agree with default ones ---

  pythia.readString("Beams:idA = 2212");                   // first beam, p = 2212, pbar = -2212
  pythia.readString("Beams:idB = 2212");                   // second beam, p = 2212, pbar = -2212
  pythia.readString("Beams:eCM = 13000.");                 // CM energy of collision

  // --- 4) Turn on the production process ---

  pythia.readString("HiddenValley:ffbar2Zv = on");
  pythia.readString("4900023:m0 = 1000");
  pythia.readString("4900023:mMin = 800");
  pythia.readString("4900023:mMax = 1200");
  pythia.readString("4900023:mWidth = 20.0");
  pythia.readString("4900023:doForceWidth = on");

  pythia.readString("4900023:oneChannel = 1 0.0750000 100 1 -1");
  pythia.readString("4900023:addChannel = 1 0.0750000 100 2 -2");
  pythia.readString("4900023:addChannel = 1 0.0750000 100 3 -3");
  pythia.readString("4900023:addChannel = 1 0.0750000 100 4 -4");
  pythia.readString("4900023:addChannel = 1 0.0750000 100 5 -5");
  pythia.readString("4900023:addChannel = 1 0.0750000 100 6 -6");
  pythia.readString("4900023:addChannel = 1 0.0250000 100 11 -11");
  pythia.readString("4900023:addChannel = 1 0.0250000 100 12 -12");
  pythia.readString("4900023:addChannel = 1 0.0250000 100 13 -13");
  pythia.readString("4900023:addChannel = 1 0.0250000 100 14 -14");
  pythia.readString("4900023:addChannel = 1 0.0250000 100 15 -15");
  pythia.readString("4900023:addChannel = 1 0.0250000 100 16 -16");

  // [SFon] With separateFlav = off, 4900101 stands for all Nf flavours. With
  // separateFlav = on it is only flavour 1, so the dark BR (0.4) is split
  // equally among the 8 dark quarks (0.4 / 8 = 0.05 each), as the Zv does by
  // default in Pythia.

  pythia.readString("4900023:addChannel = 1 0.05 100 4900101 -4900101");
  pythia.readString("4900023:addChannel = 1 0.05 100 4900102 -4900102");
  pythia.readString("4900023:addChannel = 1 0.05 100 4900103 -4900103");
  pythia.readString("4900023:addChannel = 1 0.05 100 4900104 -4900104");
  pythia.readString("4900023:addChannel = 1 0.05 100 4900105 -4900105");
  pythia.readString("4900023:addChannel = 1 0.05 100 4900106 -4900106");
  pythia.readString("4900023:addChannel = 1 0.05 100 4900107 -4900107");
  pythia.readString("4900023:addChannel = 1 0.05 100 4900108 -4900108");

  // --- 5) For MC efficiency let zprime decay only to dark quarks ---

  pythia.readString("4900023:onMode = off");
  pythia.readString("4900023:onIfAny = 4900101 4900102 4900103 4900104 "
                    "4900105 4900106 4900107 4900108");                   // [SFon]

  // --- 6) Decouple default bifundamental quarks ---

  pythia.readString("4900001:m0 = 5000");
  pythia.readString("4900002:m0 = 5000");
  pythia.readString("4900003:m0 = 5000");
  pythia.readString("4900004:m0 = 5000");
  pythia.readString("4900005:m0 = 5000");
  pythia.readString("4900006:m0 = 5000");
  pythia.readString("4900011:m0 = 5000");
  pythia.readString("4900012:m0 = 5000");
  pythia.readString("4900013:m0 = 5000");
  pythia.readString("4900014:m0 = 5000");
  pythia.readString("4900015:m0 = 5000");
  pythia.readString("4900016:m0 = 5000");

  // --- 7) Set mass of dark quarks ---
  // --- this is Lambda + current quark mass ---
  // --- [SFon] one mass per flavour, 4900101 ... 4900108 (default would be 5 GeV) ---

  for (int i = 1; i <= NFLAV; ++i)
    pythia.readString(to_string(4900100 + i) + ":m0 = " + to_string(MASS_QV));

  // --- 8) Set gauge group and number of flavours ---
  // --- caution: user should not use Nc=2 or Nf=1 without consulting theory experts ---

  pythia.readString("HiddenValley:Ngauge  = 3");
  pythia.readString("HiddenValley:Nflav  = 8");
  pythia.readString("HiddenValley:separateFlav  = on");                  // [SFon]

  // --- 9) Set other darksector parameters, must set HiddenValley:alphaOrder = 1.0 ---

  pythia.readString("HiddenValley:setLambda = on");     // 8.310+: else Lambda is ignored
  pythia.readString("HiddenValley:Lambda = 10");
  pythia.readString("HiddenValley:FSR = on");
  pythia.readString("HiddenValley:fragment = on");
  pythia.readString("HiddenValley:alphaOrder = 1");
  pythia.readString("HiddenValley:probVector = 0.75");
  pythia.readString("HiddenValley:setabsigma = 2");
  pythia.readString("HiddenValley:aLund = 0.3");
  pythia.readString("HiddenValley:bLund = 0.0077");
  pythia.readString("HiddenValley:sigmaLund = 5.1");
  pythia.readString("HiddenValley:rFact = {1.,1.,1.,1.,1.,1.,1.,1.}");

  // --- 10) Set masses of spin-0 and spin-1 mesons ---
  // --- [SFon] every species 4900ij1 (pion) and 4900ij3 (rho), 8 >= i >= j >= 1: ---
  // --- 8 diagonal + 28 off-diagonal = 36 of each (default would be 10 GeV). ---
  // --- All vector mesons have the same mass (needed in 8.313, harmless in 8.315). ---
  // This was set to 26 in dark showers paper v1 (rho)

  for (int i = 1; i <= NFLAV; ++i)
    for (int j = 1; j <= i; ++j) {
      pythia.readString(to_string(darkMeson(i, j, 1)) + ":m0 = " + to_string(MASS_PIV));
      pythia.readString(to_string(darkMeson(i, j, 3)) + ":m0 = " + to_string(MASS_RHOV));
    }

  // --- Declare additional stable state to emulate invisible fraction ---

  pythia.readString("53:m0 = 0.0");
  pythia.readString("53:isResonance = false");
  pythia.readString("53:mayDecay = off");

  // --- [stable] All dark mesons stable, so the ---
  // --- dark hadron spectrum can be studied at generator level ---
  // --- [SFon] all 72 species (the decay channels stay commented, as in the ---
  // --- original card: with mayDecay = off they would have no effect anyway) ---

  for (int i = 1; i <= NFLAV; ++i)
    for (int j = 1; j <= i; ++j) {
      pythia.readString(to_string(darkMeson(i, j, 1)) + ":mayDecay = off");
      pythia.readString(to_string(darkMeson(i, j, 3)) + ":mayDecay = off");
    }
  // pythia.readString("4900113:oneChannel = 1 1 100 4900211 -4900211");
  // pythia.readString("4900213:oneChannel = 1 1 100 4900111 4900211");
  // pythia.readString("4900111:onechannel = 1  1  91 -4 4");
  // pythia.readString("4900111:addchannel = 1  0 100 53 -53");

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // --- Print the HV parameters used in the generation: on screen ---
  // --- and in a text file next to the HepMC file ---

  printHVParameters(pythia, cout);
  ofstream parFile(outDir + outName + "_HVparameters.txt");
  printHVParameters(pythia, parFile);
  parFile.close();

  // Begin event loop. Generate event. Skip if error.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");
  int iAbort = 0;

  // --- [stable] counters: dark mesons produced and dark mesons that decayed ---
  long nDarkMesons = 0, nDarkDecayed = 0;

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // --- [stable] dark mesons 4900ijs (i, j >= 1, s = 1 or 3) must be final ---
    for (int i = 0; i < pythia.event.size(); ++i) {
      int id = pythia.event[i].idAbs();
      int s = id % 10, j = (id / 10) % 10, k = (id / 100) % 10;
      if (id / 1000 == 4900 && k >= 1 && j >= 1 && (s == 1 || s == 3)) {
        ++nDarkMesons;
        if (!pythia.event[i].isFinal()) ++nDarkDecayed;
      }
    }

    // Construct new empty HepMC event, fill it and write it out.
    toHepMC.writeNextEvent( pythia );

  // End of event loop. Statistics.
  }
  pythia.stat();

  // --- [stable] summary of the stability check ---
  cout << "\n [stable] dark mesons produced: " << nDarkMesons
       << ", dark mesons that decayed: " << nDarkDecayed
       << (nDarkDecayed == 0 ? "  (all stable, OK)" : "  (WARNING: some decayed!)") << "\n";

  // Done.
  return 0;
}
