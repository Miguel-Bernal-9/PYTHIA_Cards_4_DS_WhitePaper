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
//  Output name       : Nc3Nf8_sFoff_pp_3pi_decay_lam_10_probvec_75_v8315_stopmass9
//                      (.hepmc and _HVparameters.txt in outDir)
//  Pythia version    : 8.315. It needs the 8.313+ settings (setLambda,
//                      setabsigma, bLund, sigmaLund, rFact): with 8.307 init()
//                      fails. In 8.315 the ratio m(rho_v)/m(rho0) rescales the
//                      HV stop mass AND other HV string parameters (mJoin,
//                      tolerances); for 8.313 use the 8.313 cards
//  Process           : pp, 13 TeV -> Z' (m = 1 TeV, width 20 GeV, 800-1200 GeV)
//                      -> qv qvbar. The Z' is forced to decay to dark quarks only
//  Dark sector       : SU(3) (Nc = 3), Nf = 8, separateFlav = off
//                      Lambda = 10 GeV (setLambda = on), alphaOrder = 1,
//                      probVector = 0.75
//                      m(qv) = 10.2 GeV, m(pi_v) = 6 GeV, m(rho_v) = 25 GeV
//                      (the same mass for every flavour)
//  Hadronization     : setabsigma = 2, aLund = 0.3, bLund = 0.0077 (= 0.8/10.2^2),
//                      sigmaLund = 5.1 (= 0.5*10.2), rFact = 1: the same values as
//                      the 8.307 defaults (aLund, bmqv2, sigmamqv, rFactqv)
//  Dark shower cutoff: pTminFSR is raised internally to 1.1*Lambda = 11 GeV
//                      (the warning 'Hidden Valley pTmin too low' is expected)
//  HV stop mass      : 9 GeV = 1.5 * m(pi_v), the 8.307 value, FORCED through
//                      StringFragmentation:stopMass = 9 * 0.77549/25 = 0.27918
//                      (see [stop9]). Only the stop mass changes: mJoin and the
//                      tolerances are still rescaled with 25/0.77549.
//                      WARNING: this parameter is shared with
//                      the SM strings, so the SM hadronization (beam remnants,
//                      ISR, MPI) is also changed. Valid only for dark-sector
//                      (generator-level) studies, NOT for SM-hadron observables
//  Dark hadron decays: ALL OFF. Every dark pion and dark rho is stable and is in
//                      the final state (status 1). The '3pi_decay' in the name
//                      only refers to the original Snowmass benchmark: here NO
//                      dark hadron decays, so there are no SM particles from the
//                      dark sector and r_inv is NOT defined for these samples.
//                      Use them only for generator-level dark-hadron studies
//  PDG codes written : separateFlav = off -> only 4 labels (plus antiparticles):
//                      4900111, 4900211 (pi_v), 4900113, 4900213 (rho_v).
//                      Pythia still tracks the flavour internally: 111/113 =
//                      flavour diagonal, 211/213 = off-diagonal
//  Random seed       : Random:seed = 0 -> taken from the clock, so every run gives
//                      different events. Use a fixed seed to reproduce a sample
//  Events            : 50000 (about 10 GB of HepMC3 ASCII)
//  WARNING           : running the card again OVERWRITES the .hepmc and the
//                      _HVparameters.txt files with the same outName
//  Validation        : 1000 events, seed 12345, Pythia 8.315 (final state, per event):
//                      <N(pi_v)> = 3.89, <N(rho_v)> = 9.11,
//                      <pT(pi_v)> = 28.7 GeV, <pT(rho_v)> = 64.1 GeV
//                      (statistical uncertainty about 1-2% on each number)
// ---------------------------------------------------------------------
// Tags used in the comments of the code:
//
// [8.313 fix] NOT needed in 8.315: the 8.313 bug (the lightest vector meson
// was searched among all 4900ij3 even with separateFlav = off) was fixed in
// 8.314. In 8.315 only 4900113 is used with separateFlav = off, and all
// 4900ij3 with i, j <= Nf with separateFlav = on (all set to 25 GeV here).
// The loop over the 36 vector mesons 4900ij3 is still present: harmless,
// kept so that the 8.313 and 8.315 cards use identical particle data.
// The printout 'lightest vector used by Pythia' must say 25 GeV.
//
// [stop9] StringFragmentation:stopMass = 1.5 * m(pi_v) * m(rho0)/m(rho_v)
// = 0.27918, so that 8.315 uses an HV stop mass of 9 GeV, as 8.307 does.
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
const string outName = "Nc3Nf8_sFoff_pp_3pi_decay_lam_10_probvec_75_v8315_stopmass9";

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

  os << "\n =========== Dark sector particle data ===========\n";
  os << "  " << setw(10) << left << "id" << setw(12) << "name" << setw(12) << "m0"
     << setw(12) << "mWidth" << setw(12) << "mMin" << setw(12) << "mMax"
     << "mayDecay\n";
  vector<int> ids = {4900023, 4900101, 4900111, 4900211, 4900113, 4900213, 53};
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

  // --- Stop if the HepMC file could not be created (wrong outDir, no quota...) ---
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
  pythia.readString("4900023:addChannel = 1 0.4 100 4900101 -4900101");

  // --- 5) For MC efficiency let zprime decay only to dark quarks ---

  pythia.readString("4900023:onMode = off");
  pythia.readString("4900023:onIfAny = 4900101");

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

  // --- 7) Set mass of dark quark ---
  // --- this is Lambda + current quark mass ---

  pythia.readString("4900101:m0 = 10.2");

  // --- 8) Set gauge group and number of flavours ---
  // --- caution: user should not use Nc=2 or Nf=1 without consulting theory experts ---

  pythia.readString("HiddenValley:Ngauge  = 3");
  pythia.readString("HiddenValley:Nflav  = 8");
  pythia.readString("HiddenValley:separateFlav  = off");

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

  // --- 10) Set masses and branching ratios of spin-1 and spin-0 mesons ---
  // --- For now this treats flavor singlets on the same footing as flavor ---
  // --- adjoint multiplets, this is not entirely consistent. ---
  // --- adjustment to flavor singlet is possible but only with separateFlav=on ---

  pythia.readString("4900111:m0 = 6.0");
  pythia.readString("4900211:m0 = 6.0");
  // This was set to 26 in dark showers paper v1
  pythia.readString("4900113:m0 = 25");
  pythia.readString("4900213:m0 = 25");

  // --- [8.313 fix] same mass for ALL vector mesons 4900ij3 (8 >= i >= j >= 1): ---
  // --- not needed in 8.315, kept for identical particle data with 8.313 ---
  for (int i = 1; i <= 8; ++i)
    for (int j = 1; j <= i; ++j)
      pythia.readString(to_string(4900003 + 100 * i + 10 * j) + ":m0 = 25");

  // --- [stop9] HV string stop mass = 9 GeV, as in 8.307 (1.5 * m(pi_v)). ---
  // --- From 8.313 on Pythia uses stopMass(HV) = StringFragmentation:stopMass ---
  // --- * m(rho_v,min)/m(rho0), so the QCD parameter is set to ---
  // --- 1.5 * m(pi_v) * m(rho0) / m(rho_v) = 9 * 0.77549 / 25 = 0.2792. ---
  // --- NOTE: this parameter is shared with the SM strings, so the SM ---
  // --- hadronization (beam remnants, ISR, MPI) also changes. Fine for ---
  // --- dark-sector-only studies, not for observables with SM hadrons. ---
  // --- Must come after the dark meson masses are set (they are used here). ---

  {
    double stopMassQCD = 1.5 * pythia.particleData.m0(4900111)
                       * pythia.particleData.m0(113) / pythia.particleData.m0(4900113);
    pythia.readString("StringFragmentation:stopMass = " + to_string(stopMassQCD));
  }

  // --- Declare additional stable state to emulate invisible fraction ---

  pythia.readString("53:m0 = 0.0");
  pythia.readString("53:isResonance = false");
  pythia.readString("53:mayDecay = off");

  // --- [stable] Spin-1 mesons stable. The decay channels of the original ---
  // --- benchmark (rho_v -> pi_v pi_v) are commented out and not used. ---

  pythia.readString("4900113:mayDecay = off");
  // pythia.readString("4900113:oneChannel = 1 1 100 4900211 -4900211");

  pythia.readString("4900213:mayDecay = off");
  // pythia.readString("4900213:oneChannel = 1 1 100 4900111 4900211");

  // --- [stable] Spin-0 mesons stable. In the original benchmark the diagonal ---
  // --- pion 4900111 decays (to c cbar, the heaviest available state for a ---
  // --- pseudo-scalar) and the off-diagonal 4900211 is stable. Here BOTH are ---
  // --- stable; the decay channels are commented out and not used. ---

  pythia.readString("4900211:mayDecay = off");
  pythia.readString("4900111:mayDecay = off");
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
