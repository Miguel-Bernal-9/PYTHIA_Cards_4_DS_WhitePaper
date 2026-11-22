// --- [stop32] In 8.307 the Hidden Valley string stop mass is hard-coded as ---
// --- 1.5 * m(pi_v) inside a private class (HVStringZ) and is not a setting. ---
// --- To change it from here, the private members of Pythia are made ---
// --- accessible ONLY while Pythia.h is read. The standard library is ---
// --- included first so it is not affected. Nothing is recompiled: the ---
// --- official 8.307 library is used as it is. ---
// --- COMPILE ONLY AGAINST PYTHIA 8.307: the member names used below ---
// --- (hiddenvalleyFrag, hvStringFrag, hvZSel, mhvMeson) change between ---
// --- versions. ---

#include <bits/stdc++.h>
#define private public
#define protected public
#include "Pythia8/Pythia.h"
#undef private
#undef protected

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
//  Output name       : Nc3Nf8_sFoff_pp_3pi_decay_lam_10_probvec_75_py8307_stop32
//                      (.hepmc and _HVparameters.txt in outDir)
//  Pythia version    : 8.307 ONLY  (the [stop32] trick below uses 8.307
//                      internal class members: with another version it does
//                      not compile or does not do what it says)
//  Process           : pp, 13 TeV -> Z' (m = 1 TeV, width 20 GeV, 800-1200 GeV)
//                      -> qv qvbar. The Z' is forced to decay to dark quarks only
//  Dark sector       : SU(3) (Nc = 3), Nf = 8, separateFlav = off
//                      Lambda = 10 GeV, alphaOrder = 1, probVector = 0.75
//                      m(qv) = 10.2 GeV, m(pi_v) = 6 GeV, m(rho_v) = 25 GeV
//                      (the same mass for every flavour)
//  Dark shower cutoff: pTminFSR is raised internally to 1.1*Lambda = 11 GeV
//                      (the warning 'Hidden Valley pTmin too low' is expected)
//  HV stop mass      : 32.24 GeV, FORCED (see [stop32]). This is the 8.313+ value
//                      1.0 * m(rho_v)/m(rho0) = 25/0.77549, NOT the 8.307 default
//                      (9 GeV). Test sample to reproduce 8.313/8.315 in 8.307
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
//  Validation        : 1000 events, seed 12345, Pythia 8.307 (final state, per event):
//                      <N(pi_v)> = 3.21, <N(rho_v)> = 8.98,
//                      <pT(pi_v)> = 33.6 GeV, <pT(rho_v)> = 67.2 GeV
//                      (statistical uncertainty about 1-2% on each number)
//                      Nf = 8, sFon and sFoff with the same seed give IDENTICAL
//                      events: only the dark meson labels differ
// ---------------------------------------------------------------------
// Tags used in the comments of the code:
//
// [SFoff] this card was made from the Nf = 8 separateFlav = on card, and only
// separateFlav was changed to off. It therefore keeps the 8 Z' -> qv_i qv_i
// channels (0.05 each) and the loops that set the masses and mayDecay of all
// 36 + 36 meson species. With separateFlav = off only 4900111/211/113/213 are
// written to the output, so the extra settings are harmless; this was checked:
// with the same seed this card and the sFon card give the same kinematics.
// (The Nf = 3 sFoff card uses instead a single Z' -> 4900101 channel with
// BR 0.4, as the Snowmass .cmnd. Both choices are equivalent here because all
// dark quarks have the same mass.)
//
// [vecfix] all vector mesons 4900ij3 (8 >= i >= j >= 1, 36 species) set to
// 25 GeV. This is the fix of the 8.313 bug (it picks the lightest vector
// among ALL 4900ij3, even unused ones with the 10 GeV default). It has NO
// effect on the events in 8.307; it is kept so that the 8.307, 8.313 and
// 8.315 cards use identical particle data.
//
// [stop32] TEST: the HV string stop mass is set after pythia.init() to the
// 8.313+ prescription StringFragmentation:stopMass * m(rho_v)/m(rho0)
// = 1.0 * 25/0.77549 = 32.24 GeV, instead of the 8.307 value 1.5 * m(pi_v)
// = 9 GeV. Everything else is plain 8.307. The value really used is printed
// in the log and in _HVparameters.txt ('HV stop mass used in the generation').
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
const string outName = "Nc3Nf8_sFoff_pp_3pi_decay_lam_10_probvec_75_py8307_stop32";

// --- [SFoff] Number of dark flavours and dark sector masses (used in the loops below) ---

const int    NFLAV = 8;
const double MASS_QV  = 10.2;   // dark quarks
const double MASS_PIV  = 6.0;    // dark pions (spin 0)
const double MASS_RHOV = 25.0;   // dark rhos  (spin 1)

// --- PDG codes of the dark mesons 4900ijs, i >= j (s = 1 pion, 3 rho) ---

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

  // --- check that all masses are the ones requested ---

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

  // --- [stop32] HV string stop mass actually used, read from Pythia itself ---
  // --- (8.307 default: 1.5 * m(4900111) = 9 GeV) ---
  os << "  HV stop mass used in the generation = "
     << pythia.hadronLevel.hiddenvalleyFrag.hvStringFrag.stopMass << " GeV"
     << "   (8.307 default 1.5 * m(4900111) = "
     << 1.5 * pythia.particleData.m0(4900111) << " GeV)\n";

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

  // [SFoff] Kept from the sFon card: the dark BR (0.4) split equally among the
  // 8 dark quarks (0.4 / 8 = 0.05 each). With separateFlav = off this gives
  // the same events as a single 4900101 channel with BR 0.4 (all qv have the
  // same mass); see [SFoff] in the header.

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
                    "4900105 4900106 4900107 4900108");                   // [SFoff]

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
  // --- [SFoff] one mass per flavour, 4900101 ... 4900108 (kept from the sFon card) ---

  for (int i = 1; i <= NFLAV; ++i)
    pythia.readString(to_string(4900100 + i) + ":m0 = " + to_string(MASS_QV));

  // --- [vecfix] same mass for ALL vector mesons 4900ij3 (8 >= i >= j >= 1), ---
  // --- as in the 8.313/8.315 cards ---
  for (int i = 1; i <= 8; ++i)
    for (int j = 1; j <= i; ++j)
      pythia.readString(to_string(4900003 + 100 * i + 10 * j) + ":m0 = 25");

  // --- 8) Set gauge group and number of flavours ---
  // --- caution: user should not use Nc=2 or Nf=1 without consulting theory experts ---

  pythia.readString("HiddenValley:Ngauge  = 3");
  pythia.readString("HiddenValley:Nflav  = 8");
  pythia.readString("HiddenValley:separateFlav  = off");                  // [SFoff] off

  // --- 9) Set other darksector parameters, must set HiddenValley:alphaOrder = 1.0 ---

  pythia.readString("HiddenValley:Lambda = 10");
  pythia.readString("HiddenValley:FSR = on");
  pythia.readString("HiddenValley:fragment = on");
  pythia.readString("HiddenValley:alphaOrder = 1");
  pythia.readString("HiddenValley:probVector = 0.75");

  // --- 8.307: Lambda is used directly (setLambda does not exist yet) and the ---
  // --- hadronization uses the 8.307 defaults aLund = 0.3, bmqv2 = 0.8, ---
  // --- sigmamqv = 0.5, rFactqv = 1 (setabsigma, bLund, sigmaLund, rFact are ---
  // --- 8.313+ names and must NOT be used here: init() would fail) ---

  // --- 10) Set masses of spin-0 and spin-1 mesons ---
  // --- [SFoff] every species 4900ij1 (pion) and 4900ij3 (rho), 8 >= i >= j >= 1: ---
  // --- 36 of each, kept from the sFon card. With separateFlav = off only ---
  // --- 4900111/211/113/213 are used; the others are harmless. ---
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
  // --- [SFoff] all 72 species, kept from the sFon card (the decay channels stay commented, as in the ---
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

  // --- [stop32] Set the HV string stop mass to the 8.313+ prescription ---
  // --- stopMass(HV) = StringFragmentation:stopMass * m(rho_v,min)/m(rho0). ---
  // --- In 8.307 HVStringZ::stopMass() returns 1.5 * mhvMeson (used only for ---
  // --- the stop mass), and the HV StringFragmentation copies that value into ---
  // --- its own member 'stopMass' during init(). Both are set here, so the ---
  // --- value actually used in the HV string fragmentation is the new one. ---
  // --- Must be done AFTER pythia.init() (init() sets both). ---
  {
    double mVec = pythia.particleData.m0(4900113);
    if (pythia.settings.flag("HiddenValley:separateFlav")) {
      int nFl = pythia.settings.mode("HiddenValley:nFlav");
      for (int i = 2; i <= nFl; ++i)
        for (int j = 1; j <= i; ++j)
          mVec = min(mVec, pythia.particleData.m0(4900003 + 100 * i + 10 * j));
    }
    double stopTarget = pythia.settings.parm("StringFragmentation:stopMass")
                      * mVec / pythia.particleData.m0(113);
    HiddenValleyFragmentation& hvFrag = pythia.hadronLevel.hiddenvalleyFrag;
    double stopOld = hvFrag.hvStringFrag.stopMass;
    hvFrag.hvZSel.mhvMeson      = stopTarget / 1.5;   // HVStringZ::stopMass()
    hvFrag.hvStringFrag.stopMass = stopTarget;        // value used in fragmentation
    cout << "\n [stop32] HV string stop mass changed from " << stopOld
         << " GeV to " << hvFrag.hvStringFrag.stopMass << " GeV\n";
  }

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
