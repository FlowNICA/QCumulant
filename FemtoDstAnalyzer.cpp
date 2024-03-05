/**
 * \brief Example of how to read a file (list of files) using StFemtoEvent classes
 *
 * RunFemtoDstAnalyzer.C is an example of reading FemtoDst format.
 * One can use either FemtoDst file or a list of femtoDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 *
 * \author Grigory Nigmatkulov
 * \date May 29, 2018
 */

// This is needed for calling standalone classes
#define _VANILLA_ROOT_

// C++ headers
#include <iostream>

// ROOT headers
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TDirectoryFile.h>

// FemtoDst headers
#include <StFemtoDstReader.h>
#include <StFemtoDst.h>
#include <StFemtoEvent.h>
#include <StFemtoTrack.h>

// Flow method headers
#include <QVector.h>
#include <FlowAnalysisWithEtaSubEventPlane.h>
#include <FlowAnalysisWithThreeEtaSubEventPlane.h>
#include <FlowAnalysisWithFHCalEventPlane.h>
#include <FlowAnalysisWithLeeYangZeros.h>
#include <FlowAnalysisWithScalarProduct.h>
#include <FlowAnalysisWithQCumulant.h>
#include <FlowAnalysisWithHighOrderQCumulant.h>
#include <FlowAnalysisWithQCumulantGenericFramework.h>
#include <FlowAnalysisWithLeeYangZerosEventPlane.h>
#include <FlowAnalysisWithMCEventPlane.h>

// Load libraries (for ROOT_VERSTION_CODE >= 393215)
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
R__LOAD_LIBRARY(libStFemtoDst.so)
#endif

// Flags for flow methods, where *_1 is the first run over the data and *_2 is the second run (need to move these flags to config file!!)
Bool_t ETASUBEVENTPLANE_1 = 0;        // Eta-sub EP (first run)
Bool_t ETASUBEVENTPLANE_2 = 0;        // Eta-sub EP (second run)
Bool_t THREEETASUBEVENTPLANE_1 = 0;   // 3 eta-sub method (first run)
Bool_t THREEETASUBEVENTPLANE_2 = 0;   // 3 eta-sub method (second run)
Bool_t FHCALEVENTPLANE_1 = 0;         // FHCal EP (w.r.t. 1-st order harmonic) (first run)
Bool_t FHCALEVENTPLANE_2 = 0;         // FHCal EP (w.r.t. 1-st order harmonic) (second run)
Bool_t LYZ_SUM_1 = 0;                 // Lee-Yang Zeros using sum generating function (first run)
Bool_t LYZ_SUM_2 = 0;                 // Lee-Yang Zeros using sum generating function (second run)
Bool_t LYZ_PRODUCT_1 = 0;             // Lee-Yang Zeros using product generating function (first run) (integrated with sum GF at the moment, will be separated soon)
Bool_t LYZ_PRODUCT_2 = 0;             // Lee-Yang Zeros using product generating function (second run) (integrated with sum GF at the moment, will be separated soon)
Bool_t SCALARPRODUCT_1 = 0;           // Scalar product using eta-sub method (first run)
Bool_t SCALARPRODUCT_2 = 0;           // Scalar product using eta-sub method (second run)
Bool_t QCUMULANT = 1;                 // Q-Cumulants: 2- and 4-particle cumulants obtained by both standard and subevent methods 
Bool_t HIGHORDERQCUMULANT = 0;        // Q-Cumulants: 2- up to 8-particle cumulants using recursive algorithm
Bool_t QCUMULANTGENERICFRAMWORK = 0;  // Q-Cumulants using genegic framework 
Bool_t LYZEP = 0;                     // one needs to run LYZ_SUM_1 & 2 before set this flag to kTRUE
Bool_t MCEP = 0;                      // MC Event Plane
Int_t harmonic = 2;                   // set harmonic for eta-sub event plane, Q-Cumulants, and scalar product method
Int_t debug = 1;

//Used functions (see them below)
Bool_t isGoodEvent(StFemtoEvent *const &event, Float_t _energy);
Bool_t isGoodTrack(StFemtoTrack *const &track, Float_t _energy, const TVector3 &pVtx);
Bool_t isGoodPID(StFemtoTrack *const &track);
Int_t  GetPID(StFemtoTrack *const &track);
Float_t BBC_GetPhi(const Int_t eastWest, const Int_t tileId);
Double_t GetCent(const Int_t icent);
// Used constants
//Event cuts
const Double_t cutVtxR = 2.;                // Radial vertex cut
// const Double_t cutVtxZvpd = 3.;
// const Int_t cutNTofPoints = 2.;
const Double_t cutVpdVz = 3.;               // Difference between Vz and Vz estimated by VPD - cut apllied only for 200 GeV
const std::map<Float_t, Double_t> cutVtxZEnergy = {{7.7, 70.}, {11.5, 50.}, {19.6, 40}, {27., 40.}, {39., 40.}, {62., 40.}, {200., 30.}};
//Track cuts
// const std::map<Float_t, Double_t> cutDCA = {{7.7, 1.}, {11.5, 1.}, {19.6, 1.}, {27., 1.}, {39., 1.}, {62., 1.}, {200., 3.}};
const std::map<Float_t, Double_t> cutDCA = {{7.7, 1.}, {11.5, 1.}, {19.6, 1.}, {27., 2.}, {39., 1.}, {62., 1.}, {200., 3.}};
// const Double_t cutDCA_PID = 1.;
// const std::vector<Double_t> cutDCA_PIDsys = {0.8,2.};
const Double_t cutEta = 1.;                 // TPC eta cut 
// const std::vector<Double_t> cutEtasys = {0.75};
const Int_t    cutNhits = 15;               // TPC hits
// const std::vector<Double_t> cutNhitssys = {13,18};
const Double_t cutPtotMin = 0.1;            // global track's momentum cut
const Double_t cutNhitsPoss = 0.;           // possible number of hits in TPC
const Double_t cutNhitsRatio = 0.52;        // N_hits / N_hits_poss
const std::map<Float_t, Double_t> cutPtMin = {{7.7, 0.2}, {11.5, 0.2}, {19.6, 0.2}, {27., 0.2}, {39., 0.2}, {62., 0.2}, {200., 0.15}};
const Double_t cutPtMax = 3.6;              // max pt of primary track
const Double_t maxptRF = 3.;                // max pt for reference flow
const Double_t minptRF = 0.2;               // min pt for reference flow
const Double_t cutPMax = 10.;               // max for magnitude of primary track's momentum
const Double_t eta_gap = 0.0375;              // +-0.05, eta-gap between 2 eta sub-event
//PID cuts
const Double_t cutMass2Min = -10.;
const Double_t cutNsigPID = 2.5;



// inFile - is a name of name.FemtoDst.root file or a name
//          of a name.lis(t) files, that contains a list of
//          name1.FemtoDst.root, name2.FemtoDst.root, ... files

//_________________
void FemtoDstAnalyzer(const Char_t *inFile = "st_physics_12150008_raw_4030001.femtoDst.root", const Char_t *outFile = "femtodst_template_output.root", Float_t energy = 27.) 
{
    std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  TStopwatch timer;
  timer.Start();

  #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
    gSystem->Load("libStFemtoDst.so");
  #endif

  // Set up output file
  TFile *fo = new TFile(outFile, "recreate");

  // Set up histograms
  //-----------------------------------------------------------------------------------------------
  //Eventwise parameters
  //hHisto[0] - before event selection
  //hHisto[1] - after event selection
  TH2D *hVtxXY[2];
  TH1D *hVtxZ[2];
  TH1D *hVpdZ[2];

  //Trackwise parameters
  //hHisto[0] - before track selection (but after event selection)
  //hHisto[1] - after track selection
  TH1D *hNhits[2];
  TH1D *hNhitsFit[2];
  TH1D *hNhitsPoss[2];
  TH1D *hNhitsRatio[2];
  TH1D *hPt[2];
  TH1D *hPtot[2];
  TH1D *hPhi[2];
  TH1D *hEta[2];
  TH1D *hDCA[2];

  //PID related parameters
  //hHisto[i] is for particle species: 0-all, 1-pion, 2-kaon, 3-proton
  TH2D *hdEdxQp[4];
  TH2D *hM2Qp[4];
  TH2D *hNsigPiQp[4];
  TH2D *hNsigKaQp[4];
  TH2D *hNsigPrQp[4];
  TH2D *hM2pt[4];
  TH2D *hNsigPipt[4];
  TH2D *hNsigKapt[4];
  TH2D *hNsigPrpt[4];
  //-----------------------------------------------------------------------------------------------

  //Initialization
  //-----------------------------------------------------------------------------------------------
  for (int i = 0; i < 2; i++)
  {
    //Eventwise
    hVtxXY[i] = new TH2D(Form("hVtxXY%i", i), Form("vertex XY %i;Vtx_{X}, [cm];Vtx_{Y}, [cm]", i), 200, -5., 5., 200, -5., 5.);
    hVtxZ[i] = new TH1D(Form("hVtxZ%i", i), Form("vertex Z %i;Vtx_{Z}, [cm];N_{counts}", i), 200, -100., 100.);
    hVpdZ[i] = new TH1D(Form("hVpdZ%i", i), Form("vertex Z %i;Vpd_{Z}, [cm];N_{counts}", i), 200, -100., 100.);

    //Trackwise
    hNhits[i] = new TH1D(Form("hNhits%i", i), Form("N_{hits} %i;N_{hits};N_{counts}", i), 50, 0., 50.);
    hNhitsFit[i] = new TH1D(Form("hNhitsFit%i", i), Form("N_{hist}^{Fit} %i;N_{hist}^{Fit};N_{counts}", i), 50, 0., 50.);
    hNhitsPoss[i] = new TH1D(Form("hNhitsPoss%i", i), Form("N_{hist}^{Poss} %i;N_{hist}^{Poss};N_{counts}", i), 50, 0., 50.);
    hNhitsRatio[i] = new TH1D(Form("hNhitsRatio%i", i), Form("N_{hist}^{Fit}/N_{hist}^{Poss} %i;N_{hist}^{Fit}/N_{hist}^{Poss};N_{counts}", i), 50, 0., 1.);
    hPt[i] = new TH1D(Form("hPt%i", i), Form("p_{T} %i;p_{T}, [GeV/c];N_{counts}", i), 600, 0., 6.);
    hPtot[i] = new TH1D(Form("hPtot%i", i), Form("p_{tot} %i;p, [GeV/c];N_{counts}", i), 600, 0., 12.);
    hPhi[i] = new TH1D(Form("hPhi%i", i), Form("#varphi %i;#varphi, [rad];N_{counts}", i), 600, 0., 6.);
    hEta[i] = new TH1D(Form("hEta%i", i), Form("#eta %i;#eta;N_{counts}", i), 400, -2., 2.);
    hDCA[i] = new TH1D(Form("hDCA%i", i), Form("|DCA| %i;dca, [cm];N_{counts}", i), 200, 0., 5.);
  }
  //PID
  for (int iPID = 0; iPID < 4; iPID++)
  {
    hdEdxQp[iPID] = new TH2D(Form("hdEdxQp%i", iPID), Form("dEdx vs Q*p_{tot} %i;Q*p, [GeV/c];dEdx, [a.u.]", iPID), 1000, -5., 5., 1000, 0., 2e-5);
    hM2Qp[iPID] = new TH2D(Form("hM2Qp%i", iPID), Form("m^{2} vs Q*p_{tot} %i;Q*p, [GeV/c];m^{2}, [GeV/c^{2}]^{2}", iPID), 1000, -5., 5., 1000, -0.8, 1.8);
    hNsigPiQp[iPID] = new TH2D(Form("hNsigPiQp%i", iPID), Form("n#sigma_{#pi} vs Q*p_{tot} %i;Q*p, [GeV/c];n#sigma_{#pi}", iPID), 1000, -5., 5., 600, -3., 3.);
    hNsigKaQp[iPID] = new TH2D(Form("hNsigKaQp%i", iPID), Form("n#sigma_{K} vs Q*p_{tot} %i;Q*p, [GeV/c];n#sigma_{K}", iPID), 1000, -5., 5., 600, -3., 3.);
    hNsigPrQp[iPID] = new TH2D(Form("hNsigPrQp%i", iPID), Form("n#sigma_{p} vs Q*p_{tot} %i;Q*p, [GeV/c];n#sigma_{p}", iPID), 1000, -5., 5., 600, -3., 3.);
    hM2pt[iPID] = new TH2D(Form("hM2pt%i", iPID), Form("m^{2} vs p_{T} %i;p_{T}, [GeV/c];m^{2}, [GeV/c^{2}]^{2}", iPID), 1000, 0., 5., 1000, -0.8, 1.8);
    hNsigPipt[iPID] = new TH2D(Form("hNsigPipt%i", iPID), Form("n#sigma_{#pi} vs p_{T} %i;p_{T}, [GeV/c];n#sigma_{#pi}", iPID), 1000, 0., 5., 600, -3., 3.);
    hNsigKapt[iPID] = new TH2D(Form("hNsigKapt%i", iPID), Form("n#sigma_{K} vs p_{T} %i;p_{T}, [GeV/c];n#sigma_{K}", iPID), 1000, 0., 5., 600, -3., 3.);
    hNsigPrpt[iPID] = new TH2D(Form("hNsigPrpt%i", iPID), Form("n#sigma_{p} vs p_{T} %i;p_{T}, [GeV/c];n#sigma_{p}", iPID), 1000, 0., 5., 600, -3., 3.);
  }
  //-----------------------------------------------------------------------------------------------

  // Counters
  Int_t iPID;

  FlowAnalysisWithEtaSubEventPlane       *flowEtaSub      = nullptr; // Eta-sub Event Plane
  FlowAnalysisWithThreeEtaSubEventPlane  *flowThreeEtaSub = nullptr; // 3-Eta-sub Event Plane
  FlowAnalysisWithFHCalEventPlane        *flowFHCalEP     = nullptr; // FHCal Event Plane w.r.t 1-st harmonic
  FlowAnalysisWithLeeYangZeros           *flowLYZSUM      = nullptr; // Lee-Yang Zeros using sum generating function
  FlowAnalysisWithLeeYangZeros           *flowLYZPROD     = nullptr; // Lee-Yang Zeros using product generating function
  FlowAnalysisWithScalarProduct          *flowSP          = nullptr; // Scalar Product
  FlowAnalysisWithQCumulant              *flowQC          = nullptr; // Q-Cumulant
  FlowAnalysisWithHighOrderQCumulant     *flowHighQC      = nullptr; // 2- to 8-particle correlations using recursive algorithm
  FlowAnalysisWithQCumulantGenericFramework *flowQCGF     = nullptr; // Generic framework for Q-cumulants
  FlowAnalysisWithLeeYangZerosEventPlane *flowLYZEP       = nullptr; // Lee-Yang Zeros Event Plane
  FlowAnalysisWithMCEventPlane           *flowMCEP        = nullptr; // MC Event Plane

  if (ETASUBEVENTPLANE_1) {
    flowEtaSub = new FlowAnalysisWithEtaSubEventPlane();
    flowEtaSub->SetFirstRun(true);
    flowEtaSub->SetHarmonic(harmonic);
    flowEtaSub->SetEtaGap(eta_gap);
    flowEtaSub->Init();
  }
  if (ETASUBEVENTPLANE_2) {
    flowEtaSub = new FlowAnalysisWithEtaSubEventPlane();
    flowEtaSub->SetFirstRun(false);
    flowEtaSub->SetHarmonic(harmonic);
    flowEtaSub->SetEtaGap(eta_gap);
    flowEtaSub->SetDebugFlag(debug);
    flowEtaSub->SetInputFileFromFirstRun("FirstRun.root"); // need to be improve!!!
    flowEtaSub->Init();
  }
  if (THREEETASUBEVENTPLANE_1) {
    flowThreeEtaSub = new FlowAnalysisWithThreeEtaSubEventPlane();
    flowThreeEtaSub->SetFirstRun(true);
    flowThreeEtaSub->SetHarmonic(harmonic);
    flowThreeEtaSub->SetEtaGap(eta_gap);
    flowThreeEtaSub->Init();
  }
  if (THREEETASUBEVENTPLANE_2) {
    flowThreeEtaSub = new FlowAnalysisWithThreeEtaSubEventPlane();
    flowThreeEtaSub->SetFirstRun(false);
    flowThreeEtaSub->SetHarmonic(harmonic);
    flowThreeEtaSub->SetEtaGap(eta_gap);
    flowThreeEtaSub->SetDebugFlag(debug);
    flowThreeEtaSub->SetInputFileFromFirstRun("FirstRun.root");
    flowThreeEtaSub->Init();
  }
  if (FHCALEVENTPLANE_1) {
    flowFHCalEP = new FlowAnalysisWithFHCalEventPlane();
    flowFHCalEP->SetFirstRun(true);
    flowFHCalEP->SetHarmonic(harmonic);
    flowFHCalEP->SetEtaGap(eta_gap);
    flowFHCalEP->Init();
  }
  if (FHCALEVENTPLANE_2) {
    flowFHCalEP = new FlowAnalysisWithFHCalEventPlane();
    flowFHCalEP->SetFirstRun(false);
    flowFHCalEP->SetHarmonic(harmonic);
    flowFHCalEP->SetEtaGap(eta_gap);
    flowFHCalEP->SetDebugFlag(debug);
    flowFHCalEP->SetInputFileFromFirstRun("FirstRun.root");
    flowFHCalEP->Init();
  }

  if (LYZ_SUM_1) {
    flowLYZSUM = new FlowAnalysisWithLeeYangZeros();
    flowLYZSUM->SetDebugFlag(debug);
    flowLYZSUM->SetHarmonic(harmonic);
    flowLYZSUM->SetUseProduct(false);
    flowLYZSUM->SetFirstRun(true);
    // flowLYZSUM->SetUseMultiplicityWeight(false); // true by default
    flowLYZSUM->Init();
  }
  if (LYZ_SUM_2) {
    flowLYZSUM = new FlowAnalysisWithLeeYangZeros();
    flowLYZSUM->SetDebugFlag(debug);
    flowLYZSUM->SetHarmonic(harmonic);
    flowLYZSUM->SetUseProduct(false);
    flowLYZSUM->SetFirstRun(false);
    // flowLYZSUM->SetUseMultiplicityWeight(false); // true by default
    flowLYZSUM->SetInputFileFromFirstRun("FirstRun.root"); // need to be improve!!!
    flowLYZSUM->Init();
  }

  if (LYZ_PRODUCT_1) {
    flowLYZPROD = new FlowAnalysisWithLeeYangZeros();
    flowLYZPROD->SetDebugFlag(debug);
    flowLYZPROD->SetHarmonic(harmonic);
    flowLYZPROD->SetUseProduct(true);
    flowLYZPROD->SetFirstRun(true);
    // flowLYZPROD->SetUseMultiplicityWeight(false); // true by default
    flowLYZPROD->Init();
  }
  if (LYZ_PRODUCT_2) {
    flowLYZPROD = new FlowAnalysisWithLeeYangZeros();
    flowLYZPROD->SetDebugFlag(debug);
    flowLYZPROD->SetHarmonic(harmonic);
    flowLYZPROD->SetUseProduct(true);
    flowLYZPROD->SetFirstRun(false);
    // flowLYZPROD->SetUseMultiplicityWeight(false); // true by default
    flowLYZPROD->SetInputFileFromFirstRun("FirstRun.root"); // need to be improve!!!
    flowLYZPROD->Init();
  }
  if (SCALARPRODUCT_1) {
    flowSP = new FlowAnalysisWithScalarProduct();
    flowSP->SetFirstRun(true);
    flowSP->SetHarmonic(harmonic);
    flowSP->SetEtaGap(eta_gap);
    flowSP->Init();
  }
  if (SCALARPRODUCT_2) {
    flowSP = new FlowAnalysisWithScalarProduct();
    flowSP->SetFirstRun(false);
    flowSP->SetHarmonic(harmonic);
    flowSP->SetEtaGap(eta_gap);
    flowSP->SetInputFileFromFirstRun("FirstRun.root"); // need to be improve!!!
    flowSP->Init();
  }
  if (QCUMULANT) {
    flowQC = new FlowAnalysisWithQCumulant();
    flowQC->SetHarmonic(harmonic);
    flowQC->SetEtaGap(eta_gap);
    flowQC->Init();
  }
  if (HIGHORDERQCUMULANT) {
    flowHighQC = new FlowAnalysisWithHighOrderQCumulant();
    flowHighQC->Init();
  }
  if (QCUMULANTGENERICFRAMWORK)
  {
    flowQCGF = new FlowAnalysisWithQCumulantGenericFramework();
    flowQCGF->SetEtaGap(eta_gap);
    flowQCGF->Init();
  }
  if (LYZEP) {
    flowLYZEP = new FlowAnalysisWithLeeYangZerosEventPlane();
    flowLYZEP->SetInputFileFromFirstAndSecondRun("FirstRun.root", "SecondRun.root");
    flowLYZEP->Init();
  }
  if (MCEP) {
    flowMCEP = new FlowAnalysisWithMCEventPlane();
    flowMCEP->SetDebugFlag(debug);
    flowMCEP->SetHarmonic(harmonic);
    // flowMCEP->SetEtaGap(eta_gap);
    flowMCEP->SetEtaGap(0.);
    flowMCEP->Init();
  }

  Int_t icent, nTracks, fId;
  Double_t cent, pt, eta, phi, charge, ADC_BBC;

  StFemtoDstReader* femtoReader = new StFemtoDstReader(inFile);
  femtoReader->Init();

  // This is a way if you want to spead up IO
  std::cout << "Explicit read status for some branches" << std::endl;
  femtoReader->SetStatus("*",0);
  femtoReader->SetStatus("Event",1);
  femtoReader->SetStatus("Track",1);
  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

  if( !femtoReader->chain() ) {
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInTree = femtoReader->tree()->GetEntries();
  std::cout << "eventsInTree: "  << eventsInTree << std::endl;
  Long64_t events2read = femtoReader->chain()->GetEntries();

  std::cout << "Number of events to read: " << events2read << std::endl;

  // Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

    if (iEvent % 10000 == 0) 
      std::cout << "Working on event #[" << (iEvent+1)
      	      << "/" << events2read << "]" << std::endl;

    Bool_t readEvent = femtoReader->readFemtoEvent(iEvent);
    if( !readEvent ) {
      std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
      break;
    }

    // Retrieve femtoDst
    StFemtoDst *dst = femtoReader->femtoDst();

    // Retrieve event information
    StFemtoEvent *event = dst->event();
    if( !event ) {
      std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
      break;
    }

    // Return primary vertex position
    TVector3 pVtx = event->primaryVertex();

    hVtxXY[0]->Fill(pVtx.X(), pVtx.Y());
    hVtxZ[0]->Fill(pVtx.Z());
    hVpdZ[0]->Fill(event->vpdVz());

    // Cut for Tof Matched
    Int_t nTracksTOF = dst->numberOfTracks();
    Bool_t bCheckTOFMatch = false;
    Int_t number_tof = 0;
    for (Int_t iTrk=0; iTrk<nTracksTOF; iTrk++) {
      StFemtoTrack *femtoTrack = dst->track(iTrk);
      if (!femtoTrack) continue;
      if (femtoTrack->isTofTrack()) {
        number_tof++;
      }
      if (number_tof > 4) {
        bCheckTOFMatch = true;
        break;
      }
    }
    if(!bCheckTOFMatch) continue;

    // Event selection
    if (!isGoodEvent(event, energy)) continue;

    hVtxXY[1]->Fill(pVtx.X(), pVtx.Y());
    hVtxZ[1]->Fill(pVtx.Z());
    hVpdZ[1]->Fill(event->vpdVz());

    icent = 8 - event -> cent9();
    cent = GetCent(icent);

    if (ETASUBEVENTPLANE_1 || ETASUBEVENTPLANE_2)           flowEtaSub->Zero();
    if (THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) flowThreeEtaSub->Zero();
    if (FHCALEVENTPLANE_1 || FHCALEVENTPLANE_2)             flowFHCalEP->Zero();
    if (LYZ_SUM_1 || LYZ_SUM_2)                             flowLYZSUM->Zero();
    if (LYZ_PRODUCT_1 || LYZ_PRODUCT_2)                     flowLYZPROD->Zero();
    if (SCALARPRODUCT_1 || SCALARPRODUCT_2)                 flowSP->Zero();
    if (QCUMULANT)                                          flowQC->Zero();
    if (HIGHORDERQCUMULANT)                                 flowHighQC->Zero();
    if (QCUMULANTGENERICFRAMWORK)                           flowQCGF->Zero();
    if (LYZEP)                                              flowLYZEP->Zero();
    if (MCEP)                                             { flowMCEP->Zero(); flowMCEP->SetPsiRP(0.); }

    // TBI (BBC event-plane)
    if (FHCALEVENTPLANE_1 || FHCALEVENTPLANE_2 || THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2)
    {
      Int_t NmodulesBBC = 24;
      for (Int_t iModule = 0; iModule < NmodulesBBC; iModule++)
      {
        // West BBC
        ADC_BBC = event->bbcAdcWest(iModule);
        phi = BBC_GetPhi(1, iModule);
        eta = -3.; 
        if (FHCALEVENTPLANE_1 || FHCALEVENTPLANE_2) flowFHCalEP->ProcessFirstTrackLoop(eta, phi, ADC_BBC);
        if (THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) flowThreeEtaSub->ProcessFirstTrackLoopFHCal(eta, phi, ADC_BBC);
        // East BBC
        ADC_BBC = event->bbcAdcEast(iModule);
        phi = BBC_GetPhi(0, iModule);
        eta = 3.; 
        if (FHCALEVENTPLANE_1 || FHCALEVENTPLANE_2) flowFHCalEP->ProcessFirstTrackLoop(eta, phi, ADC_BBC);
        if (THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) flowThreeEtaSub->ProcessFirstTrackLoopFHCal(eta, phi, ADC_BBC);
      }
    }

    // Track analysis
    nTracks = dst->numberOfTracks();

    // Track loop
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

      // Retrieve i-th femto track
      StFemtoTrack *femtoTrack = dst->track(iTrk);

      if (!femtoTrack) continue;

      hNhits[0]->Fill(femtoTrack->nHits());
      hNhitsFit[0]->Fill(femtoTrack->nHitsFit());
      hNhitsPoss[0]->Fill(femtoTrack->nHitsPoss());
      if (femtoTrack->nHitsPoss() != 0)
      {
        hNhitsRatio[0]->Fill((Double_t)(femtoTrack->nHitsFit() / femtoTrack->nHitsPoss()));
      }
      hPt[0]->Fill(femtoTrack->pt());
      hPtot[0]->Fill(femtoTrack->ptot());
      hEta[0]->Fill(femtoTrack->eta());
      hPhi[0]->Fill(femtoTrack->phi());
      hDCA[0]->Fill(femtoTrack->gDCA(pVtx).Mag());

      //Track selection
      if ( !isGoodTrack(femtoTrack, energy, pVtx) ) continue;

      hNhits[1]->Fill(femtoTrack->nHits());
      hNhitsFit[1]->Fill(femtoTrack->nHitsFit());
      hNhitsPoss[1]->Fill(femtoTrack->nHitsPoss());
      if (femtoTrack->nHitsPoss() != 0)
      {
        hNhitsRatio[1]->Fill((Double_t)(femtoTrack->nHitsFit() / femtoTrack->nHitsPoss()));
      }
      hPt[1]->Fill(femtoTrack->pt());
      hPtot[1]->Fill(femtoTrack->ptot());
      hEta[1]->Fill(femtoTrack->eta());
      hPhi[1]->Fill(femtoTrack->phi());
      hDCA[1]->Fill(femtoTrack->gDCA(pVtx).Mag());

      hdEdxQp[0]->Fill(femtoTrack->charge() * femtoTrack->ptot(), femtoTrack->dEdx());
      hNsigPiQp[0]->Fill(femtoTrack->charge() * femtoTrack->ptot(), femtoTrack->nSigmaPion());
      hNsigKaQp[0]->Fill(femtoTrack->charge() * femtoTrack->ptot(), femtoTrack->nSigmaKaon());
      hNsigPrQp[0]->Fill(femtoTrack->charge() * femtoTrack->ptot(), femtoTrack->nSigmaProton());
      hNsigPipt[0]->Fill(femtoTrack->pt(), femtoTrack->nSigmaPion());
      hNsigKapt[0]->Fill(femtoTrack->pt(), femtoTrack->nSigmaKaon());
      hNsigPrpt[0]->Fill(femtoTrack->pt(), femtoTrack->nSigmaProton());

      // Check if track has TOF signal - for TOF related distributions
      if ( femtoTrack->isTofTrack() ) {
        hM2Qp[0]->Fill(femtoTrack->charge() * femtoTrack->ptot(), femtoTrack->massSqr());
        hM2pt[0]->Fill(femtoTrack->pt(), femtoTrack->massSqr());
      }

      iPID = GetPID(femtoTrack);
      if (iPID == -1) continue;

      hdEdxQp[iPID+1]->Fill(femtoTrack->charge() * femtoTrack->ptot(), femtoTrack->dEdx());
      hNsigPiQp[iPID+1]->Fill(femtoTrack->charge() * femtoTrack->ptot(), femtoTrack->nSigmaPion());
      hNsigKaQp[iPID+1]->Fill(femtoTrack->charge() * femtoTrack->ptot(), femtoTrack->nSigmaKaon());
      hNsigPrQp[iPID+1]->Fill(femtoTrack->charge() * femtoTrack->ptot(), femtoTrack->nSigmaProton());
      hNsigPipt[iPID+1]->Fill(femtoTrack->pt(), femtoTrack->nSigmaPion());
      hNsigKapt[iPID+1]->Fill(femtoTrack->pt(), femtoTrack->nSigmaKaon());
      hNsigPrpt[iPID+1]->Fill(femtoTrack->pt(), femtoTrack->nSigmaProton());

      // Check if track has TOF signal - for TOF related distributions
      if ( femtoTrack->isTofTrack() ) {
        hM2Qp[iPID+1]->Fill(femtoTrack->charge() * femtoTrack->ptot(), femtoTrack->massSqr());
        hM2pt[iPID+1]->Fill(femtoTrack->pt(), femtoTrack->massSqr());
      }

    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    if ( (LYZ_PRODUCT_1 || LYZ_PRODUCT_2) && flowLYZPROD->GetUseMultiplicityWeight() )
    { // Zero Track loop for Product LYZ
      for (Int_t iTrk = 0; iTrk < nTracks; iTrk++)
      {
        StFemtoTrack *femtoTrack = dst->track(iTrk);
        if ( !isGoodTrack(femtoTrack, energy, pVtx) ) continue;
        pt = femtoTrack->pt();
        phi = femtoTrack->phi();
        if (pt > minptRF && pt < maxptRF)
        {
          flowLYZPROD->ProcessZeroTrackLoopRP();
        }
      }
    }
    for (Int_t iTrk = 0; iTrk < nTracks; iTrk++)
    { // First Track loop

      StFemtoTrack *femtoTrack = dst->track(iTrk);
      if ( !isGoodTrack(femtoTrack, energy, pVtx) ) continue;
      pt = femtoTrack->pt();
      phi = femtoTrack->phi();
      eta = femtoTrack->phi();
      charge = femtoTrack->charge();
      fId = (charge>0) ? (GetPID(femtoTrack) + 1) : (GetPID(femtoTrack) + 5);
      
      if (pt > minptRF && pt < maxptRF)
      { // Reference Flow pt cut
        if (ETASUBEVENTPLANE_1 || ETASUBEVENTPLANE_2) flowEtaSub->ProcessFirstTrackLoop(eta, phi, pt);
        if (THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) flowThreeEtaSub->ProcessFirstTrackLoopTPC(eta, phi, pt);
        if (SCALARPRODUCT_1 || SCALARPRODUCT_2) flowSP->ProcessFirstTrackLoop(eta, phi, pt);
        if (LYZ_SUM_1 || LYZ_SUM_2) flowLYZSUM->ProcessFirstTrackLoopRP(phi, pt, icent);
        if (LYZ_PRODUCT_1 || LYZ_PRODUCT_2) flowLYZPROD->ProcessFirstTrackLoopRP(phi, pt, icent);
        if (QCUMULANT) flowQC->ProcessFirstTrackLoopRP(eta, phi);
        if (HIGHORDERQCUMULANT) flowHighQC->ProcessFirstTrackLoopRP(phi);
        if (QCUMULANTGENERICFRAMWORK) flowQCGF->ProcessFirstTrackLoopRP(eta, phi);
        if (MCEP) flowMCEP->ProcessFirstTrackLoop(eta, phi, 1.);
      }
      if (QCUMULANT)                        flowQC      -> ProcessFirstTrackLoopPOI(eta, phi, pt, fId, charge, minptRF, maxptRF);
      if (QCUMULANTGENERICFRAMWORK)         flowQCGF    -> ProcessFirstTrackLoopPOI(eta, phi, pt, fId, charge, minptRF, maxptRF);
      if (LYZ_SUM_1 || LYZ_SUM_2)           flowLYZSUM  -> ProcessFirstTrackLoopPOI(pt);
      if (LYZ_PRODUCT_1 || LYZ_PRODUCT_2)   flowLYZPROD -> ProcessFirstTrackLoopPOI(pt);
    } // end of First Track loop
    
    if (ETASUBEVENTPLANE_1 || ETASUBEVENTPLANE_2)           flowEtaSub->ProcessEventAfterFirstTrackLoop(cent);
    if (THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) flowThreeEtaSub->ProcessEventAfterFirstTrackLoop(cent);    
    if (FHCALEVENTPLANE_1 || FHCALEVENTPLANE_2)             flowFHCalEP->ProcessEventAfterFirstTrackLoop(cent);
    if (SCALARPRODUCT_1 || SCALARPRODUCT_2)                 flowSP->ProcessEventAfterFirstTrackLoop(cent);
    if (LYZ_SUM_1 || LYZ_SUM_2)                             flowLYZSUM->ProcessEventAfterFirstTrackLoop(icent);
    if (LYZ_PRODUCT_1 || LYZ_PRODUCT_2)                     flowLYZPROD->ProcessEventAfterFirstTrackLoop(icent);
    if (QCUMULANT)                                          flowQC->ProcessEventAfterFirstTrackLoop(icent);
    if (HIGHORDERQCUMULANT)                                 flowHighQC->ProcessEventAfterFirstTrackLoop(icent);
    if (QCUMULANTGENERICFRAMWORK)                           flowQCGF->ProcessEventAfterFirstTrackLoop(icent);
    if (LYZEP)                                              flowLYZEP->ProcessEventAfterFirstTrackLoop(icent);
    if (MCEP)                                               flowMCEP->ProcessEventAfterFirstTrackLoop(cent);
    
    if (ETASUBEVENTPLANE_2 || FHCALEVENTPLANE_2 || THREEETASUBEVENTPLANE_2 || LYZ_SUM_2 || LYZ_PRODUCT_2 || SCALARPRODUCT_2 || LYZEP || MCEP)
    {
      for (Int_t iTrk = 0; iTrk < nTracks; iTrk++)
      { // 2nd Track loop

        StFemtoTrack *femtoTrack = dst->track(iTrk);
        if ( !isGoodTrack(femtoTrack, energy, pVtx) ) continue;
        pt = femtoTrack->pt();
        phi = femtoTrack->phi();
        eta = femtoTrack->phi();
        charge = femtoTrack->charge();
        fId = (charge>0) ? (GetPID(femtoTrack) + 1) : (GetPID(femtoTrack) + 5);
        
        if (ETASUBEVENTPLANE_2)       flowEtaSub->ProcessSecondTrackLoop(eta, phi, pt, cent, fId, charge);
        if (THREEETASUBEVENTPLANE_2)  flowThreeEtaSub->ProcessSecondTrackLoop(eta, phi, pt, cent, fId, charge);
        if (FHCALEVENTPLANE_2)        flowFHCalEP->ProcessSecondTrackLoop(eta, phi, pt, cent, fId, charge);
        if (SCALARPRODUCT_2)          flowSP->ProcessSecondTrackLoop(eta, phi, pt, cent, fId, charge);
        if (LYZ_SUM_2)                flowLYZSUM->ProcessSecondTrackLoop(phi, pt, icent);
        if (LYZ_PRODUCT_2)            flowLYZPROD->ProcessSecondTrackLoop(phi, pt, icent);
        if (LYZEP)                    flowLYZEP->ProcessSecondTrackLoop(eta, phi, pt, cent, fId, charge);
        if (MCEP)                     flowMCEP->ProcessSecondTrackLoop(eta, phi, pt, cent, fId, charge);
      } // end of 2nd Track loop
    }
    
  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  femtoReader->Finish();

  // Writing output
  const Int_t nMethods = 11;
  TString dirNameMethod[nMethods] = {"ETASUBEP","ETA3SUBEP","FHCALEP","SP","LYZSUM","LYZPROD","QC","HQC","QCGF","LYZEP","MCEP"};
  fo->cd();
  TDirectoryFile *dirFileFinal[nMethods] = {nullptr};
  for(Int_t i=0;i<nMethods;i++)
  {
    dirFileFinal[i] = new TDirectoryFile(dirNameMethod[i].Data(),dirNameMethod[i].Data());
  }
  // if (ETASUBEVENTPLANE_1 || ETASUBEVENTPLANE_2) flowEtaSub->SaveHist();
  // if (THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) flowThreeEtaSub->SaveHist();
  // if (FHCALEVENTPLANE_1 || FHCALEVENTPLANE_2) flowFHCalEP->SaveHist();
  // if (SCALARPRODUCT_1 || SCALARPRODUCT_2) flowSP->SaveHist();
  // if (LYZ_SUM_1 || LYZ_SUM_2) flowLYZSUM->SaveHist();
  // if (LYZ_PRODUCT_1 || LYZ_PRODUCT_2) flowLYZPROD->SaveHist();
  // if (QCUMULANT) flowQC->SaveHist();
  // if (HIGHORDERQCUMULANT) flowHighQC->SaveHist();
  // if (QCUMULANTGENERICFRAMWORK) flowQCGF->SaveHist();
  // if (LYZEP) flowLYZEP->SaveHist();

  if (ETASUBEVENTPLANE_1 || ETASUBEVENTPLANE_2)           flowEtaSub->SaveHist(dirFileFinal[0]);
  if (THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) flowThreeEtaSub->SaveHist(dirFileFinal[1]);
  if (FHCALEVENTPLANE_1 || FHCALEVENTPLANE_2)             flowFHCalEP->SaveHist(dirFileFinal[2]);
  if (SCALARPRODUCT_1 || SCALARPRODUCT_2)                 flowSP->SaveHist(dirFileFinal[3]);
  if (LYZ_SUM_1 || LYZ_SUM_2)                             flowLYZSUM->SaveHist(dirFileFinal[4]);
  if (LYZ_PRODUCT_1 || LYZ_PRODUCT_2)                     flowLYZPROD->SaveHist(dirFileFinal[5]);
  if (QCUMULANT)                                          flowQC->SaveHist(dirFileFinal[6]);
  if (HIGHORDERQCUMULANT)                                 flowHighQC->SaveHist(dirFileFinal[7]);
  if (QCUMULANTGENERICFRAMWORK)                           flowQCGF->SaveHist(dirFileFinal[8]);
  if (LYZEP)                                              flowLYZEP->SaveHist(dirFileFinal[9]);
  if (MCEP)                                               flowMCEP->SaveHist(dirFileFinal[10]);

  for (int i = 0; i < 2; i++)
  {
    //Eventwise
    hVtxXY[i] -> Write();
    hVtxZ[i] -> Write();
    hVpdZ[i] -> Write();

    //Trackwise
    hNhits[i] -> Write();
    hNhitsFit[i] -> Write();
    hNhitsPoss[i] -> Write();
    hNhitsRatio[i] -> Write();
    hPt[i] -> Write();
    hPtot[i] -> Write();
    hPhi[i] -> Write();
    hEta[i] -> Write();
    hDCA[i] -> Write();
  }
  //PID
  for (int iPID = 0; iPID < 4; iPID++)
  {
    hdEdxQp[iPID] -> Write();
    hM2Qp[iPID] -> Write();
    hNsigPiQp[iPID] -> Write();
    hNsigKaQp[iPID] -> Write();
    hNsigPrQp[iPID] -> Write();
    hM2pt[iPID] -> Write();
    hNsigPipt[iPID] -> Write();
    hNsigKapt[iPID] -> Write();
    hNsigPrpt[iPID] -> Write();
  }
  fo->Close();


  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	          << std::endl;

  timer.Stop();
  timer.Print();
}

Bool_t isGoodEvent(StFemtoEvent *const &event, Float_t _energy)
{
  if (!event)                                   return false;
  if (event == nullptr)                         return false;
  // Reject vertices that are far from the central membrane along the beam
  if (event->primaryVertex().Perp() > cutVtxR)                            return false;
  if (TMath::Abs(event->primaryVertex().Z()) > cutVtxZEnergy.at(_energy)) return false;

  if ((_energy == 200.) && TMath::Abs(event->primaryVertex().Z() - event->vpdVz()) > cutVpdVz)  return false;

  return true;
}

Bool_t isGoodTrack(StFemtoTrack *const &track, Float_t _energy, const TVector3 &pVtx)
{
  if (!track)                                                           return false;

  // Primary track cut  
  if (!track->isPrimary())                                              return false;
  if ( ( track -> dEdx() ) == 0. )                                      return false; // is this one of primary track cuts?

  // Quality track cut
  if (track->nHitsFit() <= cutNhits)                                     return false; // what is the different between nHits and nHitsFit? 
  if ((Double_t)track->nHitsFit() / track->nHitsPoss() <= cutNhitsRatio) return false;
  // if (track->nHits() < cutNhits)                                        return false;
  // if ((Double_t)track->nHits() / track->nHitsPoss() < cutNhitsRatio)    return false;
  if (track->nHitsPoss() <= cutNhitsPoss)                               return false;
  if (track->gDCA(pVtx).Mag() >= cutDCA.at(_energy))                    return false;

  // Kinematic track cut
  if (TMath::Abs(track->eta()) >= cutEta)                               return false;
  if (track->pt() <= cutPtMin.at(_energy))                              return false;
  if (track->pt() > cutPtMax)                                           return false;
  if (track->ptot() > cutPMax)                                          return false;
  if (track->gMom().Mag() < cutPtotMin)                                 return false;

  return true;
}

Bool_t isGoodPID(StFemtoTrack *const &track)
{
  // Check if track has TOF signal
  if (!track->isTofTrack()) return false;
  if (track->massSqr() < cutMass2Min) return false;
  //NhitsDedx cut applied in StFemtoDst?
  // ToFYLocal cut applied in StFemtoDst?
  return true;
}

Int_t  GetPID(StFemtoTrack *const &track)
{
  Int_t i_part = -1;

  //TPC-only
  if (!track->isTofTrack())
  {
    //pion id
    if (track->ptot() >= 0.2 && track->ptot() < 0.6 &&
        TMath::Abs(track->nSigmaPion()) < 2)
    {
      i_part = 0;
    }
    // kaon id
    if (track->ptot() >= 0.2 && track->ptot() < 0.5 &&
        TMath::Abs(track->nSigmaKaon()) < 2)
    {
      i_part = 1;
    }
    // proton id
    if (track->ptot() >= 0.4 && track->ptot() < 0.9 &&
        TMath::Abs(track->nSigmaProton()) < 2)
    {
      i_part = 2;
    }
    //pion id
    if (track->ptot() >= 0.6 && track->ptot() < 0.7 &&
        TMath::Abs(track->nSigmaPion()) < 2 &&
        TMath::Abs(track->nSigmaKaon()) > 2)
    {
      i_part = 0;
    }
    // kaon id
    if (track->ptot() >= 0.5 && track->ptot() < 0.7 &&
        TMath::Abs(track->nSigmaKaon()) < 2 &&
        TMath::Abs(track->nSigmaPion()) > 3)
    {
      i_part = 1;
    }
    // proton id
    if (track->ptot() >= 0.9 && track->ptot() < 1.2 &&
        TMath::Abs(track->nSigmaProton()) < 2 &&
        TMath::Abs(track->nSigmaPion()) > 3)
    {
      i_part = 2;
    }
  }

  //TPC+TOF
  if (isGoodPID(track))
  {
    // pion id
    if (track->ptot() >= 0.2 && track->ptot() < 3.4 &&
        TMath::Abs(track->nSigmaPion()) < 3 &&
        track->massSqr() >= -0.15 && track->massSqr() < 0.1)
    {
      i_part = 0;
    }
    // kaon id
    if (track->pt() >= 0.2 && track->ptot() < 3.4 &&
        TMath::Abs(track->nSigmaKaon()) < 3 &&
        track->massSqr() >= 0.2 && track->massSqr() < 0.32)
    {
      i_part = 1;
    }
    // proton id
    if (track->ptot() >= 0.4 && track->ptot() < 3.4 &&
        TMath::Abs(track->nSigmaProton()) < 3 &&
        track->massSqr() >= 0.75 && track->massSqr() < 1.2)
    {
      i_part = 2;
    }
  }

  return i_part;
}

//--------------------------------------------------------------------------------------------------------------------------//
//this function simply connects the gain values read in to the BBC azimuthal distribution
//since tiles 7 and 9 (+ 13 and 15) share a gain value it is ambiguous how to assign the geometry here
//I prefer assigning the angle between the tiles thus "greying out" the adcs. 
//Others have assigned all of the adc to one (exclusive) or the the other. 
Float_t BBC_GetPhi(const Int_t eastWest, const Int_t tileId)
{
  //float GetPhiInBBC(int eastWest, int bbcN) { //tileId=0 to 23
  const float Pi = TMath::Pi();
  const float Phi_div = Pi / 6;
  float bbc_phi = Phi_div;
  switch(tileId) {
    case 0: bbc_phi = 3.*Phi_div;
  break;
    case 1: bbc_phi = Phi_div;
  break;
    case 2: bbc_phi = -1.*Phi_div;
  break;
    case 3: bbc_phi = -3.*Phi_div;
  break;
    case 4: bbc_phi = -5.*Phi_div;
  break;
    case 5: bbc_phi = 5.*Phi_div;
  break;
    //case 6: bbc_phi= (mRndm.Rndm() > 0.5) ? 2.*Phi_div:4.*Phi_div;	//tiles 7 and 9 are gained together we randomly assign the gain to one XOR the other
    case 6: bbc_phi = 3.*Phi_div;
  break;
    case 7: bbc_phi = 3.*Phi_div;
  break;
    case 8: bbc_phi = Phi_div;
  break;
    case 9: bbc_phi = 0.;
  break;
    case 10: bbc_phi = -1.*Phi_div;
  break;
    //case 11: bbc_phi = (mRndm.Rndm() > 0.5) ? -2.*Phi_div:-4.*Phi_div;	//tiles 13 and 15 are gained together
    case 11: bbc_phi = -3.*Phi_div;
  break;
    case 12: bbc_phi = -3.*Phi_div;
  break;
    case 13: bbc_phi = -5.*Phi_div;
  break;
    case 14: bbc_phi = Pi;
  break;
    case 15: bbc_phi = 5.*Phi_div;
  break;
  }

  //if we're looking at the east BBC we need to flip around x in the STAR coordinates, 
  //a line parallel to the beam would go through tile 1 on the W BBC and tile 3 on the 
  if(0 == eastWest){
    if (bbc_phi > -0.001){ //this is not a >= since we are talking about finite adcs -- not to important
      bbc_phi = Pi - bbc_phi;
    }
    else {
      bbc_phi= -Pi - bbc_phi;
    }
  }

  if(bbc_phi < 0.0) bbc_phi += 2.*Pi;
  if(bbc_phi > 2.*Pi) bbc_phi -= 2.*Pi;

  return bbc_phi;
}

Double_t GetCent(Int_t icent)
{
  if (icent ==  0) return 2.5;
  if (icent ==  1) return 7.5;
  if (icent ==  2) return 15.;
  if (icent ==  3) return 25.;
  if (icent ==  4) return 35.;
  if (icent ==  5) return 45.;
  if (icent ==  6) return 55.;
  if (icent ==  7) return 65.;
  if (icent ==  8) return 75.;
  return -1;
}

int main(int argc, char **argv)
{
  TString iFileName, oFileName;
  Float_t centerOfMassEnergy = 0.;

  if (argc < 5)
  {
    std::cerr << "./FemtoDstAnalyzer -i INPUT -o OUTPUT -cme COLLISION_ENERGY" << std::endl;
    return 1;
  }
  for (Int_t i = 1; i < argc; i++)
  {
    if (std::string(argv[i]) != "-i" &&
        std::string(argv[i]) != "-o" &&
        std::string(argv[i]) != "-cme")
    {
      std::cerr << "\n[ERROR]: Unknown parameter " << i << ": " << argv[i] << std::endl;
      return 2;
    }
    else
    {
      if (std::string(argv[i]) == "-i" && i != argc - 1)
      {
        iFileName = argv[++i];
        continue;
      }
      if (std::string(argv[i]) == "-i" && i == argc - 1)
      {
        std::cerr << "\n[ERROR]: Input file name was not specified " << std::endl;
        return 3;
      }
      if (std::string(argv[i]) == "-o" && i != argc - 1)
      {
        oFileName = argv[++i];
        continue;
      }
      if (std::string(argv[i]) == "-o" && i == argc - 1)
      {
        std::cerr << "\n[ERROR]: Output file name was not specified " << std::endl;
        return 4;
      }
      if (std::string(argv[i]) == "-cme" && i != argc - 1)
      {
        centerOfMassEnergy = atof(argv[++i]);
        continue;
      }
      if (std::string(argv[i]) == "-cme" && i == argc - 1)
      {
        std::cerr << "\n[ERROR]: Collision energy was not specified " << std::endl;
        return 1;
      }
    }
  }
  FemtoDstAnalyzer(iFileName, oFileName, centerOfMassEnergy);

  return 0;
}
