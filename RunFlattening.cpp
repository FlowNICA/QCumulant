// C++ headers
#include <iostream>
#include <fstream>
// ROOT headers
#include <TStopwatch.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TF2.h>
#include <TDatabasePDG.h>
#include <TProfile.h>
#include <TString.h>
#include <TEnv.h>

// PicoDst headers
#include <PicoDstMCEvent.h>
#include <PicoDstRecoEvent.h>
#include <PicoDstMCTrack.h>
#include <PicoDstRecoTrack.h>
#include <PicoDstFHCal.h>

#include <IReader.h>
#include <PicoDstReader.h>
#include <McPicoReader.h>
#include <ToyModelReader.h>
// Flow method headers
#include <QVector.h>
#include "utilities.C"

using std::cout;
using std::cerr;
using std::endl;

Bool_t readMCTracks = 0; // 0 - read reco tracks, 1 - read MC tracks
Int_t harmonic = 2; // set harmonic for eta-sub event plane, Q-Cumulants, and scalar product method
Bool_t bMotherIDcut = 1;
// Kinetic cuts by default if not using config file
Double_t maxpt = 3.0;     // max pt for differential flow
Double_t minpt = 0.;      // min pt for differential flow
Double_t maxptRF = 3.;    // max pt for reference flow
Double_t minptRF = 0.2;   // min pt for reference flow
Double_t eta_cut = 1.5;   // pseudorapidity acceptance window for flow measurements 
Double_t eta_gap = 0.05;  // +-0.05, eta-gap between 2 eta sub-event
Int_t Nhits_cut = 16;     // minimum nhits of reconstructed tracks
Double_t DCAcut = 2.0;    // 3*sigma of DCA distribution 
Double_t pid_probability = 0.9;
Long64_t Nevents = -1;

Int_t debug = 1;

std::string format = "";

void readConfig(const TString& _strFileName)
{
  if (_strFileName.Length() == 0)
  {
    return;
  }

  TEnv env(_strFileName);

  Nhits_cut = env.GetValue("Nhits_cut", 0);

  maxpt = env.GetValue("maxpt", 0.);
  minpt = env.GetValue("minpt", 0.);
  maxptRF = env.GetValue("maxptRF", 0.);
  minptRF = env.GetValue("minptRF", 0.);
  eta_cut = env.GetValue("eta_cut", 0.);
  eta_gap = env.GetValue("eta_gap", 0.);
  DCAcut = env.GetValue("DCAcut", 0.);
  pid_probability = env.GetValue("pid_probability", 0.);

  debug = env.GetValue("debug", 0);
  Nevents = env.GetValue("Nevents", 0);
  
  readMCTracks = env.GetValue("readMCTracks", 0);
  harmonic = env.GetValue("harmonic", 0);
  bMotherIDcut = env.GetValue("bMotherIDcut", 0);

}

Bool_t trackCut(PicoDstRecoTrack *const &recoTrack, TF2 *const &fDCAx, TF2 *const &fDCAy, TF2 *const &fDCAz)
{
  if (!recoTrack) { return false; }
  Double_t pt = recoTrack->GetPt();
  Double_t eta = recoTrack->GetEta();
  if (pt < minpt || pt > maxpt || fabs(eta) > eta_cut)     return false;
  // if (fabs(recoTrack->GetDCAx()) > DCAcut)        return false; // static DCAx cut
  // if (fabs(recoTrack->GetDCAy()) > DCAcut)        return false; // static DCAy cut
  // if (fabs(recoTrack->GetDCAz()) > DCAcut)        return false; // static DCAz cut
  if (recoTrack->GetNhits() < Nhits_cut)          return false; // TPC hits cut    
  if (fabs(recoTrack->GetDCAx()) > fDCAx->Eval(recoTrack->GetPt(),recoTrack->GetEta())*DCAcut) return false; // DCAx cut
  if (fabs(recoTrack->GetDCAy()) > fDCAy->Eval(recoTrack->GetPt(),recoTrack->GetEta())*DCAcut) return false; // DCAy cut
  if (fabs(recoTrack->GetDCAz()) > fDCAz->Eval(recoTrack->GetPt(),recoTrack->GetEta())*DCAcut) return false; // DCAz cut
  return true;                                   
}

Bool_t trackCut(PicoDstMCTrack *const &mcTrack)
{
  if (!mcTrack) { return false; }
  if (format == "picodst" && mcTrack->GetMotherId() != -1) { return false; }
  Double_t pt = mcTrack->GetPt();
  Double_t eta = mcTrack->GetEta();
  auto particle = (TParticlePDG*) TDatabasePDG::Instance()->GetParticle(mcTrack->GetPdg());
  if (!particle) { return false; }
  Double_t charge = 1./3.*particle->Charge();
  if (pt < minpt || pt > maxpt || fabs(eta) > eta_cut || charge == 0) { return false; }
  return true;
}

Bool_t trackCutMotherID(PicoDstRecoTrack *const &recoTrack, PicoDstMCTrack *const &mcTrack)
{
  if (!recoTrack) { return false; }
  if (!mcTrack) { return false; }
  if (mcTrack->GetMotherId() != -1) { return false; }
  Double_t pt = recoTrack->GetPt();
  Double_t eta = recoTrack->GetEta();
  if (pt < minpt || pt > maxpt || fabs(eta) > eta_cut)     return false;
  if (recoTrack->GetNhits() < Nhits_cut)                   return false; // TPC hits cut    
  return true;                                   
}

void RunFlattening(TString inputFileName, TString outputFileName, TString configFileName = "")
{
  TStopwatch timer;
  timer.Start();

  if (configFileName.Length() > 0)
  {
    readConfig(configFileName);
  }

  if (debug)
  {
    cout << "Nevents = " << Nevents << endl;
    cout << "Nhits_cut = " << Nhits_cut << endl;

    cout << "maxpt = " << maxpt << endl;
    cout << "minpt = " << minpt << endl;
    cout << "maxptRF = " << maxptRF << endl;
    cout << "minptRF = " << minptRF << endl;
    cout << "eta_cut = " << eta_cut << endl;
    cout << "eta_gap = " << eta_gap << endl;
    cout << "DCAcut = " << DCAcut << endl;
    cout << "pid_probability = " << pid_probability << endl;
    cout << "readMCTracks = " << readMCTracks << endl;
    cout << "harmonic = " << harmonic << endl;       
    cout << "readMCTracks = " << readMCTracks << endl;
    cout << "harmonic = " << harmonic << endl;
    cout << "bMotherIDcut = " << bMotherIDcut << endl;
  }

  // Configure input information
  format = GetTreeName(inputFileName);
  if (debug) { cout << "format = " << format << endl; }
  TChain *chain = initChain(inputFileName, format.c_str());

  PicoDstMCEvent *mcEvent = nullptr;

  IReader* reader = nullptr;
  if (format == "picodst")
  {
    reader = new PicoDstReader();
  }
  if (format == "mctree")
  {
    reader = new McPicoReader();
    readMCTracks = true;
  }
  if (format == "tree")
  {
    reader = new ToyModelReader();
    readMCTracks = true;
  }

  if (!reader)
  {
    cerr << "No valid format is set!" << endl;
    return;
  }

  TFile *inputDCAfile;
  TF2 *fDCAx, *fDCAy, *fDCAz;
  if (!readMCTracks && !bMotherIDcut)
  { // using DCA cuts for primary track cut
    inputDCAfile = new TFile("DCA_FIT.root","read");
    if (!inputDCAfile) 
    {
      cerr << "Cannot find DCA_FIT.root file for DCA cut. Make sure you have run the DCA correction procedure and placed DCA_FIT.root in the executable directory." << endl;
      return;
    }
    fDCAx = dynamic_cast<TF2*> (inputDCAfile->Get("f_sigma0"));
    fDCAy = dynamic_cast<TF2*> (inputDCAfile->Get("f_sigma1"));
    fDCAz = dynamic_cast<TF2*> (inputDCAfile->Get("f_sigma2"));
    if (!fDCAx || !fDCAy || !fDCAz) { cerr << "Cannot find fit function for DCA primary track cuts!" << endl; return; }
    else{ cout << "Using " << DCAcut <<" sigma of DCA distr. cut." << endl; }
  }
  TFile *inputRecentering = new TFile("Recentering.root","READ");
  if (!inputRecentering) 
  {
    cerr << "Cannot find 'Recentering.root' file. Make sure you have run the recentering procedure beforehand and placed 'Recentering.root' in the executable directory." << endl;
    return;
  }
  TProfile *pQx[neta];            // <Qx>, neta = eta+, eta-
  TProfile *pQy[neta];            // <Qy>, neta = eta+, eta-
  Double_t qxMean[ncent][neta];   // <Qx>
  Double_t qyMean[ncent][neta];   // <Qy>
  for(Int_t ieta=0; ieta<neta; ieta++)
  {
    pQx[ieta] = dynamic_cast<TProfile*> (inputRecentering->Get(Form("pQx_%i",ieta)));
    pQy[ieta] = dynamic_cast<TProfile*> (inputRecentering->Get(Form("pQy_%i",ieta)));
    if (!pQx[ieta] || !pQy[ieta]) { cerr << "Cannot find needed TProfile for flattening procedure. Exitting..." << endl; return; }
    for (Int_t icent=0; icent<ncent; icent++)
    {
      qxMean[icent][ieta] = pQx[ieta]->GetBinContent(icent+1);
      qyMean[icent][ieta] = pQy[ieta]->GetBinContent(icent+1);
    }
  }
  reader->Init(chain);

  // Configure output information
  TFile *fo = new TFile(outputFileName.Data(), "recreate");
  TProfile *pCosNPsi[neta][nharm];
  TProfile *pSinNPsi[neta][nharm];
  for (Int_t ieta(0); ieta<neta; ieta++)
  {
    for (Int_t iharm(0); iharm<nharm; iharm++){
      pSinNPsi[ieta][iharm] = new TProfile(Form("pSinNPsi_%i_%i_rec", ieta,iharm), Form("<sin(2n#Psi_{2#eta%i})> recentered, harm_%i", ieta,iharm), ncent, 0., ncent);
      pCosNPsi[ieta][iharm] = new TProfile(Form("pCosNPsi_%i_%i_rec", ieta,iharm), Form("<cos(2n#Psi_{2#eta%i})> recentered, harm_%i", ieta,iharm), ncent, 0., ncent);
    }
  }

  Int_t icent, mult;
  Double_t bimp, cent, pt, eta, phi;
  QVector *qVec[neta];
  for (Int_t ieta=0; ieta<neta; ieta++)
  {
    qVec[ieta] = new QVector(harmonic);
  }
  Double_t qxRec, qyRec, psi2;
  Long64_t chain_size = chain->GetEntries();
  Long64_t n_entries = (Nevents < chain_size && Nevents > 0) ? Nevents : chain_size;
  cout << "Hi Master, let's do some physics together..." << endl;
  for (Int_t iEv = 0; iEv < n_entries; iEv++)
  {
    if (iEv % 10000 == 0)
      std::cout << "Event [" << iEv << "/" << n_entries << "]" << std::endl;
    // chain->GetEntry(iEv);
    mcEvent = reader->ReadMcEvent(iEv);
 
    // Read MC event
    bimp = mcEvent->GetB();
    cent = CentB(bimp);
    if (cent == -1)
      continue;
    icent = GetCentBin(cent);
    for (Int_t ieta=0; ieta<neta; ieta++) qVec[ieta]->Zero();
    
    if (readMCTracks) mult = reader->GetMcTrackSize();
    else mult = reader->GetRecoTrackSize();

    for (Int_t iTrk = 0; iTrk < mult; iTrk++)
    { // Track loop
      if (readMCTracks)
      {
        auto mcTrack = (PicoDstMCTrack *) reader->ReadMcTrack(iTrk);
        if (!trackCut(mcTrack)) { continue; } // TPC cut
        pt = mcTrack->GetPt();
        eta = mcTrack->GetEta();
        phi = mcTrack->GetPhi();
      }
      else
      { // Read reco tracks
        auto recoTrack = (PicoDstRecoTrack *) reader->ReadRecoTrack(iTrk);
        auto mcTrack = (PicoDstMCTrack *) reader->ReadMcTrack(recoTrack->GetMcId());
        if (bMotherIDcut) { if (!trackCutMotherID(recoTrack, mcTrack)) { continue; } }
        else { if (!trackCut(recoTrack,fDCAx,fDCAy,fDCAz)) { continue; } }
        pt = recoTrack->GetPt();
        eta = recoTrack->GetEta();
        phi = recoTrack->GetPhi();
      }
      // ipt = findBin(pt);
      if (pt > minptRF && pt < maxptRF)
      { // Reference Flow pt cut
        if (eta <-eta_gap) qVec[0]->CalQVector(phi, pt);
        if (eta > eta_gap) qVec[1]->CalQVector(phi, pt);
      }
    } // end of Track loop
    if (qVec[0]->GetMult() != 0 && qVec[1]->GetMult() != 0)
    {
      for (Int_t ieta(0); ieta<neta; ieta++)
      {
        qxRec = qVec[ieta]->X() - qxMean[icent][ieta];
        qyRec = qVec[ieta]->Y() - qyMean[icent][ieta];
        psi2 = TMath::ATan2(qyRec, qxRec) / 2.;

        // Fill Flattening profiles
        for (Int_t iharm=0;iharm<nharm;iharm++){
          pCosNPsi[ieta][iharm]->Fill(icent, TMath::Cos( 2*(iharm+1)*psi2));
          pSinNPsi[ieta][iharm]->Fill(icent, TMath::Sin( 2*(iharm+1)*psi2));
        }
      }
    }    
    
  } // end event loop
  
  // Writing output
  fo->cd();
  for(Int_t ieta=0; ieta<neta; ieta++ ){
    pQy[ieta]->Write();
    pQx[ieta]->Write();
    for (Int_t iharm=0;iharm<nharm;iharm++){
      pCosNPsi[ieta][iharm]->Write();
      pSinNPsi[ieta][iharm]->Write();
    }
  }
  fo->Close();

  timer.Stop();
  timer.Print();
}

int main(int argc, char **argv)
{
  TString iFileName, oFileName, configFileName = "";

  if (argc < 5)
  {
    std::cerr << "./FlowQCumulant -i INPUT -o OUTPUT [OPTIONAL: -config qcumulant.cfg]" << std::endl;
    return 1;
  }
  for (Int_t i = 1; i < argc; i++)
  {
    if (std::string(argv[i]) != "-i" &&
        std::string(argv[i]) != "-o" &&
        std::string(argv[i]) != "-config")
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
      if (std::string(argv[i]) == "-config" && i != argc - 1)
      {
        configFileName = argv[++i];
        continue;
      }
      if (std::string(argv[i]) == "-config" && i == argc - 1)
      {
        std::cerr << "\n[ERROR]: Output file name was not specified " << std::endl;
        return 1;
      }
    }
  }
  RunFlattening(iFileName, oFileName, configFileName);

  return 0;
}