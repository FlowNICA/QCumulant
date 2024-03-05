#include "ToyModelReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// ClassImp(ToyModelReader);
ToyModelReader::ToyModelReader() : fChain(0)
{
}

ToyModelReader::~ToyModelReader()
{
}

void ToyModelReader::Init(TChain *chain)
{
  if (!chain)
    return;
  fChain = chain;
  // fChain->SetMakeClass(1);

  fChain->SetBranchAddress("nh", &nh, &b_nh);
  fChain->SetBranchAddress("b", &b, &b_b);
  fChain->SetBranchAddress("rp", &rp, &b_rp);
  fChain->SetBranchAddress("phi0", phi0, &b_phi0);
  fChain->SetBranchAddress("bFlow", bFlow, &b_bFlow);
  fChain->SetBranchAddress("eta", eta0, &b_eta);
  fChain->SetBranchAddress("pt", pt0, &b_pt);
}

PicoDstMCEvent *ToyModelReader::ReadMcEvent(Int_t ev_num)
{
  if (!fChain)
    return nullptr;
  fChain->GetEntry(ev_num);
  auto event = new PicoDstMCEvent();
  event->SetB(b);
  event->SetPhiRP(rp);
  // event->SetVertex(0., 0., 0.);
  return event;
}

Int_t ToyModelReader::GetMcTrackSize()
{
  return nh;
}

PicoDstMCTrack *ToyModelReader::ReadMcTrack(Int_t tr_num)
{
  if (!fChain)
  {
    return nullptr;
  }
  auto track = new PicoDstMCTrack();
  track->SetPtEtaPhi(pt0[tr_num],eta0[tr_num],phi0[tr_num]);
  track->SetPdg(211); // pions
  track->SetEnergy(0.);
  return track;
}
