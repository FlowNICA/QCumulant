#ifndef TOYMODEL_READER_H
#define TOYMODEL_READER_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <IReader.h>
#include <PicoDstMCEvent.h>
#include <PicoDstRecoEvent.h>
#include <PicoDstMCTrack.h>
#include <PicoDstRecoTrack.h>
#include <PicoDstFHCal.h>
// Header file for the classes stored in the TTree if any.

class ToyModelReader : virtual public IReader
{
public:
  TChain *fChain;  //!pointer to the analyzed TTree or TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  Int_t nh;
  Float_t b;
  Float_t rp;
  Float_t phi0[6000]; //[nh]
  Bool_t bFlow[6000]; //[nh]
  Float_t eta0[6000];  //[nh]
  Float_t pt0[6000];   //[nh]

  // List of branches
  TBranch *b_nh;    //!
  TBranch *b_b;     //!
  TBranch *b_rp;    //!
  TBranch *b_phi0;  //!
  TBranch *b_bFlow; //!
  TBranch *b_eta;   //!
  TBranch *b_pt;    //!

  ToyModelReader();
  virtual ~ToyModelReader();

  virtual void Init(TChain *chain);
  virtual PicoDstMCEvent *ReadMcEvent(Int_t ev_num);
  virtual PicoDstRecoEvent* ReadRecoEvent(Int_t ev_num) { return nullptr; }
  virtual Int_t GetMcTrackSize();
  virtual Int_t GetRecoTrackSize() { return 0; }
  virtual Int_t GetNFHCalModules() { return 0; }
  virtual PicoDstMCTrack *ReadMcTrack(Int_t tr_num);
  virtual PicoDstRecoTrack* ReadRecoTrack(Int_t tr_num) { return nullptr; }
  virtual PicoDstFHCal* ReadFHCalModule(Int_t module_num) { return nullptr; }

  // ClassDef(ToyModelReader,0);
};

#endif