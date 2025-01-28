#include "PhysicsTools/XGBoost/interface/XGBooster.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"
#include "TFile.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TKey.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include <math.h>
#include "TMath.h"
#include <limits>
#include "TSystem.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "Math/GenVector/LorentzVector.h"

using namespace std;

void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], pat::XGBooster booster);
void runBDTeval(pat::XGBooster booster, float in_Jet_pt, float in_mT_b1MET, float in_d_zeta, float in_mutau_pt, float in_bmt_dR, float& BDTval);
void CopyDir(TDirectory *source, optutl::CommandLineParser parser);
void CopyFile(const char *fname, optutl::CommandLineParser parser);
void copyFiles(optutl::CommandLineParser parser, TFile* fOld, TFile* fNew);

float compute_deltaPhi(float v1, float v2) {
  float r = std::fmod(v2 - v1, 2*M_PI);
  if (r < -M_PI) {
    r += 2*M_PI;
  }
  else if (r >= M_PI) {
    r -= 2*M_PI;
  }
  return r;
}

float compute_deltaR(float e1, float e2, float p1, float p2) {
  float deta = e1 - e2;
  float dphi = compute_deltaPhi(p1, p2);
  return std::sqrt(deta*deta + dphi*dphi);
}

float compute_mt(float pt_1, float phi_1, float pt_met, float phi_met) {
  const auto dphi = compute_deltaPhi(phi_1, phi_met);
  return std::sqrt(2.0 * pt_1 * pt_met * (1.0 - std::cos(dphi)));
}

int main (int argc, char** argv) {
  optutl::CommandLineParser parser ("Sets Event Weights in the ntuple");
  parser.addOption("newFile",optutl::CommandLineParser::kString,"newFile","newFile.root");
  parser.addOption("inputFile",optutl::CommandLineParser::kString,"input File");
  parser.addOption("newOutputFile",optutl::CommandLineParser::kDouble,"New Output File",0.0);

  parser.parseArguments (argc, argv);

  pat::XGBooster booster("Cascade_1bjet_mu-tau.model");
  std::cout<<".model loaded"<<std::endl;
  booster.addFeature("Jet_pt");
  booster.addFeature("mT_b1MET");
  booster.addFeature("d_zeta");
  booster.addFeature("mutau_pt");
  booster.addFeature("bmt_dR");

  char TreeToUse[80]="first" ;

  TFile *fProduce;

  TFile *f = TFile::Open(parser.stringValue("inputFile").c_str());
  std::cout<<"Creating new outputfile"<<std::endl;
  std::string newFileName = parser.stringValue("newFile");

  fProduce = new TFile(newFileName.c_str(),"RECREATE");
  copyFiles(parser, f, fProduce);
  fProduce = new TFile(newFileName.c_str(),"UPDATE");
  std::cout<<"listing the directories================="<<std::endl;
  fProduce->ls();
  readdir(fProduce,parser,TreeToUse,booster);
  fProduce->Close();
  f->Close();
}

void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], pat::XGBooster booster) {
  TDirectory *dirsav = gDirectory;
  TIter next(dir->GetListOfKeys());
  TKey *key;
  char stringA[80]="first";
  dir->cd();
  int k=0;
  bool got_mutau_tree(false);
  bool got_etau_tree(false);
  bool got_emu_tree(false);
  while ((key = (TKey*)next())) {
    printf("Found key=%s \n",key->GetName());

    TObject *obj = key->ReadObj();
    if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
      dir->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      sprintf(TreeToUse,"%s",key->GetName());
      if (std::string(key->GetName()).find("_tree")) {
         readdir(subdir,parser,TreeToUse,booster);
      }
      dirsav->cd();
    }
    // else if(k<1 && obj->IsA()->InheritsFrom(TTree::Class())) {
    // Check all items in the tree: do not do k<1 criteria
    else if (obj->IsA()->InheritsFrom(TTree::Class())) {
      k++;
      sprintf(TreeToUse,"%s",obj->GetName());
      if (got_mutau_tree && std::string(key->GetName()).find("mutau") != std::string::npos) continue;
      if (got_etau_tree && std::string(key->GetName()).find("etau") != std::string::npos) continue;
      if (got_emu_tree && std::string(key->GetName()).find("emu") != std::string::npos) continue;
      std::cout<<"ici!!!!"<<std::endl;
      TTree *t = (TTree*)obj;
      // add the new branch
      float BDTval_nominal = -10;
      TBranch *newBranch1 = t->Branch("bdtscore", &BDTval_nominal, "bdtscore/F");
      // read the branches for input to BDT
      float Jet_pt, mT_b1MET, d_zeta, mutau_pt, bmt_dR;
      Float_t pt_1_nominal;
      Float_t eta_1;
      Float_t phi_1;
      Float_t m_1;
      Float_t pt_2_nominal;
      Float_t eta_2;
      Float_t phi_2;
      Float_t m_2;
      Float_t bpt_deepflavour_1;
      Float_t beta_deepflavour_1;
      Float_t bphi_deepflavour_1;
      Float_t bm_deepflavour_1;
      Float_t met;
      Float_t metphi;
      Float_t D_zeta_nominal;
      t->SetBranchAddress("pt_1_nominal", &pt_1_nominal);
      t->SetBranchAddress("eta_1", &eta_1);
      t->SetBranchAddress("phi_1", &phi_1);
      t->SetBranchAddress("m_1", &m_1);
      t->SetBranchAddress("pt_2_nominal", &pt_2_nominal);
      t->SetBranchAddress("eta_2", &eta_2);
      t->SetBranchAddress("phi_2", &phi_2);
      t->SetBranchAddress("m_2", &m_2);
      t->SetBranchAddress("bpt_deepflavour_1", &bpt_deepflavour_1);
      t->SetBranchAddress("beta_deepflavour_1", &beta_deepflavour_1);
      t->SetBranchAddress("bphi_deepflavour_1", &bphi_deepflavour_1);
      t->SetBranchAddress("bm_deepflavour_1", &bm_deepflavour_1);
      t->SetBranchAddress("met", &met);
      t->SetBranchAddress("metphi", &metphi);
      t->SetBranchAddress("D_zeta_nominal", &D_zeta_nominal);

      printf("Found tree -> weighting\n");
      for (Int_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
	Jet_pt = bpt_deepflavour_1;
	mT_b1MET = compute_mt(bpt_deepflavour_1, bphi_deepflavour_1, met, metphi);
	d_zeta = D_zeta_nominal;
	ROOT::Math::PtEtaPhiMVector bjet1(bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1);
	ROOT::Math::PtEtaPhiMVector tau1(pt_1_nominal, eta_1, phi_1, m_1);
	ROOT::Math::PtEtaPhiMVector tau2(pt_2_nominal, eta_2, phi_2, m_2);
        ROOT::Math::PtEtaPhiMVector tau12 = tau1 + tau2;
        mutau_pt = tau12.Pt();
        bmt_dR = compute_deltaR(beta_deepflavour_1, tau12.Eta(), bphi_deepflavour_1, tau12.Phi());
        if(Jet_pt > 0) runBDTeval(booster, Jet_pt, mT_b1MET, d_zeta, mutau_pt, bmt_dR, BDTval_nominal); 
        newBranch1->Fill();
      }
      dir->cd();
      t->Write("",TObject::kOverwrite);
      strcpy(TreeToUse,stringA);
      if (std::string(t->GetName()).find("mutau") != std::string::npos) got_mutau_tree = true;
      if (std::string(t->GetName()).find("etau") != std::string::npos) got_etau_tree = true;
      if (std::string(t->GetName()).find("emu") != std::string::npos) got_emu_tree = true;
    }
  }
}

void runBDTeval(pat::XGBooster booster, float in_Jet_pt, float in_mT_b1MET, float in_d_zeta, float in_mutau_pt, float in_bmt_dR, float& BDTval){
  booster.set("Jet_pt", in_Jet_pt);
  booster.set("mT_b1MET", in_mT_b1MET);
  booster.set("d_zeta", in_d_zeta);
  booster.set("mutau_pt", in_mutau_pt);
  booster.set("bmt_dR", in_bmt_dR);
  BDTval = booster.predict();
}

void CopyDir(TDirectory *source, optutl::CommandLineParser parser) {
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir;
  if(source->GetName()!=parser.stringValue("inputFile")){
    adir = savdir->mkdir(source->GetName());
    std::cout<<"Source name is not outputfile name"<<std::endl;
    adir->cd();
  }
  else{
    adir->cd();
  }

  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    std::cout<<"My key is: "<<key->GetName()<<std::endl;
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())) {
      source->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      adir->cd();
      CopyDir(subdir,parser);
      adir->cd();
    } else if (cl->InheritsFrom(TTree::Class())) {
      TTree *T = (TTree*)source->Get(key->GetName());
      adir->cd();
      TTree *newT = T->CloneTree(-1,"fast");
      newT->Write();
    } else {
      source->cd();
      TObject *obj = key->ReadObj();
      adir->cd();
      obj->Write();
      delete obj;
    }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}

void CopyFile(const char *fname, optutl::CommandLineParser parser) {
  //Copy all objects and subdirs of file fname as a subdir of the current directory
  TDirectory *target = gDirectory;
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    printf("Cannot copy file: %s\n",fname);
    target->cd();
    return;
  }
  target->cd();
  CopyDir(f,parser);
  delete f;
  target->cd();
}

void copyFiles(optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) {
  //prepare files to be copied
  if(gSystem->AccessPathName(parser.stringValue("inputFile").c_str())) {
    gSystem->CopyFile("hsimple.root", parser.stringValue("inputFile").c_str());
  }

  fNew->cd();
  CopyFile(parser.stringValue("inputFile").c_str(),parser);
  fNew->ls();
  fNew->Close();
}

