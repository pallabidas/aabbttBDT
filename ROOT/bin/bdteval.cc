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

void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], pat::XGBooster boosters[]);
void runBDT_1bjet(int booster_idx, pat::XGBooster booster, float jet_pt, float jet_eta, float jet_phi, float jet_mass, float pt_1, float eta_1, float phi_1, float mass_1, float pt_2, float eta_2, float phi_2, float mass_2, float met, float met_phi, float d_zeta, float& BDTval);
void runBDT_2bjet(int booster_idx, pat::XGBooster booster, float jet_1_pt, float jet_1_eta, float jet_1_phi, float jet_1_mass, float jet_2_pt, float jet_2_eta, float jet_2_phi, float jet_2_mass, float pt_1, float eta_1, float phi_1, float mass_1, float pt_2, float eta_2, float phi_2, float mass_2, float met, float met_phi, float d_zeta, float& BDTval);
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

  pat::XGBooster mutau_1bjet_booster("mutau-1bjet.json");
  std::cout<<"mutau-1bjet.json loaded"<<std::endl;
  mutau_1bjet_booster.addFeature("mT_thMET");
  mutau_1bjet_booster.addFeature("mT_mMET");
  mutau_1bjet_booster.addFeature("d_zeta");
  mutau_1bjet_booster.addFeature("mutau_pt");
  mutau_1bjet_booster.addFeature("m_b1mt");
  mutau_1bjet_booster.addFeature("bmt_dR");
  mutau_1bjet_booster.addFeature("b1th_dR");
  mutau_1bjet_booster.addFeature("njets");

  pat::XGBooster etau_1bjet_booster("etau-1bjet.json");
  std::cout<<"etau-1bjet.json loaded"<<std::endl;
  etau_1bjet_booster.addFeature("pt_vis_nominal");
  etau_1bjet_booster.addFeature("D_zeta_nominal");
  etau_1bjet_booster.addFeature("b1e_dR");
  etau_1bjet_booster.addFeature("m_btautau_vis_nominal");
  etau_1bjet_booster.addFeature("mtMET_1_nominal");
  etau_1bjet_booster.addFeature("b1th_dR");
  etau_1bjet_booster.addFeature("mtMET_2_nominal");
  etau_1bjet_booster.addFeature("njets");

  pat::XGBooster emu_1bjet_booster("emu-1bjet.json");
  std::cout<<"emu-1bjet.json loaded"<<std::endl;
  emu_1bjet_booster.addFeature("pt_vis_nominal");
  emu_1bjet_booster.addFeature("D_zeta_nominal");
  emu_1bjet_booster.addFeature("mT_b1MET");
  emu_1bjet_booster.addFeature("pt_1_nominal");
  emu_1bjet_booster.addFeature("m_btautau_vis_nominal");
  emu_1bjet_booster.addFeature("mtMET_1_nominal");
  emu_1bjet_booster.addFeature("b1emu_dR");
  emu_1bjet_booster.addFeature("njets");

  pat::XGBooster mutau_2bjet_booster("mutau_morethan1b.json");
  std::cout<<"mutau_morethan1b.json loaded"<<std::endl;
  mutau_2bjet_booster.addFeature("m_b2mt");
  mutau_2bjet_booster.addFeature("bmt_dR");
  mutau_2bjet_booster.addFeature("b2th_dR");
  mutau_2bjet_booster.addFeature("mT_mMET");
  mutau_2bjet_booster.addFeature("m_bbmt");
  mutau_2bjet_booster.addFeature("d_ma");

  pat::XGBooster etau_2bjet_booster("etau-morethan1b.json");
  std::cout<<"etau-morethan1b.json loaded"<<std::endl;
  etau_2bjet_booster.addFeature("mbb");
  etau_2bjet_booster.addFeature("m_btautau_vis_nominal");
  etau_2bjet_booster.addFeature("mtMET_1_nominal");
  etau_2bjet_booster.addFeature("b1th_dR");
  etau_2bjet_booster.addFeature("b1e_dR");
  etau_2bjet_booster.addFeature("b2th_dR");

  pat::XGBooster emu_2bjet_booster("emu-morethan1b.json");
  std::cout<<"emu-morethan1b.json loaded"<<std::endl;
  emu_2bjet_booster.addFeature("pt_vis_nominal");
  emu_2bjet_booster.addFeature("mT_b1MET");
  emu_2bjet_booster.addFeature("mtMET_1_nominal");
  emu_2bjet_booster.addFeature("b2emu_dR");
  emu_2bjet_booster.addFeature("d_ma");

  pat::XGBooster boosters[6] = {mutau_1bjet_booster, etau_1bjet_booster, emu_1bjet_booster, mutau_2bjet_booster, etau_2bjet_booster, emu_2bjet_booster};

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
  readdir(fProduce,parser,TreeToUse,boosters);
  fProduce->Close();
  f->Close();
}

void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], pat::XGBooster boosters[]) {
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
         readdir(subdir,parser,TreeToUse,boosters);
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
      TTree *t = (TTree*)obj;

      // add the new branches
      float BDTval_nominal = -10;
      float BDTval_es1Up = -10;
      float BDTval_es1Down = -10;
      float BDTval_es2Up = -10;
      float BDTval_es2Down = -10;
      float BDTval_UESUp = -10;
      float BDTval_UESDown = -10;
      float BDTval_ResponseUp = -10;
      float BDTval_ResponseDown = -10;
      float BDTval_ResolutionUp = -10;
      float BDTval_ResolutionDown = -10;
      float BDTval_JetAbsoluteUp = -10;
      float BDTval_JetAbsoluteDown = -10;
      float BDTval_JetAbsoluteyearUp = -10;
      float BDTval_JetAbsoluteyearDown = -10;
      float BDTval_JetBBEC1Up = -10;
      float BDTval_JetBBEC1Down = -10;
      float BDTval_JetBBEC1yearUp = -10;
      float BDTval_JetBBEC1yearDown = -10;
      float BDTval_JetEC2Up = -10;
      float BDTval_JetEC2Down = -10;
      float BDTval_JetEC2yearUp = -10;
      float BDTval_JetEC2yearDown = -10;
      float BDTval_JetFlavorQCDUp = -10;
      float BDTval_JetFlavorQCDDown = -10;
      float BDTval_JetHFUp = -10;
      float BDTval_JetHFDown = -10;
      float BDTval_JetHFyearUp = -10;
      float BDTval_JetHFyearDown = -10;
      float BDTval_JetRelativeBalUp = -10;
      float BDTval_JetRelativeBalDown = -10;
      float BDTval_JetRelativeSampleUp = -10;
      float BDTval_JetRelativeSampleDown = -10;
      float BDTval_JERUp = -10;
      float BDTval_JERDown = -10;

      TBranch *newBranch1 = t->Branch("bdtscore", &BDTval_nominal, "bdtscore/F");
      TBranch *newBranch1_es1U = t->Branch("bdtscore_es1Up",   &BDTval_es1Up,   "bdtscore_es1Up/F");
      TBranch *newBranch1_es1D = t->Branch("bdtscore_es1Down",   &BDTval_es1Down,   "bdtscore_es1Down/F");
      TBranch *newBranch1_es2U = t->Branch("bdtscore_es2Up", &BDTval_es2Up, "bdtscore_es2Up/F");
      TBranch *newBranch1_es2D = t->Branch("bdtscore_es2Down",   &BDTval_es2Down,   "bdtscore_es2Down/F");
      TBranch *newBranch1UU = t->Branch("bdtscore_UESUp", &BDTval_UESUp, "bdtscore_UESUp/F");
      TBranch *newBranch1UD = t->Branch("bdtscore_UESDown", &BDTval_UESDown, "bdtscore_UESDown/F");
      TBranch *newBranch1ResponseU = t->Branch("bdtscore_responseUp", &BDTval_ResponseUp, "bdtscore_responseUp/F");
      TBranch *newBranch1ResponseD = t->Branch("bdtscore_responseDown", &BDTval_ResponseDown, "bdtscore_responseDown/F");
      TBranch *newBranch1ResolutionU = t->Branch("bdtscore_resolutionUp", &BDTval_ResolutionUp, "bdtscore_resolutionUp/F");
      TBranch *newBranch1ResolutionD = t->Branch("bdtscore_resolutionDown", &BDTval_ResolutionDown, "bdtscore_resolutionDown/F");
      TBranch *newBranch1JetAbsoluteU = t->Branch("bdtscore_JetAbsoluteUp", &BDTval_JetAbsoluteUp, "bdtscore_JetAbsoluteUp/F");
      TBranch *newBranch1JetAbsoluteD = t->Branch("bdtscore_JetAbsoluteDown", &BDTval_JetAbsoluteDown, "bdtscore_JetAbsoluteDown/F");
      TBranch *newBranch1JetAbsoluteyearU = t->Branch("bdtscore_JetAbsoluteyearUp", &BDTval_JetAbsoluteyearUp, "bdtscore_JetAbsoluteyearUp/F");
      TBranch *newBranch1JetAbsoluteyearD = t->Branch("bdtscore_JetAbsoluteyearDown", &BDTval_JetAbsoluteyearDown, "bdtscore_JetAbsoluteyearDown/F");
      TBranch *newBranch1JetBBEC1U = t->Branch("bdtscore_JetBBEC1Up", &BDTval_JetBBEC1Up, "bdtscore_JetBBEC1Up/F");
      TBranch *newBranch1JetBBEC1D = t->Branch("bdtscore_JetBBEC1Down", &BDTval_JetBBEC1Down, "bdtscore_JetBBEC1Down/F");
      TBranch *newBranch1JetBBEC1yearU = t->Branch("bdtscore_JetBBEC1yearUp", &BDTval_JetBBEC1yearUp, "bdtscore_JetBBEC1yearUp/F");
      TBranch *newBranch1JetBBEC1yearD = t->Branch("bdtscore_JetBBEC1yearDown", &BDTval_JetBBEC1yearDown, "bdtscore_JetBBEC1yearDown/F");
      TBranch *newBranch1JetEC2U = t->Branch("bdtscore_JetEC2Up", &BDTval_JetEC2Up, "bdtscore_JetEC2Up/F");
      TBranch *newBranch1JetEC2D = t->Branch("bdtscore_JetEC2Down", &BDTval_JetEC2Down, "bdtscore_JetEC2Down/F");
      TBranch *newBranch1JetEC2yearU = t->Branch("bdtscore_JetEC2yearUp", &BDTval_JetEC2yearUp, "bdtscore_JetEC2yearUp/F");
      TBranch *newBranch1JetEC2yearD = t->Branch("bdtscore_JetEC2yearDown", &BDTval_JetEC2yearDown, "bdtscore_JetEC2yearDown/F");
      TBranch *newBranch1JetFlavorQCDU = t->Branch("bdtscore_JetFlavorQCDUp", &BDTval_JetFlavorQCDUp, "bdtscore_JetFlavorQCDUp/F");
      TBranch *newBranch1JetFlavorQCDD = t->Branch("bdtscore_JetFlavorQCDDown", &BDTval_JetFlavorQCDDown, "bdtscore_JetFlavorQCDDown/F");
      TBranch *newBranch1JetHFU = t->Branch("bdtscore_JetHFUp", &BDTval_JetHFUp, "bdtscore_JetHFUp/F");
      TBranch *newBranch1JetHFD = t->Branch("bdtscore_JetHFDown", &BDTval_JetHFDown, "bdtscore_JetHFDown/F");
      TBranch *newBranch1JetHFyearU = t->Branch("bdtscore_JetHFyearUp", &BDTval_JetHFyearUp, "bdtscore_JetHFyearUp/F");
      TBranch *newBranch1JetHFyearD = t->Branch("bdtscore_JetHFyearDown", &BDTval_JetHFyearDown, "bdtscore_JetHFyearDown/F");
      TBranch *newBranch1JetRelativeBalU = t->Branch("bdtscore_JetRelativeBalUp", &BDTval_JetRelativeBalUp, "bdtscore_JetRelativeBalUp/F");
      TBranch *newBranch1JetRelativeBalD = t->Branch("bdtscore_JetRelativeBalDown", &BDTval_JetRelativeBalDown, "bdtscore_JetRelativeBalDown/F");
      TBranch *newBranch1JetRelativeSampleU = t->Branch("bdtscore_JetRelativeSampleUp", &BDTval_JetRelativeSampleUp, "bdtscore_JetRelativeSampleUp/F");
      TBranch *newBranch1JetRelativeSampleD = t->Branch("bdtscore_JetRelativeSampleDown", &BDTval_JetRelativeSampleDown, "bdtscore_JetRelativeSampleDown/F");
      TBranch *newBranch1JERU = t->Branch("bdtscore_JERUp", &BDTval_JERUp, "bdtscore_JERUp/F");
      TBranch *newBranch1JERD = t->Branch("bdtscore_JERDown", &BDTval_JERDown, "bdtscore_JERDown/F");

      // read the branches needed for input to BDT
      Int_t channel;
      Float_t pt_1_nominal;
      Float_t pt_1_es1Up;
      Float_t pt_1_es1Down;
      Float_t eta_1;
      Float_t phi_1;
      Float_t m_1;
      Float_t m_1_es1Up;
      Float_t m_1_es1Down;
      Float_t pt_2_nominal;
      Float_t pt_2_es2Up;
      Float_t pt_2_es2Down;
      Float_t eta_2;
      Float_t phi_2;
      Float_t m_2;
      Float_t m_2_es2Up;
      Float_t m_2_es2Down;
      Float_t bpt_deepflavour_1;
      Float_t bpt_deepflavour_JERUp_1;
      Float_t bpt_deepflavour_JERDown_1;
      Float_t bpt_deepflavour_JetAbsoluteUp_1;
      Float_t bpt_deepflavour_JetAbsoluteDown_1;
      Float_t bpt_deepflavour_JetAbsoluteyearUp_1;
      Float_t bpt_deepflavour_JetAbsoluteyearDown_1;
      Float_t bpt_deepflavour_JetBBEC1Up_1;
      Float_t bpt_deepflavour_JetBBEC1Down_1;
      Float_t bpt_deepflavour_JetBBEC1yearUp_1;
      Float_t bpt_deepflavour_JetBBEC1yearDown_1;
      Float_t bpt_deepflavour_JetEC2Up_1;
      Float_t bpt_deepflavour_JetEC2Down_1;
      Float_t bpt_deepflavour_JetEC2yearUp_1;
      Float_t bpt_deepflavour_JetEC2yearDown_1;
      Float_t bpt_deepflavour_JetFlavorQCDUp_1;
      Float_t bpt_deepflavour_JetFlavorQCDDown_1;
      Float_t bpt_deepflavour_JetHFUp_1;
      Float_t bpt_deepflavour_JetHFDown_1;
      Float_t bpt_deepflavour_JetHFyearUp_1;
      Float_t bpt_deepflavour_JetHFyearDown_1;
      Float_t bpt_deepflavour_JetRelativeBalUp_1;
      Float_t bpt_deepflavour_JetRelativeBalDown_1;
      Float_t bpt_deepflavour_JetRelativeSampleUp_1;
      Float_t bpt_deepflavour_JetRelativeSampleDown_1;
      Float_t beta_deepflavour_1;
      Float_t bphi_deepflavour_1;
      Float_t bm_deepflavour_1;
      Float_t bm_deepflavour_JERUp_1;
      Float_t bm_deepflavour_JERDown_1;
      Float_t bm_deepflavour_JetAbsoluteUp_1;
      Float_t bm_deepflavour_JetAbsoluteDown_1;
      Float_t bm_deepflavour_JetAbsoluteyearUp_1;
      Float_t bm_deepflavour_JetAbsoluteyearDown_1;
      Float_t bm_deepflavour_JetBBEC1Up_1;
      Float_t bm_deepflavour_JetBBEC1Down_1;
      Float_t bm_deepflavour_JetBBEC1yearUp_1;
      Float_t bm_deepflavour_JetBBEC1yearDown_1;
      Float_t bm_deepflavour_JetEC2Up_1;
      Float_t bm_deepflavour_JetEC2Down_1;
      Float_t bm_deepflavour_JetEC2yearUp_1;
      Float_t bm_deepflavour_JetEC2yearDown_1;
      Float_t bm_deepflavour_JetFlavorQCDUp_1;
      Float_t bm_deepflavour_JetFlavorQCDDown_1;
      Float_t bm_deepflavour_JetHFUp_1;
      Float_t bm_deepflavour_JetHFDown_1;
      Float_t bm_deepflavour_JetHFyearUp_1;
      Float_t bm_deepflavour_JetHFyearDown_1;
      Float_t bm_deepflavour_JetRelativeBalUp_1;
      Float_t bm_deepflavour_JetRelativeBalDown_1;
      Float_t bm_deepflavour_JetRelativeSampleUp_1;
      Float_t bm_deepflavour_JetRelativeSampleDown_1;
      Float_t met_nominal;
      Float_t met_es1Up;
      Float_t met_es1Down;
      Float_t met_es2Up;
      Float_t met_es2Down;
      Float_t met_UESUp;
      Float_t met_UESDown;
      Float_t met_ResponseUp;
      Float_t met_ResponseDown;
      Float_t met_ResolutionUp;
      Float_t met_ResolutionDown;
      Float_t met_JERUp;
      Float_t met_JERDown;
      Float_t met_JetAbsoluteUp;
      Float_t met_JetAbsoluteDown;
      Float_t met_JetAbsoluteyearUp;
      Float_t met_JetAbsoluteyearDown;
      Float_t met_JetBBEC1Up;
      Float_t met_JetBBEC1Down;
      Float_t met_JetBBEC1yearUp;
      Float_t met_JetBBEC1yearDown;
      Float_t met_JetEC2Up;
      Float_t met_JetEC2Down;
      Float_t met_JetEC2yearUp;
      Float_t met_JetEC2yearDown;
      Float_t met_JetFlavorQCDUp;
      Float_t met_JetFlavorQCDDown;
      Float_t met_JetHFUp;
      Float_t met_JetHFDown;
      Float_t met_JetHFyearUp;
      Float_t met_JetHFyearDown;
      Float_t met_JetRelativeBalUp;
      Float_t met_JetRelativeBalDown;
      Float_t met_JetRelativeSampleUp;
      Float_t met_JetRelativeSampleDown;
      Float_t metphi_nominal;
      Float_t metphi_es1Up;
      Float_t metphi_es1Down;
      Float_t metphi_es2Up;
      Float_t metphi_es2Down;
      Float_t metphi_UESUp;
      Float_t metphi_UESDown;
      Float_t metphi_ResponseUp;
      Float_t metphi_ResponseDown;
      Float_t metphi_ResolutionUp;
      Float_t metphi_ResolutionDown;
      Float_t metphi_JERUp;
      Float_t metphi_JERDown;
      Float_t metphi_JetAbsoluteUp;
      Float_t metphi_JetAbsoluteDown;
      Float_t metphi_JetAbsoluteyearUp;
      Float_t metphi_JetAbsoluteyearDown;
      Float_t metphi_JetBBEC1Up;
      Float_t metphi_JetBBEC1Down;
      Float_t metphi_JetBBEC1yearUp;
      Float_t metphi_JetBBEC1yearDown;
      Float_t metphi_JetEC2Up;
      Float_t metphi_JetEC2Down;
      Float_t metphi_JetEC2yearUp;
      Float_t metphi_JetEC2yearDown;
      Float_t metphi_JetFlavorQCDUp;
      Float_t metphi_JetFlavorQCDDown;
      Float_t metphi_JetHFUp;
      Float_t metphi_JetHFDown;
      Float_t metphi_JetHFyearUp;
      Float_t metphi_JetHFyearDown;
      Float_t metphi_JetRelativeBalUp;
      Float_t metphi_JetRelativeBalDown;
      Float_t metphi_JetRelativeSampleUp;
      Float_t metphi_JetRelativeSampleDown;
      Float_t D_zeta_nominal;
      Float_t D_zeta_es1Up;
      Float_t D_zeta_es1Down;
      Float_t D_zeta_es2Up;
      Float_t D_zeta_es2Down;
      Float_t D_zeta_UESUp;
      Float_t D_zeta_UESDown;
      Float_t D_zeta_ResponseUp;
      Float_t D_zeta_ResponseDown;
      Float_t D_zeta_ResolutionUp;
      Float_t D_zeta_ResolutionDown;
      Float_t D_zeta_JERUp;
      Float_t D_zeta_JERDown;
      Float_t D_zeta_JetAbsoluteUp;
      Float_t D_zeta_JetAbsoluteDown;
      Float_t D_zeta_JetAbsoluteyearUp;
      Float_t D_zeta_JetAbsoluteyearDown;
      Float_t D_zeta_JetBBEC1Up;
      Float_t D_zeta_JetBBEC1Down;
      Float_t D_zeta_JetBBEC1yearUp;
      Float_t D_zeta_JetBBEC1yearDown;
      Float_t D_zeta_JetEC2Up;
      Float_t D_zeta_JetEC2Down;
      Float_t D_zeta_JetEC2yearUp;
      Float_t D_zeta_JetEC2yearDown;
      Float_t D_zeta_JetFlavorQCDUp;
      Float_t D_zeta_JetFlavorQCDDown;
      Float_t D_zeta_JetHFUp;
      Float_t D_zeta_JetHFDown;
      Float_t D_zeta_JetHFyearUp;
      Float_t D_zeta_JetHFyearDown;
      Float_t D_zeta_JetRelativeBalUp;
      Float_t D_zeta_JetRelativeBalDown;
      Float_t D_zeta_JetRelativeSampleUp;
      Float_t D_zeta_JetRelativeSampleDown;

      // New Change
      Float_t bpt_deepflavour_2;
      Float_t bpt_deepflavour_JERUp_2;
      Float_t bpt_deepflavour_JERDown_2;
      Float_t bpt_deepflavour_JetAbsoluteUp_2;
      Float_t bpt_deepflavour_JetAbsoluteDown_2;
      Float_t bpt_deepflavour_JetAbsoluteyearUp_2;
      Float_t bpt_deepflavour_JetAbsoluteyearDown_2;
      Float_t bpt_deepflavour_JetBBEC1Up_2;
      Float_t bpt_deepflavour_JetBBEC1Down_2;
      Float_t bpt_deepflavour_JetBBEC1yearUp_2;
      Float_t bpt_deepflavour_JetBBEC1yearDown_2;
      Float_t bpt_deepflavour_JetEC2Up_2;
      Float_t bpt_deepflavour_JetEC2Down_2;
      Float_t bpt_deepflavour_JetEC2yearUp_2;
      Float_t bpt_deepflavour_JetEC2yearDown_2;
      Float_t bpt_deepflavour_JetFlavorQCDUp_2;
      Float_t bpt_deepflavour_JetFlavorQCDDown_2;
      Float_t bpt_deepflavour_JetHFUp_2;
      Float_t bpt_deepflavour_JetHFDown_2;
      Float_t bpt_deepflavour_JetHFyearUp_2;
      Float_t bpt_deepflavour_JetHFyearDown_2;
      Float_t bpt_deepflavour_JetRelativeBalUp_2;
      Float_t bpt_deepflavour_JetRelativeBalDown_2;
      Float_t bpt_deepflavour_JetRelativeSampleUp_2;
      Float_t bpt_deepflavour_JetRelativeSampleDown_2;
      Float_t beta_deepflavour_2;
      Float_t bphi_deepflavour_2;
      Float_t bm_deepflavour_2;
      Float_t bm_deepflavour_JERUp_2;
      Float_t bm_deepflavour_JERDown_2;
      Float_t bm_deepflavour_JetAbsoluteUp_2;
      Float_t bm_deepflavour_JetAbsoluteDown_2;
      Float_t bm_deepflavour_JetAbsoluteyearUp_2;
      Float_t bm_deepflavour_JetAbsoluteyearDown_2;
      Float_t bm_deepflavour_JetBBEC1Up_2;
      Float_t bm_deepflavour_JetBBEC1Down_2;
      Float_t bm_deepflavour_JetBBEC1yearUp_2;
      Float_t bm_deepflavour_JetBBEC1yearDown_2;
      Float_t bm_deepflavour_JetEC2Up_2;
      Float_t bm_deepflavour_JetEC2Down_2;
      Float_t bm_deepflavour_JetEC2yearUp_2;
      Float_t bm_deepflavour_JetEC2yearDown_2;
      Float_t bm_deepflavour_JetFlavorQCDUp_2;
      Float_t bm_deepflavour_JetFlavorQCDDown_2;
      Float_t bm_deepflavour_JetHFUp_2;
      Float_t bm_deepflavour_JetHFDown_2;
      Float_t bm_deepflavour_JetHFyearUp_2;
      Float_t bm_deepflavour_JetHFyearDown_2;
      Float_t bm_deepflavour_JetRelativeBalUp_2;
      Float_t bm_deepflavour_JetRelativeBalDown_2;
      Float_t bm_deepflavour_JetRelativeSampleUp_2;
      Float_t bm_deepflavour_JetRelativeSampleDown_2;

      t->SetBranchAddress("channel", &channel);
      t->SetBranchAddress("pt_1_nominal", &pt_1_nominal);
      t->SetBranchAddress("pt_1_es1Up",&pt_1_es1Up);
      t->SetBranchAddress("pt_1_es1Down",&pt_1_es1Down);
      t->SetBranchAddress("eta_1", &eta_1);
      t->SetBranchAddress("phi_1", &phi_1);
      t->SetBranchAddress("m_1", &m_1);
      t->SetBranchAddress("m_1_es1Up", &m_1_es1Up);
      t->SetBranchAddress("m_1_es1Down", &m_1_es1Down);
      t->SetBranchAddress("pt_2_nominal", &pt_2_nominal);
      t->SetBranchAddress("pt_2_es2Up", &pt_2_es2Up);
      t->SetBranchAddress("pt_2_es2Down",&pt_2_es2Down);
      t->SetBranchAddress("eta_2", &eta_2);
      t->SetBranchAddress("phi_2", &phi_2);
      t->SetBranchAddress("m_2", &m_2);
      t->SetBranchAddress("m_2_es2Up",  &m_2_es2Up);
      t->SetBranchAddress("m_2_es2Down",&m_2_es2Down);
      t->SetBranchAddress("bpt_deepflavour_1", &bpt_deepflavour_1);
      t->SetBranchAddress("bpt_deepflavour_JERUp_1", &bpt_deepflavour_JERUp_1);
      t->SetBranchAddress("bpt_deepflavour_JERDown_1", &bpt_deepflavour_JERDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteUp_1", &bpt_deepflavour_JetAbsoluteUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteDown_1", &bpt_deepflavour_JetAbsoluteDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteyearUp_1", &bpt_deepflavour_JetAbsoluteyearUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteyearDown_1", &bpt_deepflavour_JetAbsoluteyearDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1Up_1", &bpt_deepflavour_JetBBEC1Up_1);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1Down_1", &bpt_deepflavour_JetBBEC1Down_1);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1yearUp_1", &bpt_deepflavour_JetBBEC1yearUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1yearDown_1", &bpt_deepflavour_JetBBEC1yearDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetEC2Up_1", &bpt_deepflavour_JetEC2Up_1);
      t->SetBranchAddress("bpt_deepflavour_JetEC2Down_1", &bpt_deepflavour_JetEC2Down_1);
      t->SetBranchAddress("bpt_deepflavour_JetEC2yearUp_1", &bpt_deepflavour_JetEC2yearUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetEC2yearDown_1", &bpt_deepflavour_JetEC2yearDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetFlavorQCDUp_1", &bpt_deepflavour_JetFlavorQCDUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetFlavorQCDDown_1", &bpt_deepflavour_JetFlavorQCDDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetHFUp_1", &bpt_deepflavour_JetHFUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetHFDown_1", &bpt_deepflavour_JetHFDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetHFyearUp_1", &bpt_deepflavour_JetHFyearUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetHFyearDown_1", &bpt_deepflavour_JetHFyearDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeBalUp_1", &bpt_deepflavour_JetRelativeBalUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeBalDown_1", &bpt_deepflavour_JetRelativeBalDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeSampleUp_1", &bpt_deepflavour_JetRelativeSampleUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeSampleDown_1", &bpt_deepflavour_JetRelativeSampleDown_1);
      t->SetBranchAddress("beta_deepflavour_1", &beta_deepflavour_1);
      t->SetBranchAddress("bphi_deepflavour_1", &bphi_deepflavour_1);
      t->SetBranchAddress("bm_deepflavour_1", &bm_deepflavour_1);
      t->SetBranchAddress("bm_deepflavour_JERUp_1", &bm_deepflavour_JERUp_1);
      t->SetBranchAddress("bm_deepflavour_JERDown_1", &bm_deepflavour_JERDown_1);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteUp_1", &bm_deepflavour_JetAbsoluteUp_1);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteDown_1", &bm_deepflavour_JetAbsoluteDown_1);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteyearUp_1", &bm_deepflavour_JetAbsoluteyearUp_1);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteyearDown_1", &bm_deepflavour_JetAbsoluteyearDown_1);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1Up_1", &bm_deepflavour_JetBBEC1Up_1);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1Down_1", &bm_deepflavour_JetBBEC1Down_1);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1yearUp_1", &bm_deepflavour_JetBBEC1yearUp_1);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1yearDown_1", &bm_deepflavour_JetBBEC1yearDown_1);
      t->SetBranchAddress("bm_deepflavour_JetEC2Up_1", &bm_deepflavour_JetEC2Up_1);
      t->SetBranchAddress("bm_deepflavour_JetEC2Down_1", &bm_deepflavour_JetEC2Down_1);
      t->SetBranchAddress("bm_deepflavour_JetEC2yearUp_1", &bm_deepflavour_JetEC2yearUp_1);
      t->SetBranchAddress("bm_deepflavour_JetEC2yearDown_1", &bm_deepflavour_JetEC2yearDown_1);
      t->SetBranchAddress("bm_deepflavour_JetFlavorQCDUp_1", &bm_deepflavour_JetFlavorQCDUp_1);
      t->SetBranchAddress("bm_deepflavour_JetFlavorQCDDown_1", &bm_deepflavour_JetFlavorQCDDown_1);
      t->SetBranchAddress("bm_deepflavour_JetHFUp_1", &bm_deepflavour_JetHFUp_1);
      t->SetBranchAddress("bm_deepflavour_JetHFDown_1", &bm_deepflavour_JetHFDown_1);
      t->SetBranchAddress("bm_deepflavour_JetHFyearUp_1", &bm_deepflavour_JetHFyearUp_1);
      t->SetBranchAddress("bm_deepflavour_JetHFyearDown_1", &bm_deepflavour_JetHFyearDown_1);
      t->SetBranchAddress("bm_deepflavour_JetRelativeBalUp_1", &bm_deepflavour_JetRelativeBalUp_1);
      t->SetBranchAddress("bm_deepflavour_JetRelativeBalDown_1", &bm_deepflavour_JetRelativeBalDown_1);
      t->SetBranchAddress("bm_deepflavour_JetRelativeSampleUp_1", &bm_deepflavour_JetRelativeSampleUp_1);
      t->SetBranchAddress("bm_deepflavour_JetRelativeSampleDown_1", &bm_deepflavour_JetRelativeSampleDown_1);
      t->SetBranchAddress("met_nominal", &met_nominal);
      t->SetBranchAddress("met_es1Up", &met_es1Up);
      t->SetBranchAddress("met_es1Down", &met_es1Down);
      t->SetBranchAddress("met_es2Up", &met_es2Up);
      t->SetBranchAddress("met_es2Down", &met_es2Down);
      t->SetBranchAddress("met_UESUp", &met_UESUp);
      t->SetBranchAddress("met_UESDown", &met_UESDown);
      t->SetBranchAddress("met_responseUp", &met_ResponseUp);
      t->SetBranchAddress("met_responseDown", &met_ResponseDown);
      t->SetBranchAddress("met_resolutionUp", &met_ResolutionUp);
      t->SetBranchAddress("met_resolutionDown", &met_ResolutionDown);
      t->SetBranchAddress("met_JetAbsoluteUp", &met_JetAbsoluteUp);
      t->SetBranchAddress("met_JetAbsoluteDown", &met_JetAbsoluteDown);
      t->SetBranchAddress("met_JetAbsoluteyearUp", &met_JetAbsoluteyearUp);
      t->SetBranchAddress("met_JetAbsoluteyearDown", &met_JetAbsoluteyearDown);
      t->SetBranchAddress("met_JetBBEC1Up", &met_JetBBEC1Up);
      t->SetBranchAddress("met_JetBBEC1Down", &met_JetBBEC1Down);
      t->SetBranchAddress("met_JetBBEC1yearUp", &met_JetBBEC1yearUp);
      t->SetBranchAddress("met_JetBBEC1yearDown", &met_JetBBEC1yearDown);
      t->SetBranchAddress("met_JetEC2Up", &met_JetEC2Up);
      t->SetBranchAddress("met_JetEC2Down", &met_JetEC2Down);
      t->SetBranchAddress("met_JetEC2yearUp", &met_JetEC2yearUp);
      t->SetBranchAddress("met_JetEC2yearDown", &met_JetEC2yearDown);
      t->SetBranchAddress("met_JetFlavorQCDUp", &met_JetFlavorQCDUp);
      t->SetBranchAddress("met_JetFlavorQCDDown", &met_JetFlavorQCDDown);
      t->SetBranchAddress("met_JetHFUp", &met_JetHFUp);
      t->SetBranchAddress("met_JetHFDown", &met_JetHFDown);
      t->SetBranchAddress("met_JetHFyearUp", &met_JetHFyearUp);
      t->SetBranchAddress("met_JetHFyearDown", &met_JetHFyearDown);
      t->SetBranchAddress("met_JetRelativeBalUp", &met_JetRelativeBalUp);
      t->SetBranchAddress("met_JetRelativeBalDown", &met_JetRelativeBalDown);
      t->SetBranchAddress("met_JetRelativeSampleUp", &met_JetRelativeSampleUp);
      t->SetBranchAddress("met_JetRelativeSampleDown", &met_JetRelativeSampleDown);
      t->SetBranchAddress("met_JERUp", &met_JERUp);
      t->SetBranchAddress("met_JERDown", &met_JERDown);
      t->SetBranchAddress("metphi_nominal", &metphi_nominal);
      t->SetBranchAddress("metphi_es1Up", &metphi_es1Up);
      t->SetBranchAddress("metphi_es1Down",&metphi_es1Down);
      t->SetBranchAddress("metphi_es2Up", &metphi_es2Up);
      t->SetBranchAddress("metphi_es2Down", &metphi_es2Down);
      t->SetBranchAddress("metphi_UESUp", &metphi_UESUp);
      t->SetBranchAddress("metphi_UESDown", &metphi_UESDown);
      t->SetBranchAddress("metphi_responseUp", &metphi_ResponseUp);
      t->SetBranchAddress("metphi_responseDown", &metphi_ResponseDown);
      t->SetBranchAddress("metphi_resolutionUp", &metphi_ResolutionUp);
      t->SetBranchAddress("metphi_resolutionDown", &metphi_ResolutionDown);
      t->SetBranchAddress("metphi_JetAbsoluteUp", &metphi_JetAbsoluteUp);
      t->SetBranchAddress("metphi_JetAbsoluteDown", &metphi_JetAbsoluteDown);
      t->SetBranchAddress("metphi_JetAbsoluteyearUp", &metphi_JetAbsoluteyearUp);
      t->SetBranchAddress("metphi_JetAbsoluteyearDown", &metphi_JetAbsoluteyearDown);
      t->SetBranchAddress("metphi_JetBBEC1Up", &metphi_JetBBEC1Up);
      t->SetBranchAddress("metphi_JetBBEC1Down", &metphi_JetBBEC1Down);
      t->SetBranchAddress("metphi_JetBBEC1yearUp", &metphi_JetBBEC1yearUp);
      t->SetBranchAddress("metphi_JetBBEC1yearDown", &metphi_JetBBEC1yearDown);
      t->SetBranchAddress("metphi_JetEC2Up", &metphi_JetEC2Up);
      t->SetBranchAddress("metphi_JetEC2Down", &metphi_JetEC2Down);
      t->SetBranchAddress("metphi_JetEC2yearUp", &metphi_JetEC2yearUp);
      t->SetBranchAddress("metphi_JetEC2yearDown", &metphi_JetEC2yearDown);
      t->SetBranchAddress("metphi_JetFlavorQCDUp", &metphi_JetFlavorQCDUp);
      t->SetBranchAddress("metphi_JetFlavorQCDDown", &metphi_JetFlavorQCDDown);
      t->SetBranchAddress("metphi_JetHFUp", &metphi_JetHFUp);
      t->SetBranchAddress("metphi_JetHFDown", &metphi_JetHFDown);
      t->SetBranchAddress("metphi_JetHFyearUp", &metphi_JetHFyearUp);
      t->SetBranchAddress("metphi_JetHFyearDown", &metphi_JetHFyearDown);
      t->SetBranchAddress("metphi_JetRelativeBalUp", &metphi_JetRelativeBalUp);
      t->SetBranchAddress("metphi_JetRelativeBalDown", &metphi_JetRelativeBalDown);
      t->SetBranchAddress("metphi_JetRelativeSampleUp", &metphi_JetRelativeSampleUp);
      t->SetBranchAddress("metphi_JetRelativeSampleDown", &metphi_JetRelativeSampleDown);
      t->SetBranchAddress("metphi_JERUp", &metphi_JERUp);
      t->SetBranchAddress("metphi_JERDown", &metphi_JERDown);
      t->SetBranchAddress("D_zeta_nominal", &D_zeta_nominal);
      t->SetBranchAddress("D_zeta_es1Up", &D_zeta_es1Up);
      t->SetBranchAddress("D_zeta_es1Down", &D_zeta_es1Down);
      t->SetBranchAddress("D_zeta_es2Up", &D_zeta_es2Up);
      t->SetBranchAddress("D_zeta_es2Down", &D_zeta_es2Down);
      t->SetBranchAddress("D_zeta_responseUp", &D_zeta_ResponseUp);
      t->SetBranchAddress("D_zeta_responseDown", &D_zeta_ResponseDown);
      t->SetBranchAddress("D_zeta_resolutionUp", &D_zeta_ResolutionUp);
      t->SetBranchAddress("D_zeta_resolutionDown", &D_zeta_ResolutionDown);
      t->SetBranchAddress("D_zeta_UESUp", &D_zeta_UESUp);
      t->SetBranchAddress("D_zeta_UESDown", &D_zeta_UESDown);
      t->SetBranchAddress("D_zeta_JERUp", &D_zeta_JERUp);
      t->SetBranchAddress("D_zeta_JERDown", &D_zeta_JERDown);
      t->SetBranchAddress("D_zeta_JetAbsoluteUp", &D_zeta_JetAbsoluteUp);
      t->SetBranchAddress("D_zeta_JetAbsoluteDown", &D_zeta_JetAbsoluteDown);
      t->SetBranchAddress("D_zeta_JetAbsoluteyearUp", &D_zeta_JetAbsoluteyearUp);
      t->SetBranchAddress("D_zeta_JetAbsoluteyearDown", &D_zeta_JetAbsoluteyearDown);
      t->SetBranchAddress("D_zeta_JetBBEC1Up", &D_zeta_JetBBEC1Up);
      t->SetBranchAddress("D_zeta_JetBBEC1Down", &D_zeta_JetBBEC1Down);
      t->SetBranchAddress("D_zeta_JetBBEC1yearUp", &D_zeta_JetBBEC1yearUp);
      t->SetBranchAddress("D_zeta_JetBBEC1yearDown", &D_zeta_JetBBEC1yearDown);
      t->SetBranchAddress("D_zeta_JetEC2Up", &D_zeta_JetEC2Up);
      t->SetBranchAddress("D_zeta_JetEC2Down", &D_zeta_JetEC2Down);
      t->SetBranchAddress("D_zeta_JetEC2yearUp", &D_zeta_JetEC2yearUp);
      t->SetBranchAddress("D_zeta_JetEC2yearDown", &D_zeta_JetEC2yearDown);
      t->SetBranchAddress("D_zeta_JetFlavorQCDUp", &D_zeta_JetFlavorQCDUp);
      t->SetBranchAddress("D_zeta_JetFlavorQCDDown", &D_zeta_JetFlavorQCDDown);
      t->SetBranchAddress("D_zeta_JetHFUp", &D_zeta_JetHFUp);
      t->SetBranchAddress("D_zeta_JetHFDown", &D_zeta_JetHFDown);
      t->SetBranchAddress("D_zeta_JetHFyearUp", &D_zeta_JetHFyearUp);
      t->SetBranchAddress("D_zeta_JetHFyearDown", &D_zeta_JetHFyearDown);
      t->SetBranchAddress("D_zeta_JetRelativeBalUp", &D_zeta_JetRelativeBalUp);
      t->SetBranchAddress("D_zeta_JetRelativeBalDown", &D_zeta_JetRelativeBalDown);
      t->SetBranchAddress("D_zeta_JetRelativeSampleUp", &D_zeta_JetRelativeSampleUp);
      t->SetBranchAddress("D_zeta_JetRelativeSampleDown", &D_zeta_JetRelativeSampleDown);

      // New Change
      t->SetBranchAddress("bpt_deepflavour_2", &bpt_deepflavour_2);
      t->SetBranchAddress("bpt_deepflavour_JERUp_2", &bpt_deepflavour_JERUp_2);
      t->SetBranchAddress("bpt_deepflavour_JERDown_2", &bpt_deepflavour_JERDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteUp_2", &bpt_deepflavour_JetAbsoluteUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteDown_2", &bpt_deepflavour_JetAbsoluteDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteyearUp_2", &bpt_deepflavour_JetAbsoluteyearUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteyearDown_2", &bpt_deepflavour_JetAbsoluteyearDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1Up_2", &bpt_deepflavour_JetBBEC1Up_2);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1Down_2", &bpt_deepflavour_JetBBEC1Down_2);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1yearUp_2", &bpt_deepflavour_JetBBEC1yearUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1yearDown_2", &bpt_deepflavour_JetBBEC1yearDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetEC2Up_2", &bpt_deepflavour_JetEC2Up_2);
      t->SetBranchAddress("bpt_deepflavour_JetEC2Down_2", &bpt_deepflavour_JetEC2Down_2);
      t->SetBranchAddress("bpt_deepflavour_JetEC2yearUp_2", &bpt_deepflavour_JetEC2yearUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetEC2yearDown_2", &bpt_deepflavour_JetEC2yearDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetFlavorQCDUp_2", &bpt_deepflavour_JetFlavorQCDUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetFlavorQCDDown_2", &bpt_deepflavour_JetFlavorQCDDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetHFUp_2", &bpt_deepflavour_JetHFUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetHFDown_2", &bpt_deepflavour_JetHFDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetHFyearUp_2", &bpt_deepflavour_JetHFyearUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetHFyearDown_2", &bpt_deepflavour_JetHFyearDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeBalUp_2", &bpt_deepflavour_JetRelativeBalUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeBalDown_2", &bpt_deepflavour_JetRelativeBalDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeSampleUp_2", &bpt_deepflavour_JetRelativeSampleUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeSampleDown_2", &bpt_deepflavour_JetRelativeSampleDown_2);
      t->SetBranchAddress("beta_deepflavour_2", &beta_deepflavour_2);
      t->SetBranchAddress("bphi_deepflavour_2", &bphi_deepflavour_2);
      t->SetBranchAddress("bm_deepflavour_2", &bm_deepflavour_2);
      t->SetBranchAddress("bm_deepflavour_JERUp_2", &bm_deepflavour_JERUp_2);
      t->SetBranchAddress("bm_deepflavour_JERDown_2", &bm_deepflavour_JERDown_2);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteUp_2", &bm_deepflavour_JetAbsoluteUp_2);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteDown_2", &bm_deepflavour_JetAbsoluteDown_2);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteyearUp_2", &bm_deepflavour_JetAbsoluteyearUp_2);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteyearDown_2", &bm_deepflavour_JetAbsoluteyearDown_2);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1Up_2", &bm_deepflavour_JetBBEC1Up_2);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1Down_2", &bm_deepflavour_JetBBEC1Down_2);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1yearUp_2", &bm_deepflavour_JetBBEC1yearUp_2);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1yearDown_2", &bm_deepflavour_JetBBEC1yearDown_2);
      t->SetBranchAddress("bm_deepflavour_JetEC2Up_2", &bm_deepflavour_JetEC2Up_2);
      t->SetBranchAddress("bm_deepflavour_JetEC2Down_2", &bm_deepflavour_JetEC2Down_2);
      t->SetBranchAddress("bm_deepflavour_JetEC2yearUp_2", &bm_deepflavour_JetEC2yearUp_2);
      t->SetBranchAddress("bm_deepflavour_JetEC2yearDown_2", &bm_deepflavour_JetEC2yearDown_2);
      t->SetBranchAddress("bm_deepflavour_JetFlavorQCDUp_2", &bm_deepflavour_JetFlavorQCDUp_2);
      t->SetBranchAddress("bm_deepflavour_JetFlavorQCDDown_2", &bm_deepflavour_JetFlavorQCDDown_2);
      t->SetBranchAddress("bm_deepflavour_JetHFUp_2", &bm_deepflavour_JetHFUp_2);
      t->SetBranchAddress("bm_deepflavour_JetHFDown_2", &bm_deepflavour_JetHFDown_2);
      t->SetBranchAddress("bm_deepflavour_JetHFyearUp_2", &bm_deepflavour_JetHFyearUp_2);
      t->SetBranchAddress("bm_deepflavour_JetHFyearDown_2", &bm_deepflavour_JetHFyearDown_2);
      t->SetBranchAddress("bm_deepflavour_JetRelativeBalUp_2", &bm_deepflavour_JetRelativeBalUp_2);
      t->SetBranchAddress("bm_deepflavour_JetRelativeBalDown_2", &bm_deepflavour_JetRelativeBalDown_2);
      t->SetBranchAddress("bm_deepflavour_JetRelativeSampleUp_2", &bm_deepflavour_JetRelativeSampleUp_2);
      t->SetBranchAddress("bm_deepflavour_JetRelativeSampleDown_2", &bm_deepflavour_JetRelativeSampleDown_2);

      printf("Found tree -> weighting\n");
      for (Int_t i = 0; i < t->GetEntries(); ++i) { // t->GetEntries()
        t->GetEntry(i);

  int booster_idx = -1;
	if(bpt_deepflavour_1 > 0 && bpt_deepflavour_2 <= 0) {
    if (std::string(key->GetName()) == "mutau_tree") {
      booster_idx = 0;
    } else if (std::string(key->GetName()) == "etau_tree") {
      booster_idx = 1;
    } else if (std::string(key->GetName()) == "emu_tree") {
      booster_idx = 2;
    }
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_nominal, metphi_nominal, D_zeta_nominal, BDTval_nominal);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_es1Up, eta_1, phi_1, m_1_es1Up, pt_2_nominal, eta_2, phi_2, m_2, met_es1Up, metphi_es1Up, D_zeta_es1Up, BDTval_es1Up);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_es1Down, eta_1, phi_1, m_1_es1Down, pt_2_nominal, eta_2, phi_2, m_2, met_es1Down, metphi_es1Down, D_zeta_es1Down, BDTval_es1Down);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_es2Up, eta_2, phi_2, m_2_es2Up, met_es2Up, metphi_es2Up, D_zeta_es2Up, BDTval_es2Up);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_es2Down, eta_2, phi_2, m_2_es2Down, met_es2Down, metphi_es2Down, D_zeta_es2Down, BDTval_es2Down);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_UESUp, metphi_UESUp, D_zeta_UESUp, BDTval_UESUp);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_UESDown, metphi_UESDown, D_zeta_UESDown, BDTval_UESDown);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_ResponseUp, metphi_ResponseUp, D_zeta_ResponseUp, BDTval_ResponseUp);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_ResponseDown, metphi_ResponseDown, D_zeta_ResponseDown, BDTval_ResponseDown);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_ResolutionUp, metphi_ResolutionUp, D_zeta_ResolutionUp, BDTval_ResolutionUp);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_ResolutionDown, metphi_ResolutionDown, D_zeta_ResolutionDown, BDTval_ResolutionDown);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetAbsoluteUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetAbsoluteUp, metphi_JetAbsoluteUp, D_zeta_JetAbsoluteUp, BDTval_JetAbsoluteUp);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetAbsoluteDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetAbsoluteDown, metphi_JetAbsoluteDown, D_zeta_JetAbsoluteDown, BDTval_JetAbsoluteDown);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetAbsoluteyearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteyearUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetAbsoluteyearUp, metphi_JetAbsoluteyearUp, D_zeta_JetAbsoluteyearUp, BDTval_JetAbsoluteyearUp);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetAbsoluteyearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteyearDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetAbsoluteyearDown, metphi_JetAbsoluteyearDown, D_zeta_JetAbsoluteyearDown, BDTval_JetAbsoluteyearDown);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetBBEC1Up_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1Up_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetBBEC1Up, metphi_JetBBEC1Up, D_zeta_JetBBEC1Up, BDTval_JetBBEC1Up);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetBBEC1Down_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1Down_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetBBEC1Down, metphi_JetBBEC1Down, D_zeta_JetBBEC1Down, BDTval_JetBBEC1Down);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetBBEC1yearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1yearUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetBBEC1yearUp, metphi_JetBBEC1yearUp, D_zeta_JetBBEC1yearUp, BDTval_JetBBEC1yearUp);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetBBEC1yearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1yearDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetBBEC1yearDown, metphi_JetBBEC1yearDown, D_zeta_JetBBEC1yearDown, BDTval_JetBBEC1yearDown);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetEC2Up_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2Up_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetEC2Up, metphi_JetEC2Up, D_zeta_JetEC2Up, BDTval_JetEC2Up);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetEC2Down_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2Down_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetEC2Down, metphi_JetEC2Down, D_zeta_JetEC2Down, BDTval_JetEC2Down);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetEC2yearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2yearUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetEC2yearUp, metphi_JetEC2yearUp, D_zeta_JetEC2yearUp, BDTval_JetEC2yearUp);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetEC2yearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2yearDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetEC2yearDown, metphi_JetEC2yearDown, D_zeta_JetEC2yearDown, BDTval_JetEC2yearDown);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetFlavorQCDUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetFlavorQCDUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetFlavorQCDUp, metphi_JetFlavorQCDUp, D_zeta_JetFlavorQCDUp, BDTval_JetFlavorQCDUp);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetFlavorQCDDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetFlavorQCDDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetFlavorQCDDown, metphi_JetFlavorQCDDown, D_zeta_JetFlavorQCDDown, BDTval_JetFlavorQCDDown);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetHFUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetHFUp, metphi_JetHFUp, D_zeta_JetHFUp, BDTval_JetHFUp);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetHFDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetHFDown, metphi_JetHFDown, D_zeta_JetHFDown, BDTval_JetHFDown);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetHFyearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFyearUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetHFyearUp, metphi_JetHFyearUp, D_zeta_JetHFyearUp, BDTval_JetHFyearUp);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetHFyearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFyearDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetHFyearDown, metphi_JetHFyearDown, D_zeta_JetHFyearDown, BDTval_JetHFyearDown);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetRelativeBalUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeBalUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetRelativeBalUp, metphi_JetRelativeBalUp, D_zeta_JetRelativeBalUp, BDTval_JetRelativeBalUp);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetRelativeBalDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeBalDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetRelativeBalDown, metphi_JetRelativeBalDown, D_zeta_JetRelativeBalDown, BDTval_JetRelativeBalDown);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetRelativeSampleUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeSampleUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetRelativeSampleUp, metphi_JetRelativeSampleUp, D_zeta_JetRelativeSampleUp, BDTval_JetRelativeSampleUp);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetRelativeSampleDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeSampleDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetRelativeSampleDown, metphi_JetRelativeSampleDown, D_zeta_JetRelativeSampleDown, BDTval_JetRelativeSampleDown);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JERUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JERUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JERUp, metphi_JERUp, D_zeta_JERUp, BDTval_JERUp);
    runBDT_1bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JERDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JERDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JERDown, metphi_JERDown, D_zeta_JERDown, BDTval_JERDown);
	} else if (bpt_deepflavour_1 > 0 && bpt_deepflavour_2 > 0) {
    if (std::string(key->GetName()) == "mutau_tree") {
      booster_idx = 3;
    } else if (std::string(key->GetName()) == "etau_tree") {
      booster_idx = 4;
    } else if (std::string(key->GetName()) == "emu_tree") {
      booster_idx = 5;
    }
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_nominal, metphi_nominal, D_zeta_nominal, BDTval_nominal);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2, pt_1_es1Up, eta_1, phi_1, m_1_es1Up, pt_2_nominal, eta_2, phi_2, m_2, met_es1Up, metphi_es1Up, D_zeta_es1Up, BDTval_es1Up);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2, pt_1_es1Down, eta_1, phi_1, m_1_es1Down, pt_2_nominal, eta_2, phi_2, m_2, met_es1Down, metphi_es1Down, D_zeta_es1Down, BDTval_es1Down);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_es2Up, eta_2, phi_2, m_2_es2Up, met_es2Up, metphi_es2Up, D_zeta_es2Up, BDTval_es2Up);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_es2Down, eta_2, phi_2, m_2_es2Down, met_es2Down, metphi_es2Down, D_zeta_es2Down, BDTval_es2Down);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_UESUp, metphi_UESUp, D_zeta_UESUp, BDTval_UESUp);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_UESDown, metphi_UESDown, D_zeta_UESDown, BDTval_UESDown);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_ResponseUp, metphi_ResponseUp, D_zeta_ResponseUp, BDTval_ResponseUp);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_ResponseDown, metphi_ResponseDown, D_zeta_ResponseDown, BDTval_ResponseDown);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_ResolutionUp, metphi_ResolutionUp, D_zeta_ResolutionUp, BDTval_ResolutionUp);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_ResolutionDown, metphi_ResolutionDown, D_zeta_ResolutionDown, BDTval_ResolutionDown);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetAbsoluteUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteUp_1, bpt_deepflavour_JetAbsoluteUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetAbsoluteUp_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetAbsoluteUp, metphi_JetAbsoluteUp, D_zeta_JetAbsoluteUp, BDTval_JetAbsoluteUp);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetAbsoluteDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteDown_1, bpt_deepflavour_JetAbsoluteDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetAbsoluteDown_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetAbsoluteDown, metphi_JetAbsoluteDown, D_zeta_JetAbsoluteDown, BDTval_JetAbsoluteDown);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetAbsoluteyearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteyearUp_1, bpt_deepflavour_JetAbsoluteyearUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetAbsoluteyearUp_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetAbsoluteyearUp, metphi_JetAbsoluteyearUp, D_zeta_JetAbsoluteyearUp, BDTval_JetAbsoluteyearUp);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetAbsoluteyearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteyearDown_1, bpt_deepflavour_JetAbsoluteyearDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetAbsoluteyearDown_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetAbsoluteyearDown, metphi_JetAbsoluteyearDown, D_zeta_JetAbsoluteyearDown, BDTval_JetAbsoluteyearDown);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetBBEC1Up_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1Up_1, bpt_deepflavour_JetBBEC1Up_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetBBEC1Up_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetBBEC1Up, metphi_JetBBEC1Up, D_zeta_JetBBEC1Up, BDTval_JetBBEC1Up);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetBBEC1Down_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1Down_1, bpt_deepflavour_JetBBEC1Down_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetBBEC1Down_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetBBEC1Down, metphi_JetBBEC1Down, D_zeta_JetBBEC1Down, BDTval_JetBBEC1Down);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetBBEC1yearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1yearUp_1, bpt_deepflavour_JetBBEC1yearUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetBBEC1yearUp_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetBBEC1yearUp, metphi_JetBBEC1yearUp, D_zeta_JetBBEC1yearUp, BDTval_JetBBEC1yearUp);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetBBEC1yearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1yearDown_1, bpt_deepflavour_JetBBEC1yearDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetBBEC1yearDown_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetBBEC1yearDown, metphi_JetBBEC1yearDown, D_zeta_JetBBEC1yearDown, BDTval_JetBBEC1yearDown);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetEC2Up_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2Up_1, bpt_deepflavour_JetEC2Up_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetEC2Up_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetEC2Up, metphi_JetEC2Up, D_zeta_JetEC2Up, BDTval_JetEC2Up);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetEC2Down_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2Down_1, bpt_deepflavour_JetEC2Down_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetEC2Down_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetEC2Down, metphi_JetEC2Down, D_zeta_JetEC2Down, BDTval_JetEC2Down);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetEC2yearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2yearUp_1, bpt_deepflavour_JetEC2yearUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetEC2yearUp_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetEC2yearUp, metphi_JetEC2yearUp, D_zeta_JetEC2yearUp, BDTval_JetEC2yearUp);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetEC2yearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2yearDown_1, bpt_deepflavour_JetEC2yearDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetEC2yearDown_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetEC2yearDown, metphi_JetEC2yearDown, D_zeta_JetEC2yearDown, BDTval_JetEC2yearDown);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetFlavorQCDUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetFlavorQCDUp_1, bpt_deepflavour_JetFlavorQCDUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetFlavorQCDUp_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetFlavorQCDUp, metphi_JetFlavorQCDUp, D_zeta_JetFlavorQCDUp, BDTval_JetFlavorQCDUp);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetFlavorQCDDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetFlavorQCDDown_1, bpt_deepflavour_JetFlavorQCDDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetFlavorQCDDown_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetFlavorQCDDown, metphi_JetFlavorQCDDown, D_zeta_JetFlavorQCDDown, BDTval_JetFlavorQCDDown);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetHFUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFUp_1, bpt_deepflavour_JetHFUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetHFUp_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetHFUp, metphi_JetHFUp, D_zeta_JetHFUp, BDTval_JetHFUp);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetHFDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFDown_1, bpt_deepflavour_JetHFDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetHFDown_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetHFDown, metphi_JetHFDown, D_zeta_JetHFDown, BDTval_JetHFDown);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetHFyearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFyearUp_1, bpt_deepflavour_JetHFyearUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetHFyearUp_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetHFyearUp, metphi_JetHFyearUp, D_zeta_JetHFyearUp, BDTval_JetHFyearUp);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetHFyearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFyearDown_1, bpt_deepflavour_JetHFyearDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetHFyearDown_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetHFyearDown, metphi_JetHFyearDown, D_zeta_JetHFyearDown, BDTval_JetHFyearDown);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetRelativeBalUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeBalUp_1, bpt_deepflavour_JetRelativeBalUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetRelativeBalUp_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetRelativeBalUp, metphi_JetRelativeBalUp, D_zeta_JetRelativeBalUp, BDTval_JetRelativeBalUp);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetRelativeBalDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeBalDown_1, bpt_deepflavour_JetRelativeBalDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetRelativeBalDown_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetRelativeBalDown, metphi_JetRelativeBalDown, D_zeta_JetRelativeBalDown, BDTval_JetRelativeBalDown);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetRelativeSampleUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeSampleUp_1, bpt_deepflavour_JetRelativeSampleUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetRelativeSampleUp_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetRelativeSampleUp, metphi_JetRelativeSampleUp, D_zeta_JetRelativeSampleUp, BDTval_JetRelativeSampleUp);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JetRelativeSampleDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeSampleDown_1, bpt_deepflavour_JetRelativeSampleDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetRelativeSampleDown_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetRelativeSampleDown, metphi_JetRelativeSampleDown, D_zeta_JetRelativeSampleDown, BDTval_JetRelativeSampleDown);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JERUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JERUp_1, bpt_deepflavour_JERUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JERUp_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JERUp, metphi_JERUp, D_zeta_JERUp, BDTval_JERUp);
    runBDT_2bjet(booster_idx, boosters[booster_idx], bpt_deepflavour_JERDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JERDown_1, bpt_deepflavour_JERDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JERDown_2, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JERDown, metphi_JERDown, D_zeta_JERDown, BDTval_JERDown);
  } else {
    BDTval_nominal = -10;
    BDTval_es1Up = -10;
    BDTval_es1Down = -10;
    BDTval_es2Up = -10;
    BDTval_es2Down = -10;
    BDTval_UESUp = -10;
    BDTval_UESDown = -10;
    BDTval_ResponseUp = -10;
    BDTval_ResponseDown = -10;
    BDTval_ResolutionUp = -10;
    BDTval_ResolutionDown = -10;
    BDTval_JetAbsoluteUp = -10;
    BDTval_JetAbsoluteDown = -10;
    BDTval_JetAbsoluteyearUp = -10;
    BDTval_JetAbsoluteyearDown = -10;
    BDTval_JetBBEC1Up = -10;
    BDTval_JetBBEC1Down = -10;
    BDTval_JetBBEC1yearUp = -10;
    BDTval_JetBBEC1yearDown = -10;
    BDTval_JetEC2Up = -10;
    BDTval_JetEC2Down = -10;
    BDTval_JetEC2yearUp = -10;
    BDTval_JetEC2yearDown = -10;
    BDTval_JetFlavorQCDUp = -10;
    BDTval_JetFlavorQCDDown = -10;
    BDTval_JetHFUp = -10;
    BDTval_JetHFDown = -10;
    BDTval_JetHFyearUp = -10;
    BDTval_JetHFyearDown = -10;
    BDTval_JetRelativeBalUp = -10;
    BDTval_JetRelativeBalDown = -10;
    BDTval_JetRelativeSampleUp = -10;
    BDTval_JetRelativeSampleDown = -10;
    BDTval_JERUp = -10;
    BDTval_JERDown = -10;
  }


        newBranch1->Fill();
        newBranch1_es1U->Fill();
        newBranch1_es1D->Fill();
        newBranch1_es2U->Fill();
        newBranch1_es2D->Fill();
        newBranch1UU->Fill();
        newBranch1UD->Fill();
        newBranch1ResponseU->Fill();
        newBranch1ResponseD->Fill();
        newBranch1ResolutionU->Fill();
        newBranch1ResolutionD->Fill();
        newBranch1JetAbsoluteU->Fill();
        newBranch1JetAbsoluteD->Fill();
        newBranch1JetAbsoluteyearU->Fill();
        newBranch1JetAbsoluteyearD->Fill();
        newBranch1JetBBEC1U->Fill();
        newBranch1JetBBEC1D->Fill();
        newBranch1JetBBEC1yearU->Fill();
        newBranch1JetBBEC1yearD->Fill();
        newBranch1JetEC2U->Fill();
        newBranch1JetEC2D->Fill();
        newBranch1JetEC2yearU->Fill();
        newBranch1JetEC2yearD->Fill();
        newBranch1JetFlavorQCDU->Fill();
        newBranch1JetFlavorQCDD->Fill();
        newBranch1JetHFU->Fill();
        newBranch1JetHFD->Fill();
        newBranch1JetHFyearU->Fill();
        newBranch1JetHFyearD->Fill();
        newBranch1JetRelativeBalU->Fill();
        newBranch1JetRelativeBalD->Fill();
        newBranch1JetRelativeSampleU->Fill();
        newBranch1JetRelativeSampleD->Fill();
        newBranch1JERU->Fill();
        newBranch1JERD->Fill();
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

void runBDT_1bjet(int booster_idx, pat::XGBooster booster, float jet_pt, float jet_eta, float jet_phi, float jet_mass, float pt_1, float eta_1, float phi_1, float mass_1, float pt_2, float eta_2, float phi_2, float mass_2, float met, float met_phi, float d_zeta, float& BDTval) {
  ROOT::Math::PtEtaPhiMVector bjet1(jet_pt, jet_eta, jet_phi, jet_mass);
  ROOT::Math::PtEtaPhiMVector tau1(pt_1, eta_1, phi_1, mass_1);
  ROOT::Math::PtEtaPhiMVector tau2(pt_2, eta_2, phi_2, mass_2);
  ROOT::Math::PtEtaPhiMVector tau12 = tau1 + tau2;
  ROOT::Math::PtEtaPhiMVector all = tau12 + bjet1;

  if (booster_idx == 0) {
    float in_mT_thMET, in_mT_mMET, in_d_zeta, in_mutau_pt, in_m_b1mt, in_bmt_dR, in_b1th_dR;
    in_mT_thMET = compute_mt(pt_2, phi_2, met, met_phi); // DEFINE ORDER
    in_mT_mMET = compute_mt(pt_1, phi_1, met, met_phi); // DEFINE ORDER
    in_d_zeta = d_zeta;
    in_mutau_pt = tau12.Pt();
    in_m_b1mt = all.M();
    in_bmt_dR = compute_deltaR(jet_eta, tau12.Eta(), jet_phi, tau12.Phi());
    in_b1th_dR = compute_deltaR(jet_eta, eta_2, jet_phi, phi_2); // DEFINE ORDER

    booster.set("mT_thMET", in_mT_thMET);
    booster.set("mT_mMET", in_mT_mMET);
    booster.set("d_zeta", in_d_zeta);
    booster.set("mutau_pt", in_mutau_pt);
    booster.set("m_b1mt", in_m_b1mt);
    booster.set("bmt_dR", in_bmt_dR);
    booster.set("b1th_dR", in_b1th_dR);
    booster.set("njets", 1.0); // DEFINE
    BDTval = booster.predict();
  } else if (booster_idx == 1) {
    float in_pt_vis_nominal, in_D_zeta_nominal, in_b1e_dR, in_m_btautau_vis_nominal, in_mtMET_1_nominal, in_b1th_dR, in_mtMET_2_nominal;
    in_pt_vis_nominal = tau12.Pt(); // DEFINE CONSTITUENTS
    in_D_zeta_nominal = d_zeta;
    in_b1e_dR = compute_deltaR(jet_eta, eta_1, jet_phi, phi_1); // DEFINE ORDER
    in_m_btautau_vis_nominal = all.M();
    in_mtMET_1_nominal = compute_mt(pt_1, phi_1, met, met_phi);
    in_b1th_dR = compute_deltaR(jet_eta, eta_2, jet_phi, phi_2); // DEFINE ORDER
    in_mtMET_2_nominal = compute_mt(pt_2, phi_2, met, met_phi);

    booster.set("pt_vis_nominal", in_pt_vis_nominal);
    booster.set("D_zeta_nominal", in_D_zeta_nominal);
    booster.set("b1e_dR", in_b1e_dR);
    booster.set("m_btautau_vis_nominal", in_m_btautau_vis_nominal);
    booster.set("mtMET_1_nominal", in_mtMET_1_nominal);
    booster.set("b1th_dR", in_b1th_dR);
    booster.set("mtMET_2_nominal", in_mtMET_2_nominal);
    booster.set("njets", 1.0);
    BDTval = booster.predict();
  } else if (booster_idx == 2) {
    float in_pt_vis_nominal, in_D_zeta_nominal, in_mT_b1MET, in_pt_1_nominal, in_m_btautau_vis_nominal, in_mtMET_1_nominal, in_b1emu_dR;
    in_pt_vis_nominal = tau12.Pt(); // DEFINE CONSTITUENTS
    in_D_zeta_nominal = d_zeta;
    in_mT_b1MET = compute_mt(jet_pt, jet_phi, met, met_phi);
    in_pt_1_nominal = pt_1;
    in_m_btautau_vis_nominal = all.M();
    in_mtMET_1_nominal = compute_mt(pt_1, phi_1, met, met_phi);
    in_b1emu_dR = compute_deltaR(jet_eta, tau12.Eta(), jet_phi, tau12.Phi());

    booster.set("pt_vis_nominal", in_pt_vis_nominal);
    booster.set("D_zeta_nominal", in_D_zeta_nominal);
    booster.set("mT_b1MET", in_mT_b1MET);
    booster.set("pt_1_nominal", in_pt_1_nominal);
    booster.set("m_btautau_vis_nominal", in_m_btautau_vis_nominal);
    booster.set("mtMET_1_nominal", in_mtMET_1_nominal);
    booster.set("b1emu_dR", in_b1emu_dR);
    booster.set("njets", 1.0);
    BDTval = booster.predict();
  } else {
    throw std::invalid_argument( "Invalid booster_idx" );
  }
}

void runBDT_2bjet(int booster_idx, pat::XGBooster booster, float jet_1_pt, float jet_1_eta, float jet_1_phi, float jet_1_mass, float jet_2_pt, float jet_2_eta, float jet_2_phi, float jet_2_mass, float pt_1, float eta_1, float phi_1, float mass_1, float pt_2, float eta_2, float phi_2, float mass_2, float met, float met_phi, float d_zeta, float& BDTval) {
  ROOT::Math::PtEtaPhiMVector bjet1(jet_1_pt, jet_1_eta, jet_1_phi, jet_1_mass);
  ROOT::Math::PtEtaPhiMVector bjet2(jet_2_pt, jet_2_eta, jet_2_phi, jet_2_mass);
  ROOT::Math::PtEtaPhiMVector bjet12 = bjet1 + bjet2;
  ROOT::Math::PtEtaPhiMVector tau1(pt_1, eta_1, phi_1, mass_1);
  ROOT::Math::PtEtaPhiMVector tau2(pt_2, eta_2, phi_2, mass_2);
  ROOT::Math::PtEtaPhiMVector tau12 = tau1 + tau2;

  if (booster_idx == 3) {
    float in_m_b2mt, in_bmt_dR, in_b2th_dR, in_mT_mMET, in_m_bbmt, in_d_ma;
    in_m_b2mt = (bjet2 + tau12).M();
    in_bmt_dR = compute_deltaR(jet_1_eta, tau12.Eta(), jet_1_phi, tau12.Phi()); // DEFINE CONSTITUTENTS
    in_b2th_dR = compute_deltaR(jet_2_eta, eta_2, jet_2_phi, phi_2); // DEFINE ORDER
    in_mT_mMET = compute_mt(pt_1, phi_1, met, met_phi); // DEFINE ORDER
    in_m_bbmt = (bjet12 + tau12).M();
    in_d_ma = (bjet12.M() / tau12.M()) / tau12.M();

    booster.set("m_b2mt", in_m_b2mt);
    booster.set("bmt_dR", in_bmt_dR);
    booster.set("b2th_dR", in_b2th_dR);
    booster.set("mT_mMET", in_mT_mMET);
    booster.set("m_bbmt", in_m_bbmt);
    booster.set("d_ma", in_d_ma);
    BDTval = booster.predict();
  } else if (booster_idx == 4) {
    float in_mbb, in_m_btautau_vis_nominal, in_mtMET_1_nominal, in_b1th_dR, in_b1e_dR, in_b2th_dR;
    in_mbb = bjet12.M();
    in_m_btautau_vis_nominal = (bjet1 + tau12).M(); // DEFINE CONSTITUTENTS
    in_mtMET_1_nominal = compute_mt(pt_1, phi_1, met, met_phi);
    in_b1th_dR = compute_deltaR(jet_1_eta, eta_2, jet_1_phi, phi_2); // DEFINE ORDER
    in_b1e_dR = compute_deltaR(jet_1_eta, eta_1, jet_1_phi, phi_1); // DEFINE ORDER
    in_b2th_dR = compute_deltaR(jet_2_eta, eta_2, jet_2_phi, phi_2); // DEFINE ORDER

    booster.set("mbb", in_mbb);
    booster.set("m_btautau_vis_nominal", in_m_btautau_vis_nominal);
    booster.set("mtMET_1_nominal", in_mtMET_1_nominal);
    booster.set("b1th_dR", in_b1th_dR);
    booster.set("b1e_dR", in_b1e_dR);
    booster.set("b2th_dR", in_b2th_dR);
    BDTval = booster.predict();
  } else if (booster_idx == 5) {
    float in_pt_vis_nominal, in_mT_b1MET, in_mtMET_1_nominal, in_b2emu_dR, in_d_ma;
    in_pt_vis_nominal = tau12.Pt(); // DEFINE CONSTITUTENTS
    in_mT_b1MET = compute_mt(jet_1_pt, jet_1_phi, met, met_phi);
    in_mtMET_1_nominal = compute_mt(pt_1, phi_1, met, met_phi);
    in_b2emu_dR = compute_deltaR(jet_2_eta, tau12.Eta(), jet_2_phi, tau12.Phi());
    in_d_ma = (bjet12.M() / tau12.M()) / tau12.M();

    booster.set("pt_vis_nominal", in_pt_vis_nominal);
    booster.set("mT_b1MET", in_mT_b1MET);
    booster.set("mtMET_1_nominal", in_mtMET_1_nominal);
    booster.set("b2emu_dR", in_b2emu_dR);
    booster.set("d_ma", in_d_ma);
    BDTval = booster.predict();
  } else {
    throw std::invalid_argument( "Invalid booster_idx" );
  }
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

  // loop on all entries of this directory
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
  // copy all objects and subdirs of file fname as a subdir of the current directory
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
  // prepare files to be copied
  if(gSystem->AccessPathName(parser.stringValue("inputFile").c_str())) {
    gSystem->CopyFile("hsimple.root", parser.stringValue("inputFile").c_str());
  }

  fNew->cd();
  CopyFile(parser.stringValue("inputFile").c_str(),parser);
  fNew->ls();
  fNew->Close();
}

