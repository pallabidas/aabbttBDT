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
void runBDTeval(pat::XGBooster booster, float jet_pt, float jet_eta, float jet_phi, float jet_mass, float pt_1, float eta_1, float phi_1, float mass_1, float pt_2, float eta_2, float phi_2, float mass_2, float met, float met_phi, float d_zeta, float& BDTval);
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
      std::cout<<"here!!!!"<<std::endl;
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

      printf("Found tree -> weighting\n");
      for (Int_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
	if(bpt_deepflavour_1 > 0) {
	  runBDTeval(booster, bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_nominal, metphi_nominal, D_zeta_nominal, BDTval_nominal);
	  runBDTeval(booster, bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_es1Up, eta_1, phi_1, m_1_es1Up, pt_2_nominal, eta_2, phi_2, m_2, met_es1Up, metphi_es1Up, D_zeta_es1Up, BDTval_es1Up);
	  runBDTeval(booster, bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_es1Down, eta_1, phi_1, m_1_es1Down, pt_2_nominal, eta_2, phi_2, m_2, met_es1Down, metphi_es1Down, D_zeta_es1Down, BDTval_es1Down);
	  runBDTeval(booster, bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_es2Up, eta_2, phi_2, m_2_es2Up, met_es2Up, metphi_es2Up, D_zeta_es2Up, BDTval_es2Up);
	  runBDTeval(booster, bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_es2Down, eta_2, phi_2, m_2_es2Down, met_es2Down, metphi_es2Down, D_zeta_es2Down, BDTval_es2Down);
          runBDTeval(booster, bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_UESUp, metphi_UESUp, D_zeta_UESUp, BDTval_UESUp);
          runBDTeval(booster, bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_UESDown, metphi_UESDown, D_zeta_UESDown, BDTval_UESDown);
          runBDTeval(booster, bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_ResponseUp, metphi_ResponseUp, D_zeta_ResponseUp, BDTval_ResponseUp);
          runBDTeval(booster, bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_ResponseDown, metphi_ResponseDown, D_zeta_ResponseDown, BDTval_ResponseDown);
          runBDTeval(booster, bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_ResolutionUp, metphi_ResolutionUp, D_zeta_ResolutionUp, BDTval_ResolutionUp);
          runBDTeval(booster, bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_ResolutionDown, metphi_ResolutionDown, D_zeta_ResolutionDown, BDTval_ResolutionDown);
          runBDTeval(booster, bpt_deepflavour_JetAbsoluteUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetAbsoluteUp, metphi_JetAbsoluteUp, D_zeta_JetAbsoluteUp, BDTval_JetAbsoluteUp);
          runBDTeval(booster, bpt_deepflavour_JetAbsoluteDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetAbsoluteDown, metphi_JetAbsoluteDown, D_zeta_JetAbsoluteDown, BDTval_JetAbsoluteDown);
          runBDTeval(booster, bpt_deepflavour_JetAbsoluteyearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteyearUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetAbsoluteyearUp, metphi_JetAbsoluteyearUp, D_zeta_JetAbsoluteyearUp, BDTval_JetAbsoluteyearUp);
          runBDTeval(booster, bpt_deepflavour_JetAbsoluteyearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteyearDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetAbsoluteyearDown, metphi_JetAbsoluteyearDown, D_zeta_JetAbsoluteyearDown, BDTval_JetAbsoluteyearDown);
          runBDTeval(booster, bpt_deepflavour_JetBBEC1Up_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1Up_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetBBEC1Up, metphi_JetBBEC1Up, D_zeta_JetBBEC1Up, BDTval_JetBBEC1Up);
          runBDTeval(booster, bpt_deepflavour_JetBBEC1Down_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1Down_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetBBEC1Down, metphi_JetBBEC1Down, D_zeta_JetBBEC1Down, BDTval_JetBBEC1Down);
          runBDTeval(booster, bpt_deepflavour_JetBBEC1yearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1yearUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetBBEC1yearUp, metphi_JetBBEC1yearUp, D_zeta_JetBBEC1yearUp, BDTval_JetBBEC1yearUp);
          runBDTeval(booster, bpt_deepflavour_JetBBEC1yearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1yearDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetBBEC1yearDown, metphi_JetBBEC1yearDown, D_zeta_JetBBEC1yearDown, BDTval_JetBBEC1yearDown);
          runBDTeval(booster, bpt_deepflavour_JetEC2Up_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2Up_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetEC2Up, metphi_JetEC2Up, D_zeta_JetEC2Up, BDTval_JetEC2Up);
          runBDTeval(booster, bpt_deepflavour_JetEC2Down_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2Down_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetEC2Down, metphi_JetEC2Down, D_zeta_JetEC2Down, BDTval_JetEC2Down);
          runBDTeval(booster, bpt_deepflavour_JetEC2yearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2yearUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetEC2yearUp, metphi_JetEC2yearUp, D_zeta_JetEC2yearUp, BDTval_JetEC2yearUp);
          runBDTeval(booster, bpt_deepflavour_JetEC2yearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2yearDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetEC2yearDown, metphi_JetEC2yearDown, D_zeta_JetEC2yearDown, BDTval_JetEC2yearDown);
          runBDTeval(booster, bpt_deepflavour_JetFlavorQCDUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetFlavorQCDUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetFlavorQCDUp, metphi_JetFlavorQCDUp, D_zeta_JetFlavorQCDUp, BDTval_JetFlavorQCDUp);
          runBDTeval(booster, bpt_deepflavour_JetFlavorQCDDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetFlavorQCDDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetFlavorQCDDown, metphi_JetFlavorQCDDown, D_zeta_JetFlavorQCDDown, BDTval_JetFlavorQCDDown);
          runBDTeval(booster, bpt_deepflavour_JetHFUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetHFUp, metphi_JetHFUp, D_zeta_JetHFUp, BDTval_JetHFUp);
          runBDTeval(booster, bpt_deepflavour_JetHFDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetHFDown, metphi_JetHFDown, D_zeta_JetHFDown, BDTval_JetHFDown);
          runBDTeval(booster, bpt_deepflavour_JetHFyearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFyearUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetHFyearUp, metphi_JetHFyearUp, D_zeta_JetHFyearUp, BDTval_JetHFyearUp);
          runBDTeval(booster, bpt_deepflavour_JetHFyearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFyearDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetHFyearDown, metphi_JetHFyearDown, D_zeta_JetHFyearDown, BDTval_JetHFyearDown);
          runBDTeval(booster, bpt_deepflavour_JetRelativeBalUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeBalUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetRelativeBalUp, metphi_JetRelativeBalUp, D_zeta_JetRelativeBalUp, BDTval_JetRelativeBalUp);
          runBDTeval(booster, bpt_deepflavour_JetRelativeBalDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeBalDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetRelativeBalDown, metphi_JetRelativeBalDown, D_zeta_JetRelativeBalDown, BDTval_JetRelativeBalDown);
          runBDTeval(booster, bpt_deepflavour_JetRelativeSampleUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeSampleUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetRelativeSampleUp, metphi_JetRelativeSampleUp, D_zeta_JetRelativeSampleUp, BDTval_JetRelativeSampleUp);
          runBDTeval(booster, bpt_deepflavour_JetRelativeSampleDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeSampleDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JetRelativeSampleDown, metphi_JetRelativeSampleDown, D_zeta_JetRelativeSampleDown, BDTval_JetRelativeSampleDown);
          runBDTeval(booster, bpt_deepflavour_JERUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JERUp_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JERUp, metphi_JERUp, D_zeta_JERUp, BDTval_JERUp);
          runBDTeval(booster, bpt_deepflavour_JERDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JERDown_1, pt_1_nominal, eta_1, phi_1, m_1, pt_2_nominal, eta_2, phi_2, m_2, met_JERDown, metphi_JERDown, D_zeta_JERDown, BDTval_JERDown);
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

void runBDTeval(pat::XGBooster booster, float jet_pt, float jet_eta, float jet_phi, float jet_mass, float pt_1, float eta_1, float phi_1, float mass_1, float pt_2, float eta_2, float phi_2, float mass_2, float met, float met_phi, float d_zeta, float& BDTval) {
  float in_Jet_pt, in_mT_b1MET, in_d_zeta, in_mutau_pt, in_bmt_dR;
  in_Jet_pt = jet_pt;
  in_mT_b1MET = compute_mt(jet_pt, jet_phi, met, met_phi);
  in_d_zeta = d_zeta;
  ROOT::Math::PtEtaPhiMVector bjet1(jet_pt, jet_eta, jet_phi, jet_mass);
  ROOT::Math::PtEtaPhiMVector tau1(pt_1, eta_1, phi_1, mass_1);
  ROOT::Math::PtEtaPhiMVector tau2(pt_2, eta_2, phi_2, mass_2);
  ROOT::Math::PtEtaPhiMVector tau12 = tau1 + tau2;
  in_mutau_pt = tau12.Pt();
  in_bmt_dR = compute_deltaR(jet_eta, tau12.Eta(), jet_phi, tau12.Phi());
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

