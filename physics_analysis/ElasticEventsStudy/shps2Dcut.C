void shps2Dcut(const char* ana_rootfile_name="ana_rootfile.root")
{

  TString path_to_rootfile{"/volatile/halla/sbs/adr/Rootfiles/gmn_parsed/SBS4/pass0"};
  //Empty histogram to make our cut
  //TH2D *h = new TH2D("h","",100,-10,10,100,-10,10);
  TFile *file = new TFile(Form("%s/%s",path_to_rootfile,ana_rootfile_name),"READ");
  TTree* T = new TTree("");
  TCanvas *c = new TCanvas("c","",800,600);
  h->Draw();
  
  //This line waits for the user to draw a polygon and then saves it as a graphical cut (TCutG)
  //Double click to finish cutting by hand
  TCutG *cutg = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG")); // making cut, store to CUTG
  c->Update();

  //Change the line color so it is obvious that the cut is done being drawn
  cutg->SetName("cut name");
  cutg->SetLineColor(kMagenta);
  cutg->SetLineWidth(2);
  cutg->Draw("PL");
  c->Update();

  //Write the cut to a files
  //TFile *file = new TFile("test.root","recreate");

  cutg->Write();

}