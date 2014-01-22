// myCutG.C
// S. Stave
// July, 2009
//
// --------------------------------------------------
{

  TCanvas *c1 = new TCanvas("c1","c1",50,50,800,600);
  c1->SetFillColor(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  TCutG *cut;
  TFile *f2 = new TFile("myCutG.root");
  
  cut = (TCutG*)f2->Get("cut");  
  cut->SetName("cut");
  cut->SetLineColor(2);
  cut->SetLineWidth(2);

  TChain *nt = new TChain("tree");
  
  nt->Add("myNtuple.root");

  c1->Divide(2,1);
  
  c1->cd(1);
  nt->Draw("y:x");
  nt->SetMarkerColor(2);
  nt->Draw("y:x","cut","same");
  Float_t nt_x, nt_y, nt_z;
  
  nt->SetBranchAddress("x",&nt_x);
  nt->SetBranchAddress("y",&nt_y);
  nt->SetBranchAddress("z",&nt_z);

  TH2F *h2 = new TH2F("h2","2D histo",100,-10,10,100,-10,10);
  TH2F *h2_cut = new TH2F("h2_cut","2D histo with cuts",100,-10,10,100,-10,10);

  Int_t nentries,nbytes;
  nentries=(Int_t)nt->GetEntries();
  Int_t i;
  
  for (i = 0;i<nentries;i++) {
    nbytes = nt->GetEntry(i);

    //Uncut histogram
    h2->Fill(nt_x,nt_y);

    //Cut histogram
    if (cut->IsInside(nt_x,nt_y)){
      h2_cut->Fill(nt_x,nt_y);
    }
    
  }
  h2_cut->SetMarkerColor(2);
  c1->cd(2);
  h2_cut->Draw();
  h2->Draw("same");

  c1->cd();

}
