#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>

void readPartMCGM()
{
	int ncoll, npart;
	double b;

	TFile *file1 = new TFile("PartMCGMOut_AuAu_100000e_0.2TeV.root");
    TTree *tr1 = (TTree*)file1->Get("tree");

    tr1->SetBranchAddress("b",&b);
    tr1->SetBranchAddress("npart", &npart);
    tr1->SetBranchAddress("ncoll", &ncoll);

    int NEvent = (int)tr1->GetEntries();

    int nbins = 100;
    double xlow = 0., xup = 18., ylow = 0., yup = 1600.;

    TProfile *histNColl = new TProfile("histNColl", "",nbins, xlow, xup, ylow, yup);
	TProfile *histNPart = new TProfile("histNPart", "",nbins, xlow, xup, ylow, yup);

	for (int i=0; i<NEvent; i++)
    {
    	tr1->GetEntry(i);

    	histNPart->Fill(b, npart,1);
        //double rat = 2.*ncoll/npart;
    	histNColl->Fill(b, ncoll,1);
    }

    TCanvas *c1 = new TCanvas("c1","MC_Glauber_model", 800, 800);

    histNColl->GetXaxis()->SetTitle("Impact parametre b (fm)");
    histNColl->GetXaxis()->SetTitleSize(0.05);
    histNColl->GetXaxis()->SetTitleOffset(0.83);
    histNColl->GetYaxis()->SetTitle("N_{Part} and N_{Coll}");
    histNColl->GetYaxis()->SetLabelSize(0.03);
    histNColl->GetYaxis()->SetTitleSize(0.05);
    histNColl->GetYaxis()->SetTitleOffset(0.92);
    histNColl->SetLineColor(kBlue);
    histNColl->Draw();

    histNPart->SetLineColor(kRed);
    histNPart->Draw("same");

    TLegend *legend = new TLegend(100, 120);
    legend->SetHeader("Au-Au #sqrt{S_{NN}} = 0.2 TeV","C");
    legend->AddEntry("histNPart", "N_{Part}", "l");
    legend->AddEntry("histNColl", "N_{Coll}", "l");
    legend->SetTextSize(0.05);
    legend->Draw();
}




























