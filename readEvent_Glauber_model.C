#include <iostream>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>

void readEvent()
{
	double b1, Ecc1;
	int npart1, ncoll1;

    //importing root file and extracting tree
    TFile *file1 = new TFile("/Users/vedangtamhane/Projects/19Summer_VECC/RootFiles/PartMCGMOut_AuAu_100000e_0.2TeV.root");
    TTree *tr1 = (TTree*)file1->Get("tree");

    tr1->SetBranchAddress("b",&b1);
    tr1->SetBranchAddress("npart", &npart1);
    tr1->SetBranchAddress("ncoll", &ncoll1);

    int NEvent = (int)tr1->GetEntries();

    double b2, Ecc2;
	int npart2, ncoll2;

    //importing root file and extracting tree
    TFile *file2 = new TFile("/Users/vedangtamhane/Projects/19Summer_VECC/RootFiles/NuclMCGMOut_AuAu_100000e_0.2TeV.root");
    TTree *tr2 = (TTree*)file2->Get("tree");

    tr2->SetBranchAddress("b",&b2);
    tr2->SetBranchAddress("npart", &npart2);
    tr2->SetBranchAddress("ncoll", &ncoll2);

    int nbins = 100;
    double xlow = 0., xup = 18., ylow = 0., yup = 1600.;

    TProfile *histNColl1 = new TProfile("histNColl1", "",nbins, xlow, xup, ylow, yup);
	TProfile *histNPart1 = new TProfile("histNPart1", "",nbins, xlow, xup, ylow, yup);

    TProfile *histNColl2 = new TProfile("histNColl2", "",nbins, xlow, xup, ylow, yup);
	TProfile *histNPart2 = new TProfile("histNPart2", "",nbins, xlow, xup, ylow, yup);

	TProfile *histratio = new TProfile("histratio","", 50, 0., 400., 0.0, 5.); 

    for (int i=0; i<NEvent; i++)
    {
    	tr1->GetEntry(i);

        double rat1 = 2.*ncoll1/npart1;
    	histNPart1->Fill(b1, npart1,1);
    	histNColl1->Fill(b1, ncoll1,1);

    	tr2->GetEntry(i);

        double rat2 = 2.*ncoll2/npart2;
    	histNPart2->Fill(b2, npart2,1);
    	histNColl2->Fill(b2, ncoll2,1);

    	/*double rat = 1.*ncoll1/pow(npart1, 4./3.);
    	double t = npart1;
    	histratio->Fill(t, rat);*/
    } 

    //gStyle->SetPalette(51);

    TCanvas *c1 = new TCanvas("c1","MC_Glauber_model", 800, 800);

    histNPart1->GetXaxis()->SetTitle("Impact parametre(b)");
    histNPart1->GetXaxis()->SetTitleSize(0.05);
    histNPart1->GetXaxis()->SetTitleOffset(0.83);
    histNPart1->GetYaxis()->SetTitle("N_{Part}");
    histNPart1->GetYaxis()->SetLabelSize(0.03);
    histNPart1->GetYaxis()->SetTitleSize(0.05);
    histNPart1->GetYaxis()->SetTitleOffset(0.92);
    histNPart1->SetLineColor(kRed);
    histNPart1->Draw();

    histNColl1->GetXaxis()->SetTitle("Impact parametre (fm)");
    histNColl1->GetYaxis()->SetTitle("N_{Part} & N_{Coll}");
    histNColl1->GetYaxis()->SetLabelSize(0.03);
    histNColl1->GetYaxis()->SetTitleSize(0.05);
    histNColl1->GetYaxis()->SetTitleOffset(0.92);
    histNColl1->SetLineColor(kBlue);
    histNColl1->Draw("same");

    histNColl2->GetXaxis()->SetTitle("Impact parametre (fm)");
    histNColl2->GetYaxis()->SetTitle("N_{Coll} and N_{Part}");
    histNColl2->SetLineColor(kBlack);
    histNColl2->Draw("same");

    histNPart2->GetXaxis()->SetTitle("Impact parametre(b)");
    histNPart2->GetYaxis()->SetTitle("N_{Part}");
    histNPart2->SetLineColor(kBlue);
    histNPart2->Draw("same");

    TLegend *legend = new TLegend(100, 120);
    //legend->SetHeader("Au-Au Collisions at #sqrt{S_{NN}} = 0.2 TeV","C");
    legend->AddEntry("histNPart1", "N_{N-Part} Nucleons", "l");
    legend->AddEntry("histNColl1", "N_{Coll} at 200 GeV", "l");
    legend->AddEntry("histNPart2", "N_{Q-Part} Quarks", "l");
    legend->AddEntry("histNColl2", "N_{Coll} at 7.7 GeV", "l");
    legend->SetTextSize(0.05);
    legend->Draw();

}











