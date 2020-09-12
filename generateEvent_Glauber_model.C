/*
	This is the Monté Carlo Glauber model developed by Vedang Tamhane, BSMS 2017 Indian Institute of 
	Science Education and Research, Tirupati, during summer project 2019 at Variable Energy 
	Cyclotron Centre, Kolkata, under guidance of Dr Premomoy Ghosh.
*/

#include <iostream>
#include <TRandom.h>
#include <TF1.h>
#include <TDatime.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

//********************GIVEN CONSTANTS********************************************************************************************************************


//Units of R, r, etc. are in fm. Units of area are in fm^2.

//Atomic mass number of the nuclei in collision events.
const int A = 197; //Nucleus A (Gold)
const Int_t NEvent = 100000;
const double sigma = 4.2; //For 200 GeV
//const double sigma = 3.08; //for 7.7 GeV
const double a  = 0.55;
const double d0 = TMath::Sqrt(sigma/TMath::Pi());
/*
 Beam Energy in GeV and corresponding values of sigma.
 E (GeV)    sigma(fm^2)
 7.7        3.08
 11.5       2.973
 17.2       2.95
 19.6       3.008
 27         3.194
 39         3.098
 62.4       3.155
 200        4.2
 2760       6.4
 5500       7.2
*/
//******************CALCULATED CONSTANTS**************************************************************************************

const double RA = (1.12 * pow(1. * A, 1./3.));

const double rho0A = 1. * A / (4./3. * TMath::Pi() * pow(RA ,3.));

TDatime* time1 = new TDatime;
TRandom* rand1 = new TRandom(time1->GetTime());

//*******************FUNCTIONS AND CODES***************************************************************************************

/*double hardsphere(double r)
{
    if (r<=RA) {return 1;}
    else if (r>=RA) {return 0;}
}*/

TF1 *distfun = new TF1("Wood-Saxon", "x*x/(1. + exp( (x*x - [0]*[0])/([1]*[1]) ))",0., 20.);
//TF1 *distfun = new TF1("Wood-Saxon", "[2]/(1. + exp( 1. + (x-[0])/([1]) ))",0., 20.);
//TF1 *distfun = new TF1("Wood-Saxon", "hardsphere(x)",0., 20.);
//TF1 *distfun = new TF1("distfun", "([0] + [1]*x + [2]*x*x)/(1. + exp( (x - [3])/[4] ))",0., 10.);

double *event(double b)
{
    double *rawdata = new double[3];

	double xA[A], yA[A], zA[A];
    double xB[A], yB[A], zB[A];

    //generating random coordinates for a nucleon A wrt Wood-Saxon Distribution.

    //sum of coordinates to calculate COM
    double xsumA = 0, xsumB = 0, ysumA = 0, ysumB = 0, zsumA = 0, zsumB = 0;

    distfun->SetParameter(0, RA);
    distfun->SetParameter(1, a);
    distfun->SetParameter(2, rho0A);

    int i = 0;
    while (i < A)
    {
        double r = distfun->GetRandom();
        double ctheta = 2*rand1->Rndm() - 1;
        double stheta = TMath::Sqrt(1 - ctheta*ctheta);
        double phi = 2.*TMath::Pi() * rand1->Rndm();

        double xcoord = r * stheta* cos(phi);
        double ycoord = r * stheta* sin(phi);
        double zcoord = r * ctheta;

        double rho = rho0A/( 1 + TMath::Exp(1. + (r-RA)/a));
        double rho0 = rho0A/(1 + TMath::Exp(1. + (-RA)/a));
        double p1 = gRandom->Uniform(0.,1.);

        //if ((rho/rho0) < p1)
        {
            xA[i] = xcoord;
            xsumA = xsumA + xcoord;

            yA[i] = ycoord;
            ysumA = ysumA + ycoord;

            zA[i] = zcoord;
            zsumA = zsumA + zcoord;

            i++;

        }
    }

    int j = 0;
    while (j < A)
    {
        double r = distfun->GetRandom();
        double ctheta = 2*rand1->Rndm() - 1;
        double stheta = TMath::Sqrt(1 - ctheta*ctheta);
        double phi = 2.*TMath::Pi() * rand1->Rndm();

        double xcoord = r * stheta* cos(phi);
        double ycoord = r * stheta* sin(phi);
        double zcoord = r * ctheta;

        double rho = rho0A/( 1 + TMath::Exp(1. + (r-RA)/a));
        double rho0 = rho0A/(1 + TMath::Exp(1. + (-RA)/a));
        double p2 = gRandom->Uniform(0.,1.);

        //if ((rho/rho0) < p2)
        {
            xB[j] = xcoord;
            xsumB = xsumB + xcoord;

            yB[j] = ycoord;
            ysumB = ysumB + ycoord;

            zB[j] = zcoord;
            zsumB = zsumB + zcoord;

            j++;

        }
    }
//******************************************************************************************************************************************************
    
    /*for (int i = 0; i<A; i++)
    {
        //cout<<"particle = "<<i+1<<" xA="<<xA[i]<<" yA="<<yA[i]<<" zA="<<zA[i]<<"xB="<<xB[i]<<" yB="<<yB[i]<<" zB="<<zB[i]<<endl;
    }*/

//******************************************************************************************************************************************************
    //calculating COM of A
    double xCOMA = xsumA/A;
    double yCOMA = ysumA/A;

    //finding coordinates wrt COM of A
    double XA[A], YA[A];
    double XB[A], YB[A];
        
    for (int i = 0; i<A; i++)
    {
        XA[i] = xA[i] - xCOMA + b/2;
        YA[i] = yA[i] - yCOMA;
    
        XB[i] = xB[i] - xCOMA - b/2;
        YB[i] = yB[i] - yCOMA;
    }

//******************************************************************************************************************************************************
    /*for (int i = 0; i<A; i++)
    {
        cout<<"particle = "<<i+1<<" xA="<<XA[i]<<" yA="<<YA[i]<<" zA="<<ZA[i]<<" xB="<<XB[i]<<" yB="<<YB[i]<<" zB="<<ZB[i]<<endl;
    }
    */
//******************************************************************************************************************************************************

    //Calculating number of collisions and number of participants.
    int NPartA = 0, NPartB = 0, NColl = 0;

    double avgx2 = 0, avgy2 = 0, meanx = 0, meany = 0;

    for (int i = 0; i < A; i++)
    {
        int ncollisionsA = 0;

        for (int j = 0; j<A; j++)
        {
            double rr = sqrt( (XA[i]-XB[j])*(XA[i]-XB[j]) + (YA[i]-YB[j])*(YA[i]-YB[j]) );
            
            if(rr <= d0) 
            {ncollisionsA++;}
        }
        NColl = NColl + ncollisionsA;
        if (ncollisionsA>0) 
        {
            NPartA++;

            avgx2 = avgx2 + (XA[i]*XA[i]);
            avgy2 = avgy2 + (YA[i]*YA[i]);
        }
    }

    for (int i = 0; i<A; i++)
    {
        int ncollisionsB = 0;
        for (int j = 0; j<A; j++)
        {
            double rr = sqrt( (XA[j]-XB[i])*(XA[j]-XB[i]) + (YA[j]-YB[i])*(YA[j]-YB[i]) );   
            if(rr <= d0) {ncollisionsB++;}
        }
        if (ncollisionsB>0) 
        {
            NPartB++;

            avgx2 = avgx2 + (XB[i]*XB[i]);
            avgy2 = avgy2 + (YB[i]*YB[i]);
        }
        //NColl = NColl + ncollisionsB;
    }

    int NPart = NPartA + NPartB;

    double sigmasqrnum = (avgy2) - (avgx2);
    double sigmasqrden = (avgy2) + (avgx2);

    double ecc = sigmasqrnum/sigmasqrden;

    //cout<<"Eccentricity = "<<ecc<<endl;

    rawdata[0] = NPart;
    rawdata[1] = NColl;
    rawdata[2] = ecc;

    //cout<<"Eccentricity = "<<rawdata[2]<<endl;

    //cout<<"b = "<<b<<"\tNPartA = "<<NPartA<<"\tNPartB = "<<NPartB<<"\tnpart = "<<NPart<<"\tncoll1 = "<<NColl<<"\tNColl2 = "<<NColl2<<endl;
    return rawdata;
}

void generateEvent()
{
    double b, Ecc;
    int npart, ncoll;

    TFile *file = new TFile("test_eventdata100000Au-Au0.2TeV_.root", "RECREATE");
    
    TTree *tree = new TTree("tree","eventdata");
    tree->Branch("b",&b, "b/D");
    tree->Branch("npart", &npart, "npart/I");
    tree->Branch("ncoll", &ncoll, "ncoll/I");
    tree->Branch("Ecc", &Ecc,"Ecc/D");

    double fBMin = 0., fBMax = 20.;

    for (int i = 0; i < NEvent; i++)
    {
        //b = i*(15./NEvent);
        b = TMath::Sqrt((fBMax*fBMax-fBMin*fBMin )  * rand1->Rndm() + fBMin*fBMin );
        //b = 15. * rand2->Rndm();
        //cout<<"\n"<<i+1<<endl;

        double *Data = new double[3];
        Data = event(b);
        npart = Data[0];
        ncoll = Data[1];
        Ecc = Data[2];

        //cout<<"b = "<<b<<"\tnpart = "<<npart<<"\tncoll = "<<ncoll<<endl;

        tree->Fill();
        if (i%1000 == 0)
        {
            cout<<"\t Events generated: "<<i<<"\r"<<flush;
        }
    }
    cout<<"\n"<<endl;
    file->Write();
    file->Close();
}

//do your work


















