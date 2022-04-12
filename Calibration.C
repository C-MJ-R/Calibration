#include <iostream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <vector>
#include <sstream>
#include <fstream>

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TChain.h>
#include <TLegend.h>
#include <TCut.h>
#include <TROOT.h>

#include "TCanvas.h"
#include "TMath.h" 
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"


void Calibration(int nbin=300, int nbin3 = 16 , int nmin= -2000, int nmax=20000, double nmin2 = 0, double nmax2 =16, int nmin3 = 0, int nmax3 = 16)
{   
    //plotting ch_roi
    int numrun = 4;
    int run[numrun];
    double recomean[numrun];
    double mcmean[numrun];
    double reco[10000];
    double truth[10000];
    for (int lk = 1; lk < numrun+1;lk++)
    {
        run[lk] = lk;
    }

    TCanvas *c1[numrun];
    TCanvas *c2[numrun];
    TCanvas *c3[numrun];
    TCanvas *c4[numrun];
    for (int l = 1; l < numrun+1; l++)
    {
        printf("Analysing ch_roi with NPE%d\n", run[l]);
        TFile *f = new TFile(Form("outputNPE%i.root",run[l]), "read");
        TTree *data = (TTree*)f->Get("dstree");
        
        //plotting ch-roi from branch
        c1[l] = new TCanvas(Form("c1%d",run[l]), "Finger Plot",200,10,600,400);
        TH1 *h1 = new TH1F("h1",Form("ch_roi_NPE%d",run[l]), nbin, nmin, nmax);
        data->Draw("ch_roi>>h1");
        double FirstBin = h1->FindFirstBinAbove(0,1,1,-1);
        double LastBin = h1->FindLastBinAbove(0,1,1,-1);
        nmin = h1->GetXaxis()->GetBinCenter(FirstBin);
        nmax = h1->GetXaxis()->GetBinCenter(LastBin);
        h1->GetXaxis()->SetRangeUser(nmin,nmax);
        h1->Draw();
        nmin = -2000;
        nmax = 20000;

        //Peak finding ch_roi
        TSpectrum *s1 = new TSpectrum(11);
        int nfound1 = s1->Search(h1,4,"",0.005);
        printf("Found %d candidate peaks to fit\n", nfound1);
        double *x1 = s1->GetPositionX();
        vector<double> x(x1, x1 + nfound1);
        sort(x.begin(), x.end());
        //mean distance ch_roi
        float dmu = 0;
        for (int i = 0; i < nfound1-1; i++)
        {
           dmu += x[i + 1] - x[i];
        }
        dmu = dmu/(nfound1-1);

        //Fitting single guas on ch_roi
        TF1 *g[nfound1];
        for (int p=0;p<nfound1;p++) 
        {
            g[p] = new TF1("gaus","gaus",x[p] -dmu/2, x[p] + dmu/2);
	
            g[p]->SetLineWidth(2);
            g[p]->SetLineColor(kRed);
            h1->Fit(g[p],"R+Q+0");
        }
        //Fitting sum of guas on ch_roi
        string sgaus = "gaus(0) ";
        for (int s = 1; s < nfound1; s++)
        {
                sgaus += Form("+ gaus(%d) ", 3*s);
        } 
      
        TF1 *sum = new TF1("gaussum",sgaus.c_str(),x[0] - dmu/2, x[nfound1 - 1] + dmu/2);
        sum->SetNpx(1000);

        for (int k=0;k<3*nfound1;k++)
        {
	        sum->SetParameter(k,g[(k-k%3)/3]->GetParameter(k%3));
	        if(!(k-1)%3) sum->FixParameter(k,g[(k-k%3)/3]->GetParameter(k%3));
            //look in the range +/- dmu/3
	        if(!(k-1)%3) sum->SetParLimits(k, sum->GetParameter(k) - dmu/3,sum->GetParameter(k) + dmu/3);    
        }
        h1->Fit(sum,"R+Q+0");

            // this refines the fit, fixes the other paramters and refines the x
            for (int k=0;k<3*nfound1;k++)
        {
	        sum->SetParameter(k,g[(k-k%3)/3]->GetParameter(k%3));
	        if(!(k-1)%3) sum->ReleaseParameter(k);
            //look in the range =/- dmu/3
            if(!(k-1)%3) sum->SetParLimits(k, sum->GetParameter(k) - dmu/3,sum->GetParameter(k) + dmu/3);
        }
      
        sum->SetLineWidth(2);
        sum->SetLineColor(kBlack);
        h1->Fit(sum,"R+Q");

        //extracting fit parameters and finding average distance
        vector<double>peakmean;
        for (int w = 0; w< nfound1; w++)
        {
            peakmean.push_back(sum->GetParameter((3*w)+1));
        }
        sort(peakmean.begin(), peakmean.end());
        float dmuf = 0;
        for (int ww = 0; ww <nfound1 -1; ww++)
        {
            dmuf += peakmean[ww + 1] - peakmean[ww];
        }
        /////////////////////////////////////////////////
        //finding the average charge/PE
        double CpPE[nfound1];
        double Cal = 0;
        for (int i = 1; i < nfound1; i++)
        {  
            CpPE[i] = x[i]/(i);
            Cal += CpPE[i];
        }
        Cal = Cal/(nfound1-1);
        printf("The Calibration for this run is %f ADU/PE\n",Cal); 
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
        //plotting ch_roi as a function of NPE;
        
        //c2[l] = new TCanvas(Form("c2%d",run[l]), "Ch_roi as a function of PE",200,10,600,400);
        TH1 *h2 = new TH1F("h2",Form("ch_roi_NPE%d",run[l]), nbin3, nmin2, nmax2);
        
        TTreeReader myReader("dstree", f);
        TTreeReaderValue<std::vector<float>> myVectorRV(myReader,"ch_roi");

        unsigned int evtCounter = 0;
        while (myReader.Next()){
            //cout << "Event " << evtCounter++ << endl;
            for (auto&& value : *myVectorRV)
            {
                h2->Fill(value/1400.0); //Eleason - change this?
            }   
        }



        recomean[l] = h2->GetMean();
        cout << recomean[l] << endl;
        //h2->Draw();
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Loading MCNPE
        c3[l] = new TCanvas(Form("c3%d",run[l]), "MC NPE per event as a function of PE",200,10,600,400);
        TH1 *h3 = new TH1F("h3",Form("MC_NPE%d",run[l]), nbin3, nmin3, nmax3);
        data->Draw("mc_npe>>h3");
        mcmean[l] = h3->GetMean();
        cout << mcmean[l] << endl;
        h3->SetFillColor(kGreen);
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //overlaying histograms
        THStack *hs = new THStack("hs","MC and reconstrcuted number of PE");
        hs->Add(h3);
        hs->Add(h2);
        hs->Draw("nostack");
        hs->SetTitle(Form("MC and reconstrcuted number of PE for %dNPE", run[l]));
        hs->GetXaxis()->SetTitle("Number of Photoelectrons");
        hs->GetYaxis()->SetTitle("Frequency");
        TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
        leg->AddEntry(h2,"Reconstruction NPE");
        leg->AddEntry(h3,"MC NPE");
        leg->Draw();


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // delete 3 (index 2)

        for (int k = 0; k < sizeof(CpPE)/sizeof(CpPE[0]); k++)
        {
            CpPE[k] = CpPE[k + 1]; // copy next element left
        }
        
        c4[l] = new TCanvas(Form("c4%d",run[l]), "Calibration Curve",200,10,600,400);
        Int_t p = nfound1;
        Double_t xc[p], y[p];
        for (Int_t i=0; i<p; i++) 
        {
            xc[i] = 1.0*(i+1);
            y[i] = 1.0*CpPE[i];
        }
        TGraph *gr1 = new TGraph(p,xc,y);
        gr1->SetTitle(Form("Calibration Curve for %dNPE w/no secondary events",run[l]));
        gr1->GetXaxis()->SetTitle("Number of Photoelectrons");
        gr1->GetYaxis()->SetTitle("Charge/Photoelectron");
        gr1->Draw("A*");
        
      /////////////////////////////////////////////////////////////////////////////
      //OLD way to calculate number of PE  
        //CpPE = (CpPE/(nfound1-1))*1e-13;
        // printf("%g\n",CpPE);

        //dmuf = (dmuf/(nfound1-1))*1e-13;
        //printf("ch_roi 1PE corresponds to %g C\n", dmuf);

        //Finding NPE for ch_roi
        //float TCharge = (h1->Integral((peakmean[0]+dmuf/2),LastBin))*1e-13;
        //double NPEchroi = TCharge/dmuf;
       //printf("NPE for ch_roi is: %g\n", NPEchroi);
       ////////////////////////////////////////////////////////////////////////////
    }


    //Reconstruction calibration as a function of NPE generated
    TCanvas *c5;
    c5 = new TCanvas("c5", "Calibration Curve",200,10,600,400);
        Int_t j = numrun;
        Double_t xj[j], yj[j];
        for (Int_t i=1; i<j+1; i++) 
        {
            xj[i] = 1.0*(i);
            yj[i] = 1.0*(recomean[i]/mcmean[i]);
            cout << yj[i] << endl;
        }
        
        TGraph *gr2 = new TGraph(j+1,xj,yj);
        gr2->SetTitle("Reconstruction Calibration");
        gr2->GetXaxis()->SetTitle("Number of Photoelectrons");
        gr2->GetYaxis()->SetTitle("Reconstruction Calibration");
        gr2->Draw("A*");



    TCanvas *c6;
    c6 = new TCanvas("c6", "Calibration Curve", 200, 10, 600, 400);
    TH2F *h = new TH2F("h","",nbin3,nmin3,nmax3,100,nmin3,nmax3);
   
    //Open file containing tree

    TFile *f = new TFile("outputNPE5.root", "read");
       
    
   // Create a TTreeReader for the tree, for instance by passing the
   // TTree's name and the TDirectory / TFile it is in.
   TTreeReader myReader("dstree", f);

    // The branch "px" contains floats; access them as myPx.
   TTreeReaderValue<Int_t> myPx(myReader, "mc_npe");
   // The branch "py" contains floats, too; access those as myPy.
   
   TTreeReaderValue<std::vector<float>> myVectorRV(myReader,"ch_roi");
   // Loop over all entries of the TTree or TChain.
   unsigned int evtCounter = 0;
        while (myReader.Next()){
            //cout << "Event " << evtCounter++ << endl;
            for (auto&& value : *myVectorRV)
            {
                Double_t Px = *myPx;
                h->Fill(Px,value/1400.0); //Eleason - change this?
            }   
        }
    h->Draw("COLZ");
    h->GetXaxis()->SetTitle("True NPE");
    h->GetYaxis()->SetTitle("Reco NPE");
}
    



  
   



