#include <iostream>
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


void Calibration(int nbin=300, int nmin= -2000, int nmax=20000, double nmin2 = 0, double nmax2 =6)
{   
    //plotting ch_roi
    int numrun = 5;
    int run[numrun];
    vector<float> histentry;
    for (int lk = 1; lk < numrun+1;lk++)
    {
        run[lk] = lk;
    }
    TCanvas *c1[numrun];
    TCanvas *c2[numrun];
    for (int l = 1; l < numrun+1; l++)
    {
        printf("Analysing ch_roi with NPE%d\n", run[l]);
        TFile *f = new TFile(Form("outputNPE%i.root",run[l]), "read");
        TTree *data = (TTree*)f->Get("dstree");
    
        c1[l] = new TCanvas(Form("c1%d",run[l]), "Finger Plot",200,10,600,400);
        TH1 *h1 = new TH1F("h1",Form("ch_roi_NPE%d",run[l]), nbin, nmin, nmax);
        data->Draw("ch_roi>>h1");
        double FirstBin = h1->FindFirstBinAbove(0,1,1,-1);
        double LastBin = h1->FindLastBinAbove(0,1,1,-1);
        //nmin = h1->GetXaxis()->GetBinCenter(FirstBin);
        //nmax = h1->GetXaxis()->GetBinCenter(LastBin);
        //h1->GetXaxis()->SetRangeUser(nmin,nmax);
        h1->Draw();

        //Peak finding ch_roi
        TSpectrum *s1 = new TSpectrum(10);
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
        for (int i = 1; i < nfound1; i++)
        {   
            CpPE[i] = x[i]/(i);
            printf("%f\n",CpPE[i]);
        }
        // delete 3 (index 2)
        for (int k = 0; k < sizeof(CpPE)/sizeof(CpPE[0]); k++)
        {
            CpPE[k] = CpPE[k + 1]; // copy next element left
        }
        
        c2[l] = new TCanvas(Form("c2%d",run[l]), "Calibration Curve",200,10,600,400);
        Int_t p = nfound1;
        Double_t xc[p], y[p];
        for (Int_t i=0; i<p; i++) 
        {
            xc[i] = 1.0*(i+1);
            y[i] = 1.0*CpPE[i];
        }
        TGraph *gr = new TGraph(p,xc,y);
        gr->SetTitle(Form("Calibration Curve for %dNPE w/no secondary events",run[l]));
        gr->GetXaxis()->SetTitle("Number of Photoelectrons");
        gr->GetYaxis()->SetTitle("Charge/Photoelectron");
        gr->Draw("A*");
        
        //CpPE = (CpPE/(nfound1-1))*1e-13;
        // printf("%g\n",CpPE);

       dmuf = (dmuf/(nfound1-1))*1e-13;
       printf("ch_roi 1PE corresponds to %g C\n", dmuf);

        //Finding NPE for ch_roi
        float TCharge = (h1->Integral((peakmean[0]+dmuf/2),LastBin))*1e-13;
        double NPEchroi = TCharge/dmuf;
        printf("NPE for ch_roi is: %g\n", NPEchroi);
    }

}
//     printf("Analysing pk_p\n");
//     //plotting pk_p
//     TCanvas *c2;

//         c2 = new TCanvas("c2", "Finger Plot",200,10,600,400);

//         TH1 *h2 = new TH1F("h2","pk_p",nbin, nmin2, nmax2);
//         data->Draw("pk_p>>h2");
//         h2->Draw();

//     //Peak finding pk_p
//     TSpectrum *s2 = new TSpectrum(10);
//     int nfound2 = s2->Search(h2,4,"",0.005);
//     printf("Found %d candidate peaks to fit\n", nfound2);
//     double *x2 = s2->GetPositionX();
//     vector<double> X(x2, x2 + nfound2);
//     sort(X.begin(), X.end());
//     float dmu2 = 0;
//     for (int ii = 0; ii < nfound2-1; ii++)
//     {
//         dmu2 += X[ii + 1] - X[ii];
//     }

//     //mean distance pk_p
//     dmu2 = dmu2/(nfound2-1);

//     //Fitting single guas on pk_p
//     TF1 *g2[nfound2];
//     for (int pp=0;pp<nfound2;pp++) 
//     {
//         g2[pp] = new TF1("gaus2","gaus",X[pp] -dmu2/2, X[pp] + dmu2/2);
	
//         g2[pp]->SetLineWidth(2);
//         g2[pp]->SetLineColor(kRed);
//         h2->Fit(g2[pp],"R+Q+0");
//     }
//     //Fitting sum of guas on pk_p
//     string sgaus2 = "gaus(0) ";
//       for (int ss = 1; ss < nfound2; ss++)
//       {
//             sgaus2 += Form("+ gaus(%d) ", 3*ss);
//       } 
      
//       TF1 *sum2 = new TF1("gaussum2",sgaus2.c_str(),X[0] - dmu2/2, X[nfound2 - 1] + dmu2/2);
//       sum2->SetNpx(1000);
      
//       for (int kk=0;kk<3*nfound2;kk++)
//       {
// 	        sum2->SetParameter(kk,g2[(kk-kk%3)/3]->GetParameter(kk%3));
// 	        if(!(kk-1)%3) sum2->FixParameter(kk,g2[(kk-kk%3)/3]->GetParameter(kk%3));
//             //look in the range +/- dmu/3
// 	        if(!(kk-1)%3) sum2->SetParLimits(kk, sum2->GetParameter(kk) - dmu2/3,sum2->GetParameter(kk) + dmu2/3);
//       }
    
//        // this refines the fit, fixes the other paramters and refines the x
//        for (int kk=0;kk<3*nfound2;kk++)
//       {
// 	        sum2->SetParameter(kk,g2[(kk-kk%3)/3]->GetParameter(kk%3));
// 	        if(!(kk-1)%3) sum2->ReleaseParameter(kk);
//             //look in the range =/- dmu/3
//             if(!(kk-1)%3) sum2->SetParLimits(kk, sum2->GetParameter(kk) - dmu2/3,sum2->GetParameter(kk) + dmu2/3);
//       }
      
//        sum2->SetLineWidth(2);
//        sum2->SetLineColor(kBlack);
//        h2->Fit(sum2,"R+Q");


    
//     //extracting fit parameters and finding average distance pk_p
//     vector<double>peakmean2;
//     for (int l = 0; l< nfound2; l++)
//     {
//         peakmean2.push_back(sum2->GetParameter((3*l)+1));
//     }
//     sort(peakmean2.begin(), peakmean2.end());
//     float dmuf2 = 0;
//     for (int ll = 0; ll <nfound2 -1; ll++)
//     {
//         dmuf2 += peakmean2[ll + 1] - peakmean2[ll];
//     }

//     dmuf2 = (dmuf2/(nfound2-1));
//     printf("pk_p 1PE corresponds to %g C \n", dmuf2);

//     //Finding NPE pk_p
//     double Tpkp = 0;
//     for (int i = 0; i < nfound2 -1; i++)
//     {
//         Tpkp += peakmean2[i];
//     }
//     double NPEpkp = Tpkp/peakmean[1];
    
//     printf("NPE for pk_p is: %g\n", NPEpkp);
// }