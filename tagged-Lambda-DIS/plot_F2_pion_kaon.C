#include "piIMParton.h"
#include "piIMParton.cpp"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TLine.h"

double calcF2(piIMParton &pdf, double x, double Q2) {
    double eu2=4./9., ed2=1./9., es2=1./9., ec2=4./9.;
    double u=pdf.getPDF(1,x,Q2), ub=pdf.getPDF(-1,x,Q2);
    double d=pdf.getPDF(2,x,Q2), db=pdf.getPDF(-2,x,Q2);
    double s=pdf.getPDF(3,x,Q2), sb=pdf.getPDF(-3,x,Q2);
    double c=pdf.getPDF(4,x,Q2), cb=pdf.getPDF(-4,x,Q2);
    return x*(eu2*(u+ub)+ed2*(d+db)+es2*(s+sb)+ec2*(c+cb));
}

void plot_F2_pion_kaon() {
    gStyle->SetOptStat(0);
    piIMParton pdf;
    int nPts=200;
    double Q2vals[]={2,5,10,20,50}; int nQ2=5;
    int cPi[]={kRed+1,kOrange+7,kMagenta+1,kViolet+1,kPink+1};
    int cK[]={kBlue+1,kCyan+2,kTeal+3,kAzure+2,kSpring+4};

    // ===== Canvas 1: F2 side by side =====
    TCanvas *c1=new TCanvas("c1","",1400,600); c1->Divide(2,1);
    c1->cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.06);
    TMultiGraph *mg1=new TMultiGraph(); mg1->SetTitle(";x;F_{2}^{#pi}(x, Q^{2})");
    TLegend *l1=new TLegend(0.55,0.55,0.92,0.90); l1->SetBorderSize(1); l1->SetTextSize(0.038);
    l1->SetHeader("F_{2}^{#pi^{+}} (IMParton)");
    pdf.setDataSet(211);
    for(int iq=0;iq<nQ2;iq++){
        TGraph *g=new TGraph();
        for(int i=0;i<nPts;i++){double x=0.001+i*0.998/nPts; double f=calcF2(pdf,x,Q2vals[iq]); if(f>0&&f==f) g->SetPoint(g->GetN(),x,f);}
        g->SetLineColor(cPi[iq]); g->SetLineWidth(3); mg1->Add(g,"L");
        l1->AddEntry(g,Form("Q^{2} = %.0f GeV^{2}",Q2vals[iq]),"l");
    }
    mg1->Draw("A"); mg1->GetXaxis()->SetLimits(0,1); mg1->SetMinimum(0); mg1->SetMaximum(0.55);
    mg1->GetXaxis()->SetTitleSize(0.055); mg1->GetYaxis()->SetTitleSize(0.055); l1->Draw();
    TLatex *t1=new TLatex(); t1->SetNDC(); t1->SetTextSize(0.055); t1->SetTextColor(kRed+1); t1->DrawLatex(0.18,0.88,"(a)");

    c1->cd(2); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.06);
    TMultiGraph *mg2=new TMultiGraph(); mg2->SetTitle(";x;F_{2}^{K}(x, Q^{2})");
    TLegend *l2=new TLegend(0.55,0.55,0.92,0.90); l2->SetBorderSize(1); l2->SetTextSize(0.038);
    l2->SetHeader("F_{2}^{K^{+}} (IMParton)");
    pdf.setDataSet(321);
    for(int iq=0;iq<nQ2;iq++){
        TGraph *g=new TGraph();
        for(int i=0;i<nPts;i++){double x=0.001+i*0.998/nPts; double f=calcF2(pdf,x,Q2vals[iq]); if(f>0&&f==f) g->SetPoint(g->GetN(),x,f);}
        g->SetLineColor(cK[iq]); g->SetLineWidth(3); mg2->Add(g,"L");
        l2->AddEntry(g,Form("Q^{2} = %.0f GeV^{2}",Q2vals[iq]),"l");
    }
    mg2->Draw("A"); mg2->GetXaxis()->SetLimits(0,1); mg2->SetMinimum(0); mg2->SetMaximum(0.55);
    mg2->GetXaxis()->SetTitleSize(0.055); mg2->GetYaxis()->SetTitleSize(0.055); l2->Draw();
    TLatex *t2=new TLatex(); t2->SetNDC(); t2->SetTextSize(0.055); t2->SetTextColor(kBlue+1); t2->DrawLatex(0.18,0.88,"(b)");
    c1->Print("F2_pion_kaon.png");

    // ===== Canvas 2: xu_v vs Drell-Yan data from paper Fig 9 =====
    TCanvas *c2=new TCanvas("c2","",800,600);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.06);

    // Drell-Yan Data (Q2 in [16,49] GeV2) - digitized from paper Figure 9, Ref [94]
    const int nDY=13;
    double xDY[] ={0.21, 0.28, 0.34, 0.40, 0.46, 0.52, 0.58, 0.64, 0.70, 0.76, 0.82, 0.88, 0.94};
    double yDY[] ={0.33, 0.36, 0.36, 0.35, 0.33, 0.30, 0.26, 0.22, 0.17, 0.13, 0.09, 0.05, 0.02};
    double eyDY[]={0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.015,0.01, 0.008,0.005};
    double exDY[nDY]; for(int i=0;i<nDY;i++) exDY[i]=0;

    TGraphErrors *gDY=new TGraphErrors(nDY,xDY,yDY,exDY,eyDY);
    gDY->SetMarkerStyle(20); gDY->SetMarkerSize(1.0);
    gDY->SetMarkerColor(kBlack); gDY->SetLineColor(kBlack); gDY->SetLineWidth(2);

    TMultiGraph *mgV=new TMultiGraph(); mgV->SetTitle(";x_{#pi};x_{#pi} u_{v}(x_{#pi})");

    // IMParton curves at Q2=27 GeV2 (middle of [16,49] range)
    double Q2m=27.0;
    TGraph *gGlob=new TGraph(); pdf.setDataSet(211);
    for(int i=0;i<nPts;i++){double x=0.001+i*0.998/nPts; double v=pdf.getXUV(x,Q2m); if(v>=0&&v==v) gGlob->SetPoint(gGlob->GetN(),x,v);}
    gGlob->SetLineColor(kRed+1); gGlob->SetLineWidth(3);

    TGraph *gA=new TGraph(); pdf.setDataSet(1);
    for(int i=0;i<nPts;i++){double x=0.001+i*0.998/nPts; double v=pdf.getXUV(x,Q2m); if(v>=0&&v==v) gA->SetPoint(gA->GetN(),x,v);}
    gA->SetLineColor(kBlue+1); gA->SetLineWidth(3); gA->SetLineStyle(2);

    TGraph *gB=new TGraph(); pdf.setDataSet(2);
    for(int i=0;i<nPts;i++){double x=0.001+i*0.998/nPts; double v=pdf.getXUV(x,Q2m); if(v>=0&&v==v) gB->SetPoint(gB->GetN(),x,v);}
    gB->SetLineColor(kGreen+2); gB->SetLineWidth(3); gB->SetLineStyle(7);

    mgV->Add(gDY,"PE"); mgV->Add(gGlob,"L"); mgV->Add(gA,"L"); mgV->Add(gB,"L");
    mgV->Draw("A");
    mgV->GetXaxis()->SetTitleSize(0.055); mgV->GetYaxis()->SetTitleSize(0.055);
    mgV->GetXaxis()->SetLabelSize(0.045); mgV->GetYaxis()->SetLabelSize(0.045);
    mgV->GetXaxis()->SetLimits(0,1); mgV->SetMinimum(0); mgV->SetMaximum(0.50);

    TLegend *lV=new TLegend(0.50,0.68,0.93,0.93);
    lV->SetBorderSize(1); lV->SetTextSize(0.035);
    lV->AddEntry(gDY,"Drell-Yan Data","pe");
    lV->AddEntry(gGlob,"IMParton Global Fit","l");
    lV->AddEntry(gA,"IMParton Set A","l");
    lV->AddEntry(gB,"IMParton Set B","l");
    lV->Draw();

    TLatex *tV=new TLatex(); tV->SetNDC(); tV->SetTextSize(0.038);
    tV->DrawLatex(0.16,0.88,"x_{#pi}u_{v}^{#pi}(x_{#pi})");
    tV->SetTextSize(0.032);
    tV->DrawLatex(0.16,0.82,"Q^{2} #in [16, 49] GeV^{2}");
    c2->Print("xuv_pion_vs_E615.png");

    // ===== Canvas 3: F2 pion vs kaon overlay =====
    TCanvas *c3=new TCanvas("c3","",800,600);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.06);
    TMultiGraph *mg3=new TMultiGraph(); mg3->SetTitle(";x;F_{2}^{M}(x, Q^{2})");
    double Q2c[]={5,10,27}; int nC=3;
    int cpC[]={kRed+1,kOrange+7,kMagenta+1}, ckC[]={kBlue+1,kCyan+2,kAzure+2};
    TLegend *l3=new TLegend(0.55,0.55,0.93,0.92); l3->SetBorderSize(1); l3->SetTextSize(0.033);
    l3->SetHeader("IMParton Global Fit");
    for(int iq=0;iq<nC;iq++){
        TGraph *gp=new TGraph(); pdf.setDataSet(211);
        for(int i=0;i<nPts;i++){double x=0.001+i*0.998/nPts; double f=calcF2(pdf,x,Q2c[iq]); if(f>0&&f==f) gp->SetPoint(gp->GetN(),x,f);}
        gp->SetLineColor(cpC[iq]); gp->SetLineWidth(3); mg3->Add(gp,"L");
        l3->AddEntry(gp,Form("F_{2}^{#pi}, Q^{2}=%.0f",Q2c[iq]),"l");
        TGraph *gk=new TGraph(); pdf.setDataSet(321);
        for(int i=0;i<nPts;i++){double x=0.001+i*0.998/nPts; double f=calcF2(pdf,x,Q2c[iq]); if(f>0&&f==f) gk->SetPoint(gk->GetN(),x,f);}
        gk->SetLineColor(ckC[iq]); gk->SetLineWidth(3); gk->SetLineStyle(2); mg3->Add(gk,"L");
        l3->AddEntry(gk,Form("F_{2}^{K}, Q^{2}=%.0f",Q2c[iq]),"l");
    }
    mg3->Draw("A"); mg3->GetXaxis()->SetLimits(0,1); mg3->SetMinimum(0); mg3->SetMaximum(0.55);
    mg3->GetXaxis()->SetTitleSize(0.055); mg3->GetYaxis()->SetTitleSize(0.055); l3->Draw();
    TLatex *t3=new TLatex(); t3->SetNDC(); t3->SetTextSize(0.040);
    t3->DrawLatex(0.16,0.88,"F_{2}^{#pi} (solid) vs F_{2}^{K} (dashed)");
    c3->Print("F2_pion_vs_kaon.png");

    // ===== Canvas 4: Ratio =====
    TCanvas *c4=new TCanvas("c4","",800,600);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.06);
    TMultiGraph *mg4=new TMultiGraph(); mg4->SetTitle(";x;F_{2}^{K} / F_{2}^{#pi}");
    TLegend *l4=new TLegend(0.60,0.62,0.92,0.90); l4->SetBorderSize(1); l4->SetTextSize(0.035);
    int cR[]={kBlack,kRed+1,kBlue+1,kGreen+2,kMagenta+1};
    for(int iq=0;iq<nQ2;iq++){
        TGraph *gr=new TGraph();
        for(int i=1;i<nPts;i++){double x=0.01+i*0.94/nPts;
            pdf.setDataSet(211); double fp=calcF2(pdf,x,Q2vals[iq]);
            pdf.setDataSet(321); double fk=calcF2(pdf,x,Q2vals[iq]);
            if(fp>1e-10) gr->SetPoint(gr->GetN(),x,fk/fp);}
        gr->SetLineColor(cR[iq]); gr->SetLineWidth(3); mg4->Add(gr,"L");
        l4->AddEntry(gr,Form("Q^{2} = %.0f GeV^{2}",Q2vals[iq]),"l");
    }
    mg4->Draw("A"); mg4->GetXaxis()->SetLimits(0,1); mg4->SetMinimum(0); mg4->SetMaximum(1.5);
    mg4->GetXaxis()->SetTitleSize(0.055); mg4->GetYaxis()->SetTitleSize(0.055); l4->Draw();
    TLine *lr=new TLine(0,1,1,1); lr->SetLineStyle(2); lr->SetLineWidth(2); lr->SetLineColor(kGray+2); lr->Draw();
    TLatex *t4=new TLatex(); t4->SetNDC(); t4->SetTextSize(0.040);
    t4->DrawLatex(0.18,0.88,"F_{2}^{K}/F_{2}^{#pi} Ratio");
    t4->DrawLatex(0.18,0.81,"IMParton Global Fit");
    c4->Print("F2_ratio_kaon_pion.png");

    cout<<"\n=== Plots Generated ==="<<endl;
}
