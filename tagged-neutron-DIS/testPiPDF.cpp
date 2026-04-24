#include"piIMParton.h"

#include"piIMParton.cpp"


void testPiPDF(){

	piIMParton pipdf;

	TGraph *gxuv = new TGraph();
	TGraph *gxu = new TGraph();
	TGraph *gxubar = new TGraph();
	TGraph *gxgluon = new TGraph();

	double step = (log(1)-log(1e-3)) / 400.0;
	for(int i=0; i<395; i++){
		double xpi = exp( log(1e-3) + i*step );
		///get xuv at Q2 = 10 GeV2.
		gxuv->SetPoint(i,  xpi,  xpi*(pipdf.getPDF(1, xpi, 10.0) - pipdf.getPDF(-1, xpi, 10.0)) );
		///get xu at Q2 = 10 GeV2.
		gxu->SetPoint(i,  xpi,  xpi*pipdf.getPDF(1, xpi, 10.0) );
		///get xubar at Q2 = 10 GeV2.
		gxubar->SetPoint(i,  xpi,  xpi*pipdf.getPDF(-1, xpi, 10.0) );
		///get xgluon at Q2 = 10 GeV2.
		gxgluon->SetPoint(i,  xpi,  xpi*pipdf.getPDF(0, xpi, 10.0) );
	}


	/// plotting
	TCanvas *c = new TCanvas("c","c",800,600);
	c->Divide(2,2);
	c->cd(1);
	gxuv->Draw("a l");
	gxuv->SetTitle("xu_{V} of #pi^{+}");
	c->cd(2);
	gxu->Draw("a l");
	gxu->SetTitle("xu of #pi^{+}");
	c->cd(3);
	gxubar->Draw("a l");
	gxubar->SetTitle("x#bar{u} of #pi^{+}");
	c->cd(4);
	gxgluon->Draw("a l");
	gxgluon->SetTitle("xgluon of #pi^{+}");

}


