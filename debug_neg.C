void debug_neg() {
    TFile *f = TFile::Open("tagged-neutron-DIS/TaggedNeutron-DIS-EicC.root");
    TTree *t = (TTree*)f->Get("tree");
    double d4sigma;
    t->SetBranchAddress("d4sigma", &d4sigma);
    int neg10=0, neg20=0, neg30=0, pos=0;
    for(int i=0; i<t->GetEntries(); i++){
        t->GetEntry(i);
        if(d4sigma > 0) pos++;
        else if(d4sigma > -15) neg10++;
        else if(d4sigma > -25) neg20++;
        else neg30++;
    }
    cout << "Positive: " << pos << endl;
    cout << "d4sigma=-10 (sigma<0): " << neg10 << endl;
    cout << "d4sigma=-20 (NaN): " << neg20 << endl;
    cout << "d4sigma=-30 (xpi>0.999): " << neg30 << endl;
}
