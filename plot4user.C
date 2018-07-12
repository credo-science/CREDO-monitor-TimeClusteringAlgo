using namespace std; 

void plot4user(string user_name)
{
  //------------------Extracting parameters from analysis output file-------------------
  vector<double> time;
  vector<double> meanTS_SC; 
  vector<double> testStatistic_data;
  vector<double> threeSig_TS;
  vector<double> fiveSig_TS;
  vector<double> meanNumOfDoublets_SC;
  vector<double> counterNumOfDoublets_Data;
  vector<double> threeSig_doublet;
  vector<double> fiveSig_doublet;

  int numOfPoints = 0;
  double a,b,c,d,e,f,g,h,i;

  ifstream input;
  char data_in [256];
  sprintf(data_in,"output_parameters_%s.txt",user_name.c_str());
  string line;

  input.open(data_in);
  if (!input)
 	{
   	cout << "Unable to open file." << endl;
   	return 0;
 	}

  while (input >> a >> b >> c >> d >> e >> f >> g >> h >> i)
 	{
 		time.push_back(a);
   	meanTS_SC.push_back(b);
    testStatistic_data.push_back(c);
   	threeSig_TS.push_back(d);
   	fiveSig_TS.push_back(e);
   	meanNumOfDoublets_SC.push_back(f);
    counterNumOfDoublets_Data.push_back(g);
    threeSig_doublet.push_back(h);
    fiveSig_doublet.push_back(i);

    numOfPoints++;   		     	   	
  }
  input.close();

  //-------------------------------------------------------------------------------

  //------------------------------------PLOT---------------------------------------
	TString plot_name = TString::Format("plot_user-%s.root",user_name.c_str());
	TFile file1(plot_name, "RECREATE");

  TString canvas_name = TString::Format("Live plot of smartphone data for user %s",user_name.c_str());
	TCanvas *c1 = new TCanvas("c1",canvas_name);
	c1->Divide(1,2);

  TGraph* gr2 = new TGraph(numOfPoints,&time[0],&meanTS_SC[0]);
  gr2->SetMarkerStyle(3);
  gr2->SetMarkerColor(kBlue+3);
  gr2->SetLineStyle(1);
  gr2->SetLineColor(kBlue+3);

  TGraph* gr3 = new TGraph(numOfPoints,&time[0],&testStatistic_data[0]);
  gr3->SetMarkerStyle(3);
  gr3->SetMarkerColor(2);
  gr3->SetLineStyle(1);
  gr3->SetLineColor(2);

  TGraph* gr4 = new TGraph(numOfPoints,&time[0],&threeSig_TS[0]);
  gr4->SetMarkerStyle(8);
  gr4->SetMarkerColor(kGreen+2);
  gr4->SetLineStyle(10);
  gr4->SetLineColor(kGreen+2);
  gr4->SetFillStyle(3002);
  gr4->SetFillColor(kGreen+2);

  TGraph* gr5 = new TGraph(numOfPoints,&time[0],&fiveSig_TS[0]);
  gr5->SetMarkerStyle(8);
  gr5->SetMarkerColor(7);
  gr5->SetLineStyle(10);
  gr5->SetLineColor(7);

  c1->cd(1);

  TMultiGraph *mg2 = new TMultiGraph();

  mg2->Add(gr2);
  mg2->Add(gr3);
	mg2->Add(gr4);
	mg2->Add(gr5);
  TString mg2_title = TString::Format("Test statistic and 3-5 sigma limits - User %s",user_name.c_str());
	mg2->SetTitle(mg2_title);
	mg2->Draw("ALP");
	mg2->GetXaxis()->SetTitle("Date [dd/mm/yy HH:MM]");
  mg2->GetXaxis()->SetTimeFormat("%d/%m/%y %H:%M %F1970-01-01 00:00:01");
  mg2->GetXaxis()->SetTimeDisplay(1);
	mg2->GetXaxis()->SetLabelSize(0.04);
	mg2->GetXaxis()->SetTitleSize(0.06);
	mg2->GetXaxis()->SetTitleOffset(0.75);
	mg2->GetYaxis()->SetTitle("log10(Doublet duration [sec])");
	mg2->GetYaxis()->SetLabelSize(0.05);
	mg2->GetYaxis()->SetTitleSize(0.06);
	mg2->GetYaxis()->SetTitleOffset(0.5);

	auto Legend2 = new TLegend(0.1,1,0.3,0.9);
  Legend2->AddEntry(gr2,"Mean test statictic of background (expected)","l");
  Legend2->AddEntry(gr3,"Test statistic of data (observed)","l");
  Legend2->Draw();

  auto Legend3 = new TLegend(0.8,1,0.9,0.9);
	Legend3->AddEntry(gr4,"3-sigma limit","l");
  Legend3->AddEntry(gr5,"5-sigma limit","l");
  Legend3->Draw();

	gPad->Modified();

  TGraph* gr6 = new TGraph(numOfPoints,&time[0],&meanNumOfDoublets_SC[0]);
  gr6->SetMarkerStyle(3);
  gr6->SetMarkerColor(kBlue+3);
  gr6->SetLineStyle(1);
  gr6->SetLineColor(kBlue+3);

  TGraph* gr7 = new TGraph(numOfPoints,&time[0],&counterNumOfDoublets_Data[0]);
  gr7->SetMarkerStyle(3);
  gr7->SetMarkerColor(2);
  gr7->SetLineStyle(1);
  gr7->SetLineColor(2);

  TGraph* gr8 = new TGraph(numOfPoints,&time[0],&threeSig_doublet[0]);
  gr8->SetMarkerStyle(8);
  gr8->SetMarkerColor(kGreen+2);
  gr8->SetLineStyle(10);
  gr8->SetLineColor(kGreen+2);
  gr8->SetFillStyle(3002);
  gr8->SetFillColor(kGreen+2);

  TGraph* gr9 = new TGraph(numOfPoints,&time[0],&fiveSig_doublet[0]);
  gr9->SetMarkerStyle(8);
  gr9->SetMarkerColor(7);
  gr9->SetLineStyle(10);
  gr9->SetLineColor(7);

  c1->cd(2);

  TMultiGraph *mg3 = new TMultiGraph();

  mg3->Add(gr6);
  mg3->Add(gr7);
  mg3->Add(gr8);
  mg3->Add(gr9);
  TString mg3_title = TString::Format("# of doublets in 5 min. time windows - User %s",user_name.c_str());
  mg3->SetTitle(mg3_title);  
  mg3->Draw("ACP");
  mg3->GetXaxis()->SetTitle("Date [dd/mm/yy HH:MM]");
  mg3->GetXaxis()->SetTimeFormat("%d/%m/%y %H:%M %F1970-01-01 00:00:01");
  mg3->GetXaxis()->SetTimeDisplay(1);  
  mg3->GetXaxis()->SetLabelSize(0.04);
  mg3->GetXaxis()->SetTitleSize(0.06);
  mg3->GetXaxis()->SetTitleOffset(0.75);
  mg3->GetYaxis()->SetTitle("Number of doublets");
  mg3->GetYaxis()->SetLabelSize(0.05);
  mg3->GetYaxis()->SetTitleSize(0.06);
  mg3->GetYaxis()->SetTitleOffset(0.5);

  auto Legend4 = new TLegend(0.1,1,0.3,0.9);
  Legend4->AddEntry(gr6,"Mean number of doublets of background (expected)","l");
  Legend4->AddEntry(gr7,"Number of doublets in data (observed)","l");
  Legend4->Draw();

  auto Legend5 = new TLegend(0.8,1,0.9,0.9);
  Legend5->AddEntry(gr8,"3-sigma limit","l");
  Legend5->AddEntry(gr9,"5-sigma limit","l");
  Legend5->Draw();

  gPad->Modified();

	c1->Write();

  TString plot_name_png = TString::Format("plot_user-%s.png",user_name.c_str());
  c1->Print(plot_name_png);
}
