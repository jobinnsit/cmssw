void CPECPU_plot() {
   //Draw a simple graph

   //TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
    TCanvas  *c1 = new TCanvas("c1", "c1",800,800);
   c1->SetFillColor(42);
   c1->SetGrid();

   const Int_t n = 101;
   Double_t x[n], y[n];
   // x is no of pixel and y is time
   // read the file
   FILE *fp;
   //fp = fopen( "RawToDigiCPUTime.txt", "r" );
   fp  = fopen("CPE_CPU_Time.txt", "r");
   Int_t i=0, event =0,xc=0, yc=0;
   char line[100];
   fgets(line, sizeof line, fp);

while(!feof(fp)){

     fscanf(fp,"%i %i %i",&event, &xc, &yc );
     y[i] = yc;//1000;
     x[i] = xc;
     i++;
   }
   fclose(fp);
   TGraph *gr = new TGraph(n,x,y);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->SetTitle("CPE CPU Time vs Number of Clusters");
   gr->GetXaxis()->SetTitle("Number of clusters");
   gr->GetYaxis()->SetTitle("CPU Time(#mu seconds)");
   gr->Draw("AP");

   // TCanvas::Update() draws the frame, after which one can change it
   c1->Update();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   c1->Modified();
}
