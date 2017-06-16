void CPEGPU_plot() {
   //Draw a simple graph

   //TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
    TCanvas  *c1 = new TCanvas("c1", "c1",800,800);
   c1->SetFillColor(42);
   c1->SetGrid();

   const Int_t n = 102;
   Double_t x[n], y[n];
   // x is no of pixel and y is time
   // read the file
   FILE *fp;
   //fp = fopen( "RawToDigiCPUTime.txt", "r" );
   fp  = fopen("CPE_GPU_Time.txt", "r");
   Int_t i=0, event =0,xc=0;
   Float_t yc=0;
   char line[120];
   fgets(line, sizeof line, fp);
   puts(line);
  // skip the first entry launching overhead problem
  fscanf(fp, "%i  %i  %f", &event, &xc, &yc);
  
  while(!feof(fp)){

     fscanf(fp,"%i %i %f",&event, &xc, &yc );
     //printf("%i  %i   %i\n", event, xc,  yc);
     y[i] = yc;
     x[i] = xc;
     i++;
     //printf("i %i\n", i);
   
   }
   fclose(fp);
   TGraph *gr = new TGraph(n,x,y);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->SetTitle("CPE GPU Time vs Number of Clusters");
   gr->GetXaxis()->SetTitle("Number of Cluster");
   gr->GetYaxis()->SetTitle("GPU Time(#mu Seconds)");
   gr->Draw("AP");

   // TCanvas::Update() draws the frame, after which one can change it
   c1->Update();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   c1->Modified();
}
