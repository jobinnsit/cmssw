void Cluster_CPU_GPU_gain() {
   //Draw a simple graph
   
   //TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
    TCanvas  *c1 = new TCanvas("c1", "c1",800,800);
   c1->SetFillColor(42);
   c1->SetGrid();

   const Int_t n = 101;
   Double_t x_ratio[n], y_ratio[n], y_cpu, x_cpu;
   // x is no of pixel and y is time
   // read the file
   FILE *fp_cpu, *f_gpu;
   fp_cpu = fopen( "Cluster_CPU_Time.txt", "r");
   fp_gpu = fopen( "Cluster_GPU_Time.txt", "r" );
   
   Double_t event =0,xc=0, yc=0;
   Int_t i =0;
   // read the header
   char line[100];
   fgets(line, sizeof line, fp_cpu);
   fgets(line, sizeof line, fp_gpu);
   //skip the first entry 
   fscanf(fp_cpu, "%lf %lf %lf", &event, &xc, &yc);
   fscanf(fp_gpu, "%lf %lf %lf", &event, &xc, &yc);

   while(!feof(fp_gpu)){

     fscanf(fp_cpu,"%lf %lf %lf",&event, &xc, &yc );
     y_cpu = yc;
     fscanf(fp_gpu,"%lf %lf %lf",&event, &xc, &yc);
     //yc = yc*1000;
     y_ratio[i] = y_cpu/yc;
     x_ratio[i] = xc;///1000;
     printf("time_cpu:  %lf   time_gpu:  %lf    no_pixels:  %lf    ratio:  %lf\n",y_cpu,yc,x_ratio[i], y_ratio[i]);
     i++;
     if(i==99) break; // to avoid reading last new line


}
   fclose(fp_cpu);
   fclose(fp_gpu);
   TGraph *gr = new TGraph(n,x_ratio,y_ratio);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->SetTitle("Cluster Gain plot CPU vs GPU ");
   gr->GetXaxis()->SetTitle("Number of Pixels");
   gr->GetYaxis()->SetTitle("Gain");
   gr->Draw("AP");

   // TCanvas::Update() draws the frame, after which one can change it
   c1->Update();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   c1->Modified();
}
