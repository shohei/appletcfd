/**=====================================================================
   Metric transformation of general curvilinear coordinate systems    
           Metric.Java  Metric for staggered mesh                    
                      All rights reserved, Copyright (C) 2001-2003,
           Ver.2.0    Last update: December 25, 2002, Kiyoshi Minemura
=======================================================================*/

public class Metric {
   public double XIxW[][], XIyW[][],  ETAxW[][], JW[][], ETAyW[][];
   public double XIyS[][], ETAxS[][], ETAyS[][], JS[][];
   public double D1u[][], D2u[][], D3u[][], D4u[][], D5u[][], D6u[][], D7u[][];
   public double D8u[][], D0u[][], D1v[][], D2v[][], D3v[][], D4v[][], D5v[][];
   public double D6v[][], D7v[][], D8v[][], D0v[][], D1P[][], D2P[][], D3P[][];
   public double D4P[][], D5P[][], D6P[][], D7P[][], D8P[][], D0P[][];  
   public double c1[][], c2[][], c3[][], C1w[][], C2w[][], C2s[][], C3s[][];
   private double x[][], y[][];
   private int mx,my;
   
   // access mesthods
   public void setGrid( int co[] ){
      mx = co[0];  my = co[1];
   }
   public void setXYarray( double x[][], double y[][] ){
      this.x = x;   this.y = y;
   }
   // method for calculating metrics 
   public void metric(){
      int    i, j, ip, jp, im, jm, mxm=mx-1, mym=my-1;
      double jacob, Xxi, Yxi, Xeta, Yeta;
      XIxW =  new double[mx][mym]; XIyW = new double[mx][mym];
      ETAxW = new double[mx][mym]; JW =   new double[mx][mym]; 
      XIyS =  new double[mxm][my]; JS =   new double[mxm][my];
      ETAxS = new double[mxm][my]; ETAyS = new double[mxm][my];
      ETAyW = new double[mx][mym];
      D1u = new double[mx][mym];  D2u = new double[mx][mym]; 
      D3u = new double[mx][mym];  D4u = new double[mx][mym]; 
      D5u = new double[mx][mym];  D6u = new double[mx][mym];
      D7u = new double[mx][mym];  D8u = new double[mx][mym]; 
      D0u = new double[mx][mym];  D1v = new double[mxm][my]; 
      D2v = new double[mxm][my];  D3v = new double[mxm][my];
      D4v = new double[mxm][my];  D5v = new double[mxm][my]; 
      D6v = new double[mxm][my];  D7v = new double[mxm][my];
      D8v = new double[mxm][my];  D0v = new double[mxm][my];
      D1P = new double[mxm][mym]; D2P = new double[mxm][mym];
      D3P = new double[mxm][mym]; D4P = new double[mxm][mym];
      D5P = new double[mxm][mym]; D6P = new double[mxm][mym];
      D7P = new double[mxm][mym]; D8P = new double[mxm][mym];
      D0P = new double[mxm][mym];
	  
      c1 = new double[mx][my];  c2 = new double[mx][my];  
      c3 = new double[mx][my];
      C1w = new double[mx][mym]; C2w = new double[mx][mym];
      C2s = new double[mxm][my]; C3s = new double[mxm][my];
	  
      double XIx[][] = new double[mx][my], XIy[][] = new double[mx][my];
      double ETAx[][] = new double[mx][my], ETAy[][] = new double[mx][my];
      double J[][] = new double [mx][my];  // double c1[][]=new double [mx][my];
      // double c2[][]=new double[mx][my];    double c3[][]=new double [mx][my];
      double c1w, c1e, c2w, c2e, c2s, c2n, c3s, c3n;
      int    k, k1;
	  
      // g.drawString("　　Calculation start of metrics", 30, 95);

      // Calculate metrixes at nodal points
      for(i = 0; i<mx; i++){
         for(j = 0; j<my; j++){
            if(i == 0){
               Xxi = 0.5*(-3.0*x[0][j]+4.0*x[1][j]-x[2][j]);
               Yxi = 0.5*(-3.0*y[0][j]+4.0*y[1][j]-y[2][j]);
            }else if(i == mxm){
               Xxi = 0.5*(x[mx-3][j]-4.0*x[mx-2][j]+3.0*x[mxm][j]);
               Yxi = 0.5*(y[mx-3][j]-4.0*y[mx-2][j]+3.0*y[mxm][j]);
            }else{
               Xxi = -0.5*(x[i-1][j]-x[i+1][j]);
               Yxi = -0.5*(y[i-1][j]-y[i+1][j]);
            }
            if(j == 0){
               Xeta=0.5*(-3.0*x[i][0]+4.0*x[i][1]-x[i][2]);
               Yeta=0.5*(-3.0*y[i][0]+4.0*y[i][1]-y[i][2]);
            }else if(j == mym){
               Xeta=0.5*(x[i][my-3]-4.0*x[i][my-2]+3.0*x[i][mym]);
               Yeta=0.5*(y[i][my-3]-4.0*y[i][my-2]+3.0*y[i][mym]);
            }else{
               Xeta=-0.5*(x[i][j-1]-x[i][j+1]);
               Yeta=-0.5*(y[i][j-1]-y[i][j+1]);
            }
            jacob=Xxi*Yeta-Xeta*Yxi;   J[i][j]=jacob; 
            XIx[i][j] = Yeta/jacob;    ETAx[i][j]=-Yxi/jacob;
            XIy[i][j] =-Xeta/jacob;    ETAy[i][j]= Xxi/jacob;
         }
      }
      //  Calculate corresponding value for contravariant velocitie U 
      for(i = 0; i<mx; i++){
         for(j = 0; j<mym; j++){  jp=j+1;
            XIxW[i][j]= 0.5*(XIx[i][j]+XIx[i][jp]);
            ETAxW[i][j]=0.5*(ETAx[i][j]+ETAx[i][jp]);
            XIyW[i][j]= 0.5*(XIy[i][j]+XIy[i][jp]);
            JW[i][j] =  0.5*(J[i][j]+J[i][jp]);
            ETAyW[i][j]=0.5*(ETAy[i][j]+ETAy[i][jp]);
         }
      }	  
      //                          value for contravariant velocitie V 
      for(i = 0; i<mxm; i++){    ip=i+1;
         for(j = 0; j<my; j++){ 
            ETAxS[i][j]=0.5*(ETAx[i][j]+ETAx[ip][j]);
            XIyS[i][j]= 0.5*(XIy[i][j]+XIy[ip][j]);
            ETAyS[i][j]=0.5*(ETAy[i][j]+ETAy[ip][j]);
            JS[i][j] =  0.5*(J[i][j]+J[ip][j]);
         }
      } 
      // g.drawString("　　　　Laplacian coefficients", 30, 110);

      // Calculate coefficients of diffusion terms
      for(i = 0; i<mx; i++){
         for(j = 0; j<my; j++){
            c1[i][j]=J[i][j]*(XIx[i][j]*XIx[i][j]+XIy[i][j]*XIy[i][j]);
            c3[i][j]=J[i][j]*(ETAx[i][j]*ETAx[i][j]+ETAy[i][j]*ETAy[i][j]);
            c2[i][j]=J[i][j]*(XIx[i][j]*ETAx[i][j]+XIy[i][j]*ETAy[i][j]);
         }
      }
      for(i = 0; i<mx; i++){
         for(j = 0; j<mym; j++){  jp=j+1;
            C1w[i][j]= 0.5*(c1[i][j]+c1[i][jp]);
            C2w[i][j]= 0.5*(c2[i][j]+c2[i][jp]);
         }
      }	  
      //                          value for contravariant velocitie V 
      for(i = 0; i<mxm; i++){    ip=i+1;
         for(j = 0; j<my; j++){ 
            C2s[i][j]= 0.5*(c2[i][j]+c2[ip][j]);
            C3s[i][j]= 0.5*(c3[i][j]+c3[ip][j]);
         }
     }
	  
     //    coefficients for velocity u
     for(i = 1; i<mxm; i++){     im=i-1;  ip=i+1;  
        for(j = 0; j<mym; j++){  jp=j+1;
           c1w=0.25*(c1[im][j]+c1[im][jp]+c1[i][j]+c1[i][jp]);  
           c2w=0.25*(c2[im][j]+c2[im][jp]+c2[i][j]+c2[i][jp]);
           c1e=0.25*(c1[i][j]+c1[i][jp]+c1[ip][j]+c1[ip][jp]);
           c2e=0.25*(c2[i][j]+c2[i][jp]+c2[ip][j]+c2[ip][jp]);
           c2s=c2[i][j];  c3s=c3[i][j]; c2n=c2[i][jp];  c3n=c3[i][jp];
           D1u[i][j]=0.25*(c2w+c2s);    D2u[i][j]=c1w-0.25*(c2n-c2s);
           D3u[i][j]=-0.25*(c2w+c2n);    D4u[i][j]=c3s-0.25*(c2e-c2w); 
           D5u[i][j]=c3n+0.25*(c2e-c2w); D6u[i][j]=-0.25*(c2e+c2s);
           D7u[i][j]=c1e+0.25*(c2n-c2s); D8u[i][j]=0.25*(c2e+c2n);
           D0u[i][j]=-c1e-c1w-c3n-c3s; 
        }
     }
     //    coefficients for velocity v
     for(i = 0; i<mxm; i++){     ip=i+1;    
        for(j = 0; j<mym; j++){  jp=j+1; jm=j-1;
           c1w=c1[i][j]; c1e=c1[ip][j]; c2w=c2[i][j]; c2e=c2[ip][j];
           if(j == 0){ k=mxm-i;  k1=k-1;
              c2s=0.25*(c2[i][0]+c2[k][1]+c2[ip][0]+c2[k1][1]);
              c3s=0.25*(c3[i][0]+c3[k][1]+c3[ip][0]+c3[k1][1]);
           }else{ 						  
              c2s=0.25*(c2[i][jm]+c2[i][j]+c2[ip][jm]+c2[ip][j]);
              c3s=0.25*(c3[i][jm]+c3[i][j]+c3[ip][jm]+c3[ip][j]);
           }
           c2n=0.25*(c2[i][j]+c2[i][jp]+c2[ip][j]+c2[ip][jp]);
           c3n=0.25*(c3[i][j]+c3[i][jp]+c3[ip][j]+c3[ip][jp]);
           D1v[i][j]=0.25*(c2w+c2s);      D2v[i][j]=c1w-0.25*(c2n-c2s);
           D3v[i][j]=-0.25*(c2w+c2n);     D4v[i][j]=c3s-0.25*(c2e-c2w);
           D5v[i][j]=c3n+(c2e-c2w);       D6v[i][j]=-0.25*(c2e+c2s); 
           D7v[i][j]=c1e+0.25*(c2n-c2s);  D8v[i][j]=0.25*(c2e+c2n);
           D0v[i][j]=-c1e-c1w-c3n-c3s;  
        }
      }

      //  coefficients for solving Poisson equation of pressure
      for(i = 0; i<mxm; i++){     ip=i+1;         
         for(j = 0; j<mym; j++){  jp=j+1;
            c1w=0.5*(c1[i][j]+c1[i][jp]); c1e=0.5*(c1[ip][j]+c1[ip][jp]);
            c2w=0.5*(c2[i][j]+c2[i][jp]); c2e=0.5*(c2[ip][j]+c2[ip][jp]);
            c2s=0.5*(c2[i][j]+c2[ip][j]); c2n=0.5*(c2[i][jp]+c2[ip][jp]);
            c3s=0.5*(c3[i][j]+c3[ip][j]); c3n=0.5*(c3[i][jp]+c3[ip][jp]);
            D1P[i][j]=0.25*(c2w+c2s);     D2P[i][j]=c1w-0.25*(c2n-c2s);
            D3P[i][j]=-0.25*(c2w+c2n);    D4P[i][j]=c3s-0.25*(c2e-c2w); 
            D5P[i][j]=c3n+0.25*(c2e-c2w); D6P[i][j]=-0.25*(c2e+c2s);
            D7P[i][j]=c1e+0.25*(c2n-c2s); D8P[i][j]=0.25*(c2e+c2n);
            D0P[i][j]=-c1e-c1w-c3n-c3s;
         }
      } 
      // g.drawString("END of Metrix calculation", 30, 185);
   }
}
