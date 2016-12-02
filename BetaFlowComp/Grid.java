/*********************************************************************/
/*   Input data for simulating incompressible fluid flow             */
/*        GRID class:  Mesh generation and its display               */
/*                              Last update November 9,2001.         */
/*         All rights reserved, Copyright (C) 2001, Kiyoshi Minemura */ 
/*********************************************************************/
import java.awt.*;

//  ===  Class of grid generation  ===================================
public class Grid {  
   private int   mx, my;       //  number of mesh points
   public double x[][], y[][]; //  nodal coordinates
   private float xp[], yp[];   //  one dimensional nodal coordinates
   public int    mxRin, mxRout, mxhalf, mxL1, mxL2, mxU1, mxU2; 
   private int    E[][];       //  node data of triangular mesh 
   private float xrange [] = new float [2];
   private float yrange [] = new float [2];
   private String flow;
   
   //    Method for selecting flow object
   public void grid_select(String flow){
      this.flow = flow;
      if(flow == "duct")      grid_duct();
      else if(flow == "bend") grid_bend();
      else if(flow == "cylinder"){
                              grid_cylinder();   grid_E();  grid_EP();
      }
      oneDim(); 
   }
   
   //    EXAMPLE No. 1 (flow = "duct")
   //   ====== Mesh for calculating a straight pipe flow =========
   private void grid_duct(){
      mx=17;  my=17;           // node numbers specified
      x = new double [mx][my]; y = new double [mx][my];
      double a, co, b, sj, s0, s1, s2; 
      int i, j;
	  
      // perpendicular direction ( duct width= 1)
      s0=2.0/(double)(my-1);
      a=1.12;  b=(a+1.0)/(a-1.0);  co=1.0/(b-1.0);
      for(j=0; j<my; j++){
         sj=s0*(double)j;               // equal ratio
         s1=Math.pow(b, sj); s2=Math.pow(b, sj-1.0);
         y[0][j]=co*(s1-1.0)/(s2+1.0); x[0][j]=0.0;
      }
      // axial direction  ( duct length=10) 
      a=1.6;  b=(a+1.0)/(a-1.0);   co=10.0/(b-1.0);
      for(i=0; i<mx; i++){
         sj=(double)i/(double)(mx-1);   // equal ratio
         s1=Math.pow(b, sj); 
         x[i][0]=co*(s1-1.0);
      }  
      //  
      for(i=1 ; i<mx; i++){
         for(j=0; j<my; j++){ x[i][j]=x[i][0]; y[i][j]=y[0][j];}
      }
      //  graphic output range
      xrange[0]=1.5f; xrange[1]=9.5f; yrange[0]=-4.0f; yrange[1]=4.0f;
   }
 
   //  Example No. 2 (flow ="bend")    
   // ====@Mesh for calculating a bend flow =======================
   private void grid_bend(){
      mx=33;  my=21;  
      x = new double [mx][my];    y = new double [mx][my];
      double r2el=1.5;    // dimensionless pipe radius (center)
      double L1=Math.PI*r2el;  // length of straight pipes before and after bend
      double a, b, co, xi, sj, s1, s2, theta, pi=Math.PI; 
      int    i, j;

      // inlet part of stright pipe
      a=1.1;  b=(a+1.0)/(a-1.0);  co=1.0/(b-1.0); 
      mxRin=35*(mx-1)/100;    // node numbers of inlet pipe=35% of mx        
      for(j=0 ; j<my; j++){
         sj=2.0*(double)j/(double)(my-1);                 // equal ratio
         s1=Math.pow(b, sj); s2=Math.pow(b, sj-1.0);
         y[0][j]=-r2el+co*(s1-1.0)/(s2+1.0);  x[0][j]=0.0;	
      }
      a=1.18;  b=(a+1.0)/(a-1.0);  co=L1/(b-1.0);
      for(i=1; i<=mxRin; i++){
         sj=2.0*(double)i/(double)mxRin;   // equal ratio
         s1=Math.pow(b, sj);   s2=Math.pow(b, sj-1.0); 
         xi=co*(s1-1.0)/(s2+1.0); 
         for(j=0; j<my; j++){ 
            x[i][j]=xi;    y[i][j]=y[0][j]; 
         }		 
      }
  
      //    bend part
      mxRout=mxRin+3*(mx-1)/10; // node numbers of bend = 30% of mx 
      co=0.5*pi/(double)(mxRout-mxRin);  
      for(i=mxRin+1; i<mxRout; i++){
         for(j=0; j<my; j++){
            theta=co*(double)(i-mxRin);
            x[i][j]=-y[0][j]*Math.sin(theta)+L1;
            y[i][j]=y[0][j]*Math.cos(theta);
         }
      }

      //   outlet part of stright pipe
      a=1.1;   b=(a+1.0)/(a-1.0);  co=1.0/(b-1.0); 
      for(j=0; j<my; j++){
         sj=2.0*(double)j/(double)(my-1);   // equal ratio
         s1=Math.pow(b, sj);  s2=Math.pow(b, sj-1.0);
         x[mxRout][j]=L1+r2el-co*(s1-1.0)/(s2+1.0);
      }
      co=L1/Math.pow((double)(mx-mxRout-1), 1.18);
      for(i=mxRout+1; i<mx; i++){
         xi=co*Math.pow((double)(i-mxRout), 1.18);
         for(j=0; j<my; j++){
            y[i][j]=xi;  x[i][j]=x[mxRout][j]; 
         }
      }  
      // graphic output range
      xrange[0]=0.0f;              yrange[0]=-(float)r2el; 
      xrange[1]=(float)(L1+r2el);  yrange[1]=(float)L1;   
   }

   // Example No.3 (flow = "cylinder")
   // ========== C-type  "Flow around cylindar with slip boundary"
   private void grid_cylinder(){
      mx=121; my=24;  // mx=grid number in x direction, must be odd number. 
      //mx=99, my=18;  mx=75, my=17; mx=71; my=12;
      x = new double [mx][my];    y = new double [mx][my];
      double a, b, co, sj, s0, s1, s2, pi=Math.PI; 
      int i, j, k;   
      a=1.12;  b=(a+1.0)/(a-1.0);  co=1.0*(a-1.0); 
	   
      double rr=1.0;        //  radisu of Cylinder
      double bhalf=9.0*rr;  //  half width of channel
      double Lds=20.0*rr;   //  distance between cylinder center and outlet
      double Lus=8.0*rr-bhalf*(2.0-Math.sqrt(3.0)); // distance between edge and center 
      // dependent parameters
      mxhalf=mx/2; mxL1=mxhalf-18; mxL2=mxhalf+18; mxU1=mxhalf-5; mxU2=mxhalf+5;
      double ang10=pi/18.0, rateI, rateO, Lx, yy;
   
      //  inlet and outlet stright parts
      rateI=(Lds-rr)/(Math.pow(1.04, (double)(mxL1))-1.0);
      rateO=(Lds-2.0*rr)/(Math.pow(1.02, (double)(mxL1))-1.0);
      for(i=0; i<=mxL1; i++){
         x[i][0]=rr+rateI*(Math.pow(1.04, (double)(mxL1-i))-1.0);
         y[i][0]=0.0;
         x[i][my-1]=2.0*rr+rateO*(Math.pow(1.02, (double)(mxL1-i))-1.0);
         y[i][my-1]=-bhalf;
      }
      // periphery of cylinder
      rateO=(2.0*rr+Lus)/(double)(mxU1-mxL1);
      for(i=mxL1+1; i<=mxU1; i++){
         x[i][0]= rr*Math.cos(ang10*(double)(i-mxL1)); 
         y[i][0]=-rr*Math.sin(ang10*(double)(i-mxL1));
         x[i][my-1]=2.0*rr-rateO*(double)(i-mxL1);
         y[i][my-1]=-bhalf;
      }
      // upstream periphery of cyliner
      rateO=bhalf/(Math.pow(1.15, (double)(mxhalf-mxU1))-1.0);
      for(i=mxU1+1; i<=mxhalf; i++){
         x[i][0]= rr*Math.cos(ang10*(double)(i-mxL1));
         y[i][0]=-rr*Math.sin(ang10*(double)(i-mxL1));
         yy=-bhalf+rateO*(Math.pow(1.15, (double)(i)-mxU1)-1.0);
         Lx=Math.sqrt(4.0-yy*yy/(bhalf*bhalf));
         x[i][my-1]=-(bhalf*(Lx-Math.sqrt(3.0))+Lus);
         y[i][my-1]=yy;
      }
      // left side of channel
      for(i=mxhalf+1; i<mx; i++){
         j=i-2*(i-mxhalf);
         x[i][0]=x[j][0];       y[i][0]=-y[j][0];
         x[i][my-1]=x[j][my-1]; y[i][my-1]=-y[j][my-1];
      }
      //  inlet and outlet boundary
      a=1.5;  b=(a+1.0)/(a-1.0);  co=1.0/(b-1.0); 
      for(j=0; j<my; j++){
         x[0][j]=Lds+rr;
         sj=(double)j/(double)(my-1);     // equal ratio
         s1=Math.pow(b, sj);         // s2=Math.pow(b, sj-1.0);
         y[0][j]=-bhalf*co*(s1-1.0); // -bhalf*aa*(s1-1.0)/(s2+1.0);
         x[mx-1][j]=x[0][j];   y[mx-1][j]=-y[0][j];
      }
	  
      // internal mesh points
      k=(int)(0.5*pi/ang10)+mxL1;
      for(i=1; i<=k; i++){  
         for(j=1; j<=my-2; j++){
            x[i][j]=0.9*Lds*(double)(k-i)/(double)k; 
            y[i][j]=-1.8*rr*(1.0+1.8*rr*(double)j/(double)(my-1));
         }
      }
      for(i=k+1; i<=mxhalf; i++){
         for(j=1; j<=my-2; j++){
            x[i][j]=y[k][j]*Math.sin(ang10*(i-k));
            y[i][j]=y[k][j]*Math.cos(ang10*(i-k));
         }
      }
      for(i=mxhalf+1; i<=mx-2; i++){
         k=2*mxhalf-i;
         for(j=1; j<=my-2; j++){
            x[i][j]= x[k][j];
            y[i][j]=-y[k][j];
         }
      } 
      xrange[0]=-(float)(rr*2.0);  xrange[1]=(float)(Lds*0.85);  
      yrange[0]=-(float)(0.9*(Lds+rr)*0.5); yrange[1]=-(float)yrange[0];   
   }  // End of grid_C
   
   // Grid generation using Laplace equation
   private void grid_E(){
      double P[][] = new double[mx][my];
      double Q[][] = new double[mx][my];
      int i, j, k=0, kmax=500, k3=kmax/10, k6=k3*2;
      double eps1=5.0E-3, resid=0.0, rate;
      double Xxi, Yxi, Xeta, Yeta, alph, beta, gamma, jac, ag, pp, qq, sx, sy;
	  
      while(k < kmax){
         resid=0.0;
         if(k <= k3) rate=0.25; else rate=1.0;
         for(i=1; i<mx-1; i++){
            for(j=1; j<my-1; j++){
               Xxi=(x[i+1][j]-x[i-1][j])*0.5;
               Yxi=(y[i+1][j]-y[i-1][j])*0.5;
               Xeta=(x[i][j+1]-x[i][j-1])*0.5;
               Yeta=(y[i][j+1]-y[i][j-1])*0.5;
               alph=Xeta*Xeta+Yeta*Yeta;
               beta=Xxi*Xeta+Yxi*Yeta;
               gamma=Xxi*Xxi+Yxi*Yxi;
               jac=Xxi*Yeta-Xeta*Yxi;
               ag=0.5/(alph+gamma);   pp=P[i][j]; qq=Q[i][j];
               sx=(alph*(x[i+1][j]+x[i-1][j])
                  -0.5*beta*(x[i+1][j+1]-x[i-1][j+1]-x[i+1][j-1]+x[i-1][j-1])
                  +gamma*(x[i][j+1]+x[i][j-1])
                  +jac*jac*(pp*Xxi+qq*Xeta))*ag-x[i][j];
               sy=(alph*(y[i+1][j]+y[i-1][j])
                  -0.5*beta*(y[i+1][j+1]-y[i-1][j+1]-y[i+1][j-1]+y[i-1][j-1])
                  +gamma*(y[i][j+1]+y[i][j-1])
                  +jac*jac*(pp*Yxi+qq*Yeta))*ag-y[i][j];
               x[i][j]=x[i][j]+rate*sx;
               y[i][j]=y[i][j]+rate*sy;
               resid=resid+sx*sx+sy*sy;
            }
         }
         k++;
         if(resid < eps1) break;
      }
      //  g.drawString("  k="+k, 10,10); g.drawString("  residual = "+resid,10,25);
   }

   // Grid generation using Poisson equation
   private void grid_EP(){
      int i, j;
      double P[][] = new double[mx][my];
      double Q[][] = new double[mx][my];
      int    kmax=1000, k=0, k1=kmax/10;
      double eps1=1.0E-4, resid=0.0, rate;
      double Xxi, Yxi, Xeta, Yeta, alph, beta, gamma, jac, ag, pp, qq, sx, sy;
	  
      for(i=0; i<mx; i++){
         for(j=0; j<my-4; j++){
            if(i<=mxU1 || i>= mxU2)
               Q[i][j]=-5.0*Math.exp(-0.3*j);//
               // Q[i][j]=25*Math.exp(-0.7*(my-j));
               // Q[i][j]=-45.0*Math.exp(-0.3*j);
         }
      }

      while(k < kmax){
         resid=0.0;
         if(k <= k1) rate=0.2; else rate=1.0;  
         for(i=1; i<mx-1; i++){
            for(j=1; j<my-1; j++){
               Xxi=(x[i+1][j]-x[i-1][j])*0.5;
               Yxi=(y[i+1][j]-y[i-1][j])*0.5;
               Xeta=(x[i][j+1]-x[i][j-1])*0.5;
               Yeta=(y[i][j+1]-y[i][j-1])*0.5;
               alph=Xeta*Xeta+Yeta*Yeta;
               beta=Xxi*Xeta+Yxi*Yeta;
               gamma=Xxi*Xxi+Yxi*Yxi;
               jac=Xxi*Yeta-Xeta*Yxi;
               ag=0.5/(alph+gamma);   
               pp=0.0;  qq=Q[i][j];
               sx=(alph*(x[i+1][j]+x[i-1][j])
                  -0.5*beta*(x[i+1][j+1]-x[i-1][j+1]-x[i+1][j-1]+x[i-1][j-1])
                  +gamma*(x[i][j+1]+x[i][j-1])
                  +jac*jac*(pp*Xxi+qq*Xeta))*ag-x[i][j];
               sy=(alph*(y[i+1][j]+y[i-1][j])
                  -0.5*beta*(y[i+1][j+1]-y[i-1][j+1]-y[i+1][j-1]+y[i-1][j-1])
                  +gamma*(y[i][j+1]+y[i][j-1])
                  +jac*jac*(pp*Yxi+qq*Yeta))*ag-y[i][j];
               x[i][j]=x[i][j]+rate*sx;
               y[i][j]=y[i][j]+rate*sy;
               resid=resid+sx*sx+sy*sy;
            }
         }
         k++;
         if(resid < eps1) break;
      }
      //  g.drawString("  k="+k, 10,40); 
      //  g.drawString("  residual = "+resid,10,55);
   }   
   
   // accessor methods
   private void oneDim(){
      int k, kmax=mx*my;
      xp=new float [kmax]; yp=new float [kmax];
	  
      for(int i=0; i<mx; i++){
         for(int j=0; j<my; j++){
            k=my*i+j; xp[k]=(float)x[i][j]; yp[k]=(float)y[i][j];
         }
      }
   }
   public float [] getXrange(){ return xrange; }
   public float [] getYrange(){ return yrange; }
   public float [] getXarray(){ return xp; }
   public float [] getYarray(){ return yp; }
   public int getMx(){  return mx; }
   public int getMy(){  return my; }
   public String getFlow(){ return flow; }
}
