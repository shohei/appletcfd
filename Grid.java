/**====================================================================
   Java Applet for generating Mesh
        Grid.java   ( Grid generation )
             All Rights Reserved, Copyright (C) 2001-2003, K. MINEMURA
       　　　　　　 Last updated by K. Minemura on May 24, 2003.
======================================================================*/
import java.awt.*;

//  ===  Grid Class for Grid Generation   ===================
public class Grid {
   private int  mx, my;      //  number of mesh
   private double x[][], y[][], Q[][], P[][],xx[],yy[];
   private final double PI = Math.PI,  eps = 1.0e-4;
   private double xmin, ymin, xmax, ymax;
   private int    mxRin, mxRout, mxhalf, mxL1, mxL2, mxU1, mxU2;
   private int    E[][];   //  node data of triangular mesh
   private double resid;   //  residual
   private int    iter;    //  iteration numbers

   // EXAMPLE No. 1 (flow = "Duct")
   //   ====== Mesh for calculating a straight pipe flow =========
   public void gridDuct(){
      mx = 17;    my = 17;
      x = new double [mx][my];    y = new double [mx][my];

      double s0 = 2.0/(double)(my-1), sj, s1, s2;
      double a0 = 1.12, aa = 0.5*(a0-1.0), b0 = (a0+1.0)/(a0-1.0);
      for(int j=0; j<my; j++){
         sj = s0*(double)j;               // equal ratio
         s1 = Math.pow(b0, sj);   s2 = Math.pow(b0, sj-1.0);
         y[0][j] = aa*(s1-1.0)/(s2+1.0);
         x[0][j] = 0.0;
      }
      s0 = 10.0/(double)(mx-1);
      for(int i=1 ; i<mx; i++){
         sj = s0*(double)i;               // equal distance
         for(int j=0; j<my; j++){ x[i][j] = sj; y[i][j] = y[0][j]; }
      }
      // range for displaying the physical space
      xmin = 0.0;  ymin = -3.0;  xmax = 10.0;  ymax = 3.0;
   }

   //  Example No. 2 (flow ="Elbow")
   // ====　Mesh for calculating a elbow flow =======================
   public void gridElbow(){
      mx = 31;   my = 17;
      x = new double [mx][my]; y = new double [mx][my];

      double s0 = 2.0/(double)(my-1), sj, s1, s2;
      double a0 = 1.12, aa = 0.5*(a0-1.0), b0 = (a0+1.0)/(a0-1.0);
      double rate, L1, theta, xy;
      double r2el = 1.5;            //  dimensionless pipe radius

      // inlet part of straight pipe
      mxRin = 35*(mx-1)/100;    L1 = PI*r2el;
      for(int j=0 ; j<my; j++){
         sj = s0*(double)j;                 // equal ratio
         s1 = Math.pow(b0, sj);   s2 = Math.pow(b0, sj-1.0);
         y[0][j] = -r2el+aa*(s1-1.0)/(s2+1.0);
         x[0][j] = 0.0;
      }
      rate = L1/(double)mxRin;
      for(int i=1 ; i<=mxRin; i++){
         xy = rate*(double)i;
         for(int j=0; j<my; j++){ x[i][j]=xy; y[i][j]=y[0][j];}
      }

      //    bend part
      mxRout = mxRin+3*(mx-1)/10;
      rate = 0.5*PI/(double)(mxRout-mxRin);
      for(int i=mxRin+1; i<mxRout; i++){
         for(int j=0; j<my; j++){
            theta = rate*(double)(i-mxRin);
            x[i][j] = -y[0][j]*Math.sin(theta)+L1;
            y[i][j] = y[0][j]*Math.cos(theta);
         }
      }

      //   outlet part of straight pipe
      for(int j=0; j<my; j++){
         sj=s0*(double)j;                 // equal ratio
         s1=Math.pow(b0, sj); s2=Math.pow(b0, sj-1.0);
         x[mxRout][j]=L1+r2el-aa*(s1-1.0)/(s2+1.0);
      }
      rate=L1/Math.pow((double)(mx-mxRout-1), 1.18);
      for(int i=mxRout+1; i<mx; i++){
         xy=rate*Math.pow((double)(i-mxRout), 1.18);
         for(int j=0; j<my; j++){
            y[i][j]=xy;  x[i][j]=x[mxRout][j];
         }
      }

      // range for displaying the physical space
      xmin = 0.0; ymin = -r2el; xmax = L1+r2el;  ymax = 1.1*L1;
   }

   // ==== O-type  "Flow around cylinder in infinitely large space =========
   public void gridO(){
      mx = 50;    my = 25;
      x = new double [mx][my];    y = new double [mx][my];

      double r1 = 0.5,  r2 = 20.0;
      double a0 = 1.12, aa = 1.0*(a0-1.0), b0 = (a0+1.0)/(a0-1.0);
      double s0 = 1.0/((double)my-1.0), sj;
      aa = (r2-r1)/(b0-1.0);
      double theta, r, rate = 2.0*PI/(double)(mx-1);

      for(int j=0; j<my; j++){
         sj = s0*(double)j;                 // equal ratio
         r = r1+aa*(Math.pow(b0, sj)-1.0);
         for(int i=0; i<mx; i++){
            theta = -PI+rate*(double)(i);
            x[i][j] = r*Math.cos(theta);
            y[i][j] = r*Math.sin(theta);
         }
      }
      // range for displaying the physical space
      xmin = -r2; ymin = -r2; xmax = r2; ymax = r2*1.1;
   }

   // Example No.3 (flow = "cylinder")
   // ========== C-type  "Flow around cylinder with slip boundary" =====
   public void gridCylinder(){
      mx = 121;    my = 24;  // mx=grid number in x direction, must be odd number.
      //  mx=73, my=17;  mx=71; my=12; mx=121, my=24;
      x = new double [mx][my];    y = new double [mx][my];
      P = new double [mx][my];    Q = new double [mx][my];

      int i, j, k;
      double a0 = 1.12, aa = 1.0*(a0-1.0);
      double s0, s1, s2;
      double rr = 1.0;        //  radius of Cylinder
      double bhalf = 9.0*rr;  //  half width of channel
      double Lds = 21.0*rr;   //  distance between cylinder center and outlet
      double Lus = 8.0*rr-bhalf*(2.0-Math.sqrt(3.0)); // distance between edge and center
      // dependent parameters
      mxhalf = mx/2; mxL1 = mxhalf-18; mxL2 = mxhalf+18;
      mxU1 = mxhalf-5; mxU2 = mxhalf+5;
      double ang10 = PI/18.0, rateI, rateO, Lx, yy0, sj, ix;

      //  inlet and outlet straight parts
      double a=1.5,  b=(a+1.0)/(a-1.0), b0;
                a=2.7;  b0=(a+1.0)/(a-1.0);
      rateI = (Lds-rr)/(b-1.0);// (Math.pow(1.02, (double)(mxL1))-1.0);
      rateO = (Lds-1.5*rr)/(b0-1.0);// (Lds-2.0*rr)/(Math.pow(1.02, (double)(mxL1))-1.0);
      for(i=0; i<=mxL1; i++){
		  ix=(double)(mxL1-i)/(double)(mxL1);
         x[i][0] = rr+rateI*(Math.pow( b, ix )-1.0);// (double)(mxL1-i))-1.0);
         y[i][0] = 0.0;
         x[i][my-1] = 1.5*rr+rateO*(Math.pow(b0,ix)-1.0);//(Math.pow(1.02, (double)(mxL1-i))-1.0);
         y[i][my-1] = -bhalf;
      }
      // periphery of cylinder
      rateO = (1.5*rr+Lus)/(double)(mxU1-mxL1);
      for(i=mxL1+1; i<=mxU1; i++){
         x[i][0] =  rr*Math.cos(ang10*(double)(i-mxL1));
         y[i][0] = -rr*Math.sin(ang10*(double)(i-mxL1));
         x[i][my-1] = 1.5*rr-rateO*(double)(i-mxL1);
         y[i][my-1] = -bhalf;
      }
      // upstream periphery of cylinder
      rateO = bhalf/(Math.pow(1.15, (double)(mxhalf-mxU1))-1.0);
      for(i=mxU1+1; i<=mxhalf; i++){
         x[i][0] = rr*Math.cos(ang10*(double)(i-mxL1));
         y[i][0] =-rr*Math.sin(ang10*(double)(i-mxL1));
         yy0 = -bhalf+rateO*(Math.pow(1.15, (double)(i)-mxU1)-1.0);
         Lx = Math.sqrt(4.0-yy0*yy0/(bhalf*bhalf));
         x[i][my-1] = -(bhalf*(Lx-Math.sqrt(3.0))+Lus);
         y[i][my-1] = yy0;
      }
      // left side of channel
      for(i=mxhalf+1; i<mx; i++){
         j = i-2*(i-mxhalf);
         x[i][0] = x[j][0];       y[i][0] = -y[j][0];
         x[i][my-1] = x[j][my-1]; y[i][my-1] = -y[j][my-1];
      }
      //  inlet and outlet boundary
      a0 = 6.0;              // factors of mesh distance
      b0 = (a0+1.0)/(a0-1.0); s0 = 1.0/(double)(my-1); // s0=2.0/(double)(my-1);
      aa = 1.0/(b0-1);
      for(j=0; j<my; j++){
         x[0][j] = Lds;
         sj = s0*(double)j;     // equal ratio
         s1 = Math.pow(b0, sj);        // s2=Math.pow(b0, sj-1.0);
         y[0][j] = -bhalf*aa*(s1-1.0); // -bhalf*aa*(s1-1.0)/(s2+1.0);
         x[mx-1][j] = x[0][j];   y[mx-1][j] = -y[0][j];
      }

      // internal mesh points
      k = (int)(0.5*PI/ang10)+mxL1;
      for(i=1; i<=k; i++){
         for(j=1; j<=my-2; j++){
            x[i][j] = 0.9*Lds*(double)(k-i)/(double)k;
            y[i][j] = -1.8*rr*(1.0+1.8*rr*(double)j/(double)(my-1));
         }
      }
      for(i=k+1; i<=mxhalf; i++){
         for(j=1; j<=my-2; j++){
            x[i][j] = y[k][j]*Math.sin(ang10*(i-k));
            y[i][j] = y[k][j]*Math.cos(ang10*(i-k));
         }
      }
      for(i=mxhalf+1; i<=mx-2; i++){
         k = 2*mxhalf-i;
         for(j=1; j<=my-2; j++){
            x[i][j] = x[k][j];
            y[i][j] =-y[k][j];
         }
      }
	  // range for displaying the physical space
      xmin = -6*rr;  xmax = Lds; ymin = -(Lds+6*rr)*0.5; ymax = -ymin;
   }

   // Nodal points decided by Laplace equation
   public void gridE(){
      int i, j, k=0, kmax=500, k3=kmax/10, k6=k3*2;
      double eps1 = 50.0*eps, rate;
      double Xxi, Yxi, Xeta, Yeta, alph, beta, gamma, jac, ag, pp, qq, sx, sy;

      while(k < kmax){
         resid = 0.0;
         if(k <= k3) rate = 0.25; else rate = 1.0;
         for(i=1; i<mx-1; i++){
            for(j=1; j<my-1; j++){
               Xxi = (x[i+1][j]-x[i-1][j])*0.5;
               Yxi = (y[i+1][j]-y[i-1][j])*0.5;
               Xeta = (x[i][j+1]-x[i][j-1])*0.5;
               Yeta = (y[i][j+1]-y[i][j-1])*0.5;
               alph = Xeta*Xeta+Yeta*Yeta;
               beta = Xxi*Xeta+Yxi*Yeta;
               gamma = Xxi*Xxi+Yxi*Yxi;
               jac = Xxi*Yeta-Xeta*Yxi;
               ag = 0.5/(alph+gamma);   pp = P[i][j]; qq = Q[i][j];
               sx = (alph*(x[i+1][j]+x[i-1][j])
                   -0.5*beta*(x[i+1][j+1]-x[i-1][j+1]-x[i+1][j-1]+x[i-1][j-1])
                   +gamma*(x[i][j+1]+x[i][j-1])
                   +jac*jac*(pp*Xxi+qq*Xeta))*ag-x[i][j];
               sy = (alph*(y[i+1][j]+y[i-1][j])
                   -0.5*beta*(y[i+1][j+1]-y[i-1][j+1]-y[i+1][j-1]+y[i-1][j-1])
                   +gamma*(y[i][j+1]+y[i][j-1])
                   +jac*jac*(pp*Yxi+qq*Yeta))*ag-y[i][j];
               x[i][j] = x[i][j]+rate*sx;
               y[i][j] = y[i][j]+rate*sy;
               resid = resid+sx*sx+sy*sy;
            }
         }
         k++;
         if(resid < eps1) break;
      }
      iter = k;
   }

   // Nodal points decided by Poisson equation
   public void gridEP(){
      int    k = 0, kmax = 1000, k1 = kmax/10;
      double rate;
      double Xxi, Yxi, Xeta, Yeta, alph, beta, gamma, jac, ag, pp, qq, sx, sy;

      for(int i=5; i<mx-5; i++){
         for(int j=0; j<my-4; j++){
            if(i<=mxU1 || i>= mxU2)
               Q[i][j] = -5.0*Math.exp(-0.3*j);
               // Q[i][j] = -5*Math.exp(-0.3*j);
               // Q[i][j] = -15.0*Math.exp(-0.3*j);
         }
      }

      while(k < kmax){
         resid = 0.0;
         if(k <= k1) rate = 0.2; else rate = 1.0;
         for(int i=1; i<mx-1; i++){
            for(int j=1; j<my-1; j++){
               Xxi = (x[i+1][j]-x[i-1][j])*0.5;
               Yxi = (y[i+1][j]-y[i-1][j])*0.5;
               Xeta = (x[i][j+1]-x[i][j-1])*0.5;
               Yeta = (y[i][j+1]-y[i][j-1])*0.5;
               alph = Xeta*Xeta+Yeta*Yeta;
               beta = Xxi*Xeta+Yxi*Yeta;
               gamma = Xxi*Xxi+Yxi*Yxi;
               jac = Xxi*Yeta-Xeta*Yxi;
               ag = 0.5/(alph+gamma);
               pp = 0.0;    qq = Q[i][j];
               sx = (alph*(x[i+1][j]+x[i-1][j])
                   -0.5*beta*(x[i+1][j+1]-x[i-1][j+1]-x[i+1][j-1]+x[i-1][j-1])
                   +gamma*(x[i][j+1]+x[i][j-1])
                   +jac*jac*(pp*Xxi+qq*Xeta))*ag-x[i][j];
               sy = (alph*(y[i+1][j]+y[i-1][j])
                   -0.5*beta*(y[i+1][j+1]-y[i-1][j+1]-y[i+1][j-1]+y[i-1][j-1])
                   +gamma*(y[i][j+1]+y[i][j-1])
                   +jac*jac*(pp*Yxi+qq*Yeta))*ag-y[i][j];
               x[i][j] = x[i][j]+rate*sx;
               y[i][j] = y[i][j]+rate*sy;
               resid = resid+sx*sx+sy*sy;
            }
         }
         k++;
         if(resid < eps) break;
      }
      iter=k;
   }
   // === method for numbering node ====
   public void linkNode(){
      int  nmax = 2*(mx-1)*(my-1);
      E = new int [nmax][3];
      oneDim();
      // number all the sides of triangular meshes
      int k = -1;
      for(int i=0; i<mx-1; i++){
         for(int j=0; j<my-1; j++){
            k++; E[k][0] = my*i+j;   E[k][1] = my*i+j+1;   E[k][2] = my*(i+1)+j;
            k++; E[k][0] = my*i+j+1; E[k][1] = my*(i+1)+j; E[k][2] = my*(i+1)+j+1;
         }
      }

   }
   //  === method for converting two-dimensional data to one dimensional one ===
   public void oneDim(){
      int k, kmax = mx*my;
      xx = new double [kmax];  yy = new double [kmax];
      for(int i=0; i<mx; i++){
         for(int j=0; j<my; j++){
            k = my*i+j;  xx[k] = x[i][j];  yy[k] = y[i][j];
         }
      }
   }


   // ===  access methods for exchanging date with other classes ===
   public double[] getRange(){
      double range[] = {xmin, xmax, ymin, ymax};
      return range;
   }
   public int getMx(){ return mx; }
   public int getMy(){ return my; }
   public double [] getXarray(){ return xx; }
   public double [] getYarray(){ return yy; }
   public double getResid(){     return resid; }
   public int    getIter(){  return iter;  }
   public int[][] getLink(){ return E; }
}
