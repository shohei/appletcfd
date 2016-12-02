/**===================================================================
    Simulation of a Convection Equation
    　 ConvSolver.java  ( Solver part )
             All Rights Reserved, Copyright (C) 2001-2002, K. MINEMURA
       　　　　　　 Last updated by K. Minemura on November 4, 2002.
=====================================================================*/
import java.awt.*;

//     Solver Class for Convection Equation
public class ConvSolver{
   //  Declaration of instance variables related finite difference method
   public static double t;         // elapsed time
   public static double dt=0.02;   // time step for iteration
   public static double Cour;      // Courant number
   public static double diff;      // diffusion number
   private double f0[],ftcs[],fp[],fpc[],fquick[],fkk[]; // dependent variables
   private double fc[],fc0[],g[],g0[],fs[],fa[];
   private double x[];             // grid coordinate
   private int    Nx=100;          // number of grid
   private double u0 = 1.0;        // velocity
   private double xmin, xmax, ymin, ymax; // physical drawing range for applet
   private double dx;              // grid length
   private double alpha = 0.05;    // diffusion coefficient

   // Constructor for initialization of generated variables
   public ConvSolver(){
      initFlow();
   }

   //  ===  Setting up method for initial values ===
   public void initFlow(){
      xmin=0.0; xmax=5.0; ymin=-0.5; ymax=1.5; // range of presentation
      x = new double [Nx];     // grid coordinate of x-direction
      dx=(xmax-xmin)/(double)(Nx-1);
      for(int i=0; i<Nx; i++)   x[i] = dx*(double)i;
      ftcs  = new double [Nx]; // value for FTCS scheme
      f0 = new double [Nx];    // value before one step (FTCS)
      fa = new double [Nx];    // same with f0
      g  = new double [Nx];
      fc = new double [Nx];
      fc0= new double [Nx];
      g0 = new double [Nx];
      fp  = new double [Nx];
      fquick=new double [Nx];
      fkk  = new double [Nx];
      fs = new double [Nx];    // initial value
      fpc= new double [Nx];


      Cour = u0*dt/dx;          // Courant number
      diff = alpha*dt/(dx*dx);  // diffusion number
      t=0.0;                    // initialization of elapsed time
      for(int i=0; i<Nx; i++){  // install initial value
         f0[i] = 1.0;
         if(i < (int)(0.4/dx) || i> (int)(1.2/dx)) f0[i] = 0.0;
      }
         fc0=(double[])f0.clone();   // copy of array
         fa=(double[])f0.clone();    fs=(double[])f0.clone();
   }

   // ===  Calculating method using FTCS scheme ===
   public void solvFTCS(){
      for(int i=1; i<Nx-1; i++)  ftcs[i] = f0[i]-0.5*Cour*(f0[i+1]-f0[i-1])
                                          +diff*(f0[i+1]-2.0*f0[i]+f0[i-1]);
      setBoundary( ftcs );          //  install boundary value
      f0 = (double[]) ftcs.clone(); //  re-substitution of variable
      t = t + dt;
   }

   // ===  Install method of boundary value ===
   public void setBoundary( double h[] ){
      h[0] = 0.0;    h[Nx-1] = h[Nx-2];
   }

   // ==========================-
   // ===  Calculating method by the second predictor-corrector scheme ===
   public void P_C(){
      for(int i=1; i<Nx-1; i++){
         fp[i] = fa[i]-0.5*Cour*(fa[i+1]-fa[i-1])
                 +diff*(fa[i+1]-2.0*fa[i]+fa[i-1]);
      }
      setBoundary( fp );
      for(int i=1; i<Nx-1; i++){
         fpc[i] = 0.5*(fa[i]+fp[i]-0.5*Cour*(fp[i+1]-fp[i-1])
                  +diff*(fp[i+1]-2.0*fp[i]+fp[i-1]));
      }
      setBoundary( fpc );
      fa=(double[])fpc.clone();
   }

   // === Calculating method using Quick scheme ===
   public void QUICK(){
      for(int i=2; i<Nx-1; i++){
         fquick[i] = fa[i]-Cour*(3.0*fa[i+1]+3.0*fa[i]
                    -7.0*fa[i-1]+fa[i-2])/8.0
                    +diff*(fa[i+1]-2.0*fa[i]+fa[i-1]);
      }
      fquick[1] = 0.0;  setBoundary( fquick );
      fa = (double []) fquick.clone();
   }

   // ===  Calculating method using Kawamara-Kuwahara scheme ===
   public void K_K(){
      for(int i=2; i<Nx-2; i++){
         fkk[i] = fa[i]-Cour*(fa[i+2]-2.0*fa[i+1]+9.0*fa[i]
                 -10.0*fa[i-1]+2.0*fa[i-2])/6.0
                 +diff*(fa[i+1]-2.0*fa[i]+fa[i-1]);
      }
      fkk[1] = 0.0;  fkk[Nx-2] = fkk[Nx-3];  setBoundary( fkk );
      fa = (double[])fkk.clone();
   }

   // ===  Calculating method using CIP scheme ===
   public void solvCip(){
      double xx,fdif,xam1,xbm1;
      // advection phase of Fractional Step
      for(int i=1; i<Nx-1; i++){
         xx = -u0*dt;
         fdif = (fc0[i]-fc0[i-1])/dx;
         xam1 = (g0[i]+g0[i-1]-2.0*fdif)/(dx*dx);
         xbm1 = (-3.0*fdif+2.0*g0[i]+g0[i-1])/dx;
         fc[i] = ((xam1*xx+xbm1)*xx+g0[i])*xx+fc0[i];
         g[i] = (3.0*xam1*xx+2.0*xbm1)*xx+g0[i];
      }
      setBoundary( fc );  setBoundary( g );
      fc0 = (double[])fc.clone(); 	  g0 = (double[])g.clone();
      //  non-advective phase of Fractional Step
      for(int i=1; i<Nx-1; i++){
         fc[i] = fc0[i]+diff*(fc0[i+1]-2.0*fc0[i]+fc0[i-1]);
         g[i]  = g0[i] +diff*(g0[i+1] -2.0*g0[i] +g0[i-1]);
      }
      setBoundary( fc ); setBoundary( g );
      fc0 = (double[])fc.clone();     g0 = (double[])g.clone();
   }

   //  ===  Access methods  ===
   public double[] getRange(){
      double range[]={xmin, xmax, ymin, ymax};
      return range;
   }
   // public double[] getArray( double a[]){ return a; }
   public double[] getXarray(){  return x;  }
   public double[] getFtcsArray(){ return ftcs;  }
   public double[] getCipArray(){ return fc; }
   public double[] getPCarray(){ return fpc; }
   public double[] getQuickArray(){ return fquick; }
   public double[] getKKarray(){ return fkk; }
   public double[] getYsArray(){ return fs; }
   public void setDt( double tt ){ dt = tt; }
   public double getDiff(){ diff = alpha*dt/(dx*dx); return diff; }
   public double getCourant(){ Cour=u0*dt/dx; return Cour; }
}
