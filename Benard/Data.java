/**=====================================================================
    Input data for simulating Benard natural heat convection             
          　Data.java :  boudary and initial condition　　         　
                     All rights reserved, Copyright (C) 2001-2003,  
         Ver.2.0     Last update: December 25, 2002, Kiyoshi Minemura
=======================================================================*/

public class Data {
   public static double Re=75.0;    // Reynolds number(Re)
   public static double Ra=10000.0; // Rayleigh number 
   public static double Pr=1.0;     // Prandtle number

   private double u[][], v[][];    // u,v=velocity in x, y direction
   private double dt=0.01;         // time step
   private double Tsouth=1.0;      // temparature at bottom
   private double repL=1.0;        // representative length
   // public double rho=998.2;     // density [Kg/m^3]
   // public double mu=1.0e-3;     // viscosity
   // public double gravity=9.8;   // gravity
   // public double alp, sph=4.1816e+3; // specific heat
   // public double tcn=5.94;      // thermal coefficient of nutural convection
   // public double Ymax=0.1, Xmax=1.0, deltaT=40.0;
   // public double beta=0.207e-3; // expansion coefficient [1/K]
   private String BClef[], BCrig[], BClow[], BCupp[]; // boundary condition
   private int  NSmax;            // maximum number of calculation step
   private int  mx, my;           // number of nodes
   private float up[], vp[];      // velocity
   private String flow;
	
   public void selectData( String flow ){ 
      this.flow = flow;
   } 
   public void setGrid( int co[] ){
      mx = co[0];  my = co[1];
      u = new double[mx][my];   v = new double[mx][my]; 
      BClef = new String[my];   BCrig = new String[my];
      BClow = new String[mx];   BCupp = new String[mx];
	   
      if(flow == "duct")        BC_duct(  );
      oneDim();
   }
	
   //  When flow = "duct", input the boundary condition.
   private void BC_duct(){
      int i, j;		
      NSmax = (int)(5000/dt);  repL = 0.1;
		
      // boundary conditions of lower side of computational domain
      for(i = 0; i<mx; i++)   BClow[i] = "wall";  
      // boundary conditions of upper side 
      for(i = 0; i<mx; i++)   BCupp[i] = "wall";
      // boundary conditions of left side
      for(j = 0; j<my; j++)   BClef[j] = "wall";
      // boundary conditions of right side
      for(j = 0; j<my; j++)   BCrig[j] = "wall";
   }

   // access methods
   private void oneDim(){
      int k, kmax = mx*my;
      up = new float [kmax];  vp = new float [kmax];
	  
      for(int i=0; i<mx; i++){
         for(int j=0; j<my; j++){
            k=my*i+j; 
            up[k] = (float)u[i][j];  vp[k] = (float)v[i][j];
         }
      }
   }
   public float [] getUParray(){ return up; }
   public float [] getVParray(){ return vp; }
   public int  getNSmax(){  return NSmax;  }
   public String [] getBClow(){ return BClow; }
   public String [] getBCupp(){ return BCupp; }
   public String [] getBClef(){ return BClef; }
   public String [] getBCrig(){ return BCrig; } 
   public double [] getParameters(){
      double Co[] = { dt, repL, Tsouth };
      return Co; 
   }
}
