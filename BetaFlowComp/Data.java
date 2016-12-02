/************************************************************************/
/*    Input data for simulating incompressible fluid flow               */
/*          　DATA class: boudary and initial condition　　         　　*/
/*                           Last update: Novembe 9, 2001.              */
/*          All rights reserved, Copyright (C) 2001, Kiyoshi Minemura   */
/************************************************************************/
import java.awt.*;

//  ===== Set boundary and initial conditions ==========================　
public class Data {
   public static double Re=300.0; // nu=inverse of Reynolds number(Re)
   public double u[][], v[][];    // u,v=velocity in x, y direction
   private double dt=0.005;       // time step
   private double Umean=1.0;       // dimensionless inlet flow velocity
   private double p0=1.0;          // inlet pressure
   private double nu=1.04e-6;      // 
   private double repL=1.0;        //
   public String BClef[], BCrig[], BClow[], BCupp[]; // boundary condition
   private int  NSmax;            // maximum number of calculation step
   private int    mx, my, mx1, my1, mx2, my2, mxm, mym; // number of nodes
   private float uz[], vz[];           // velocity
   private float xrange [] = new float [2];
   private float yrange [] = new float [2]; // range of drawing
	
   public void Data_select(Grid gr, String flow){  
      mx = gr.getMx();  my = gr.getMy();
      mxm = mx-1;  mym = my-1;  mx1 = mx+1;  my1 = my+1; 
      mx2 = mx+2;  my2 = my+2;
      u = new double[mx][my];   v = new double[mx][my]; 
      BClef = new String[my];   BCrig = new String[my];
      BClow = new String[mx];   BCupp = new String[mx];
	   
      if(flow == "duct")          BC_duct( gr );
      else if(flow == "bend")     BC_bend( gr );
      else if(flow == "cylinder") BC_cylinder( gr );
      oneDim();
   } 
	
   //   Example No. 1 (flow="pipe") 
   private void BC_duct(Grid gr){
      int i, j;		
      NSmax=(int)(20/dt);  repL=0.1;
		
      // boundary conditions of lower side of computational domain
      for(i = 0; i<mx; i++)   BClow[i]="wall";  
      // boundary conditions of upper side 
      for(i = 0; i<mx; i++)   BCupp[i]="wall";
      // boundary conditions of left side
      for(j = 0; j<my; j++)   BClef[j]="in";
      // boundary conditions of right side
      for(j = 0; j<my; j++)   BCrig[j]="out";
	
      // initial conditions
      //for(j = 1; j <mym; j++){  
      //   u[0][j]=Umean; //2.0*Umean*(1.0-4.0*(gr.y[0][j]-0.5)*(gr.y[0][j]-0.5));
      //}
      for(i = 0; i < mx; i++){   
         for(j = 1; j<mym; j++)  u[i][j]=Umean;
      }
      // xmin=1.6; ymin=-3; xmax=8.6; ymax=3;
   } // end of BC_pipe
	
   //   Example No. 2 (flow="bend") 
   private void BC_bend(Grid gr){
      int i, j;
      double theta, rate, um, vm, pi=Math.PI;
      NSmax=(int)(30.0/dt);  repL=0.1;
	
      // boundary conditions of lower side of computational domain
      for(i = 0; i<mx; i++)   BClow[i]="wall";  
      // boundary conditions on upper side 
      for(i = 0; i<mx; i++)   BCupp[i]="wall";
      // boundary conditions on left side
      for(j = 0; j<my; j++)   BClef[j]="in";
      // boundary conditions of right side
      for(j = 0; j<my; j++)   BCrig[j]="out";
		
      // initial conditions; all those on other sides are zero.      
      rate=0.5*pi/(double)(gr.mxRout-gr.mxRin);  
      for(i = 0; i <= gr.mxRout; i++){      
         if(i <= gr.mxRin){   
            for(j = 1; j < mym; j++) 	  u[i][j]=Umean;
            // u[i][j]=2.0*Umean*(1.0-4.0*(gr.y[1][j]+1.0)*(gr.y[1][j]+1.0)); 
         }else{
            theta=rate*(double)(i-gr.mxRin);
            for(j = 1; j< mym; j++)  u[i][j]=u[0][j]*Math.cos(theta);
         }
      }
		
      for(i = gr.mxRin+1; i < mx; i++){  
         if(i <= gr.mxRout){
            theta=rate*(double)(i-gr.mxRin);
            for(j = 1; j < mym; j++)  v[i][j]=u[0][j]*Math.sin(theta);
         }else{
            for(j = 1; j < mym; j++)  v[i][j]=u[0][j];
         }
      }
   } // end of BC_bend 

   //   Example No. 3 (flow="cylinder") 
   private void BC_cylinder(Grid gr){
      int i, j, j1, j2;
      NSmax = (int)(500.0/dt); repL=0.1;
      // Boundary conditions of lower side of computational domain
      for(i = 0; i<gr.mxL1; i++)        BClow[i]="periodic"; 
      for(i =gr.mxL1; i<gr.mxL2; i++)   BClow[i]="wall";
      for(i =gr.mxL2; i<mx; i++)        BClow[i]="periodic";
      // Boundary conditions on upper side 
      for(i = 0; i<gr.mxU1; i++)        BCupp[i]="slip";
      for(i =gr.mxU1; i<gr.mxU2; i++)   BCupp[i]="in";
      for(i =gr.mxU2; i<mx; i++)        BCupp[i]="slip";
      // Boundary conditions on left side
      for(j = 0; j<my; j++)        BClef[j]="out";
      // Boundary conditions of right side
      for(j = 0; j<my; j++)        BCrig[j]="out";		

      // Initial conditions;       
      for(i=0; i<mx; i++){         
         if(BClow[i]=="wall") j1=1;  else j1=0;
         for(j=j1; j<my; j++) u[i][j]=Umean; 
      }
      u[gr.mxL2][0]=0.0;
      xrange[0]=1.0f; yrange[0]=-4.5f; xrange[1]=10.0f; yrange[1]=4.5f;     

   } // end of BC_cylinder  	

   private void oneDim(){
      int k, kmax=mx*my;
      uz=new float [kmax]; vz=new float [kmax];
	  
      for(int i=0; i<mx; i++){
         for(int j=0; j<my; j++){
            k=my*i+j; uz[k]=(float)u[i][j]; vz[k]=(float)v[i][j];
         }
      }
   }
   public float [] getXrange(){ return xrange; }
   public float [] getYrange(){ return yrange; }
   public float [] getUarray(){ return uz; }
   public float [] getVarray(){ return vz; }
   public int  getNSmax(){  return NSmax;  }
   public String [] getBClow(){ return BClow; }
   public String [] getBCupp(){ return BCupp; }
   public String [] getBClef(){ return BClef; }
   public String [] getBCrig(){ return BCrig; } 
   public double [] getConsts(){
      double Co[]=new double [3];
      Co[0]=nu; Co[1]=dt; Co[2]=repL; 
      return Co; 
   }
}