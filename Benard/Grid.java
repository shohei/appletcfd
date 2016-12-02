/**====================================================================
   Input data for simulating Benard natural heat convection              
        Grid.java :  Grid generation                
                  All rights reserved, Copyright (C) 2001-2003,  
        Ver.2.0     Last update; December 25,2002, Kiyoshi Minemura  
======================================================================*/

public class Grid {  
   private int    mx, my;       // number of mesh points
   private double x[][], y[][]; // nodal coordinates
   private float xp[], yp[];    // one dimensional nodal coordinates
   private double xmin,xmax,ymin,ymax;
   private String flow;
   
   //    Method for selecting flow object
   public void selectGrid( String flow ){
      this.flow = flow;
      if(flow == "duct")      grid_duct();
      oneDim(); 
   } 
   
   // When flow = "duct", generate the grid for Benard heat convection.
   public void grid_duct(){
      mx=70;  my=16;    // node numbers specified
      x = new double [mx][my]; y = new double [mx][my];
      double a, co, b, sj, s0, s1, s2; 
      int i, j;
	  
      // perpendicular direction ( duct width = 1)
      s0 = 2.0/(double)(my-1);
      a = 1.15;   b = (a+1.0)/(a-1.0);  co = 1.0/(b-1.0);
      for(j=0; j<my; j++){
         sj = s0*(double)j;               // equal ratio
         s1 = Math.pow(b, sj);  s2 = Math.pow(b, sj-1.0);
         y[0][j] = co*(s1-1.0)/(s2+1.0); x[0][j]=0.0;
      }
      // axial direction  ( duct length = 10) 
      s0 = 2.0/(double)(mx-1);
      a = 1.3;  b = (a+1.0)/(a-1.0);   co = 10.0/(b-1.0);
      for(i=0; i<mx; i++){
         sj = s0*(double)i;   // equal ratio
         s1 = Math.pow(b, sj);  s2 = Math.pow(b, sj-1.0);
         x[i][0] = co*(s1-1.0)/(s2+1.0);
      }
      //  
      for(i=1 ; i<mx; i++){
         for(j=0; j<my; j++){ x[i][j]=x[i][0]; y[i][j]=y[0][j];}
      }
      //  graphic output range
      xmin=-0.25; xmax=10.25; ymin=-0.5; ymax=1.5;
   }
   
   // access methods
   private void oneDim(){
      int k, kmax = mx*my;
      xp = new float [kmax];  yp = new float [kmax];
  
      for(int i=0; i<mx; i++){
         for(int j=0; j<my; j++){
            k=my*i+j; 
            xp[k] = (float)x[i][j];  yp[k] = (float)y[i][j];
         }
      }
   }
   public double [] getRange(){
      double ra[] = { xmin, xmax, ymin, ymax };
      return ra;
   }
   public double[][] getXarray(){ return x; }
   public double[][] getYarray(){ return y; }
   public float [] getXpArray(){ return xp; }
   public float [] getYpArray(){ return yp; }
   public int [] getGrid(){ 
      int co[]={ mx, my };  return co;
   }
}      
