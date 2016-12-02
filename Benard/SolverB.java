/**=====================================================================
     Solver for simulating Benard natural heat convection                                 　SolverB.java : fractional step method     　　         　　
     All rights reserved, Copyright (C) 2001-2003, Kiyoshi Minemura 
     Ver.2.0     Last update: December 25, 2002. 
=======================================================================*/
import java.awt.*;

public class SolverB implements Runnable {
   private Metric me;     private DrawCanvasB drw;  
   private Graphics g,bg;  
   private String flow, graph="temper", scheme="QUICK";
   private float up[], vp[];
   private String BClow[], BCupp[], BClef[], BCrig[];
   private Thread solver_thread;   

   private double U[][], V[][], Ua[][], Va[][], ua[][], va[][], u[][], v[][];
   private double P[][], D[][], Gx[][], Gy[][], nu, dt, Uout, repL, Lap;
   private double T[][], Gt[][], Tsouth, invRe, alph;
   private double x[][], y[][];
   private int    mx, my, mx1, my1, mx2, my2, mxm, mym, kmax;
   private int    NSmax, NStemp, IterNum;   
   
   // Constructor
   public SolverB( DrawCanvasB drw, Metric me ){
      this.drw = drw;    this.me = me; 
   }

   //  Access methods for initialization
   public void setGrid( int co[] ){
      mx = co[0]; my = co[1];
      mx1=mx+1;  my1=my+1;  mx2=mx+2;  my2=my+2; 
      mxm=mx-1;  mym=my-1;  kmax=mx*my;  
      U =new double[mx][mym]; V =new double[mxm][my]; // contravariant velocity
      Ua=new double[mx][mym]; Va=new double[mxm][my]; // auxiliary contravariant velocity
      ua=new double[mx2][my1];va=new double[mx1][my2]; // auxiliary velocity	  
      u =new double[mx2][my1];v =new double[mx1][my2];// physical velocity（staggared)	  
      D =new double[mxm][mym];P =new double[mx1][my1];// residual of continuity 
      up = new float[kmax];   vp= new float[kmax];   // physical velocity at node
      Gx = new double[mx][my];Gy = new double [mx][my];
      Gt = new double[mx][my];T  = new double [mx1][my1]; 
   }
   public void setNSmax(int N ){  NSmax = N; }
   public void setBC( String up[], String low[],
                      String rig[], String lef[]){
      BCupp = up;  BClow = low;  BCrig = rig;  BClef = lef;
   }
   public void setParameters( double Co[] ){
      dt=Co[0]; repL=Co[1]; Tsouth=Co[2];
   }
   
   //  Initialization	
   public void initSolver( String flow){
      this.g = drw.getGraphics();  this.bg = drw.bg;
      NStemp = 0;  IterNum = 0;  

      // initilaization of flow field  
      for(int i=0; i<mx1; i++) T[i][0] = Tsouth;

      drw.deleteBuffer();
      solver_thread = new Thread(this);
   }

   public void startSolver(){
      IterNum = NStemp;
      solver_thread = new Thread(this);
      solver_thread.start();  // compute method "run"
   }

   public void stopSolver(){   // Interruption of computation
      if( solver_thread != null ){
         NStemp = IterNum;
         solver_thread = null;
         if(IterNum < NSmax ) g.drawImage(drw.getBuffer(), 0, 0, drw);
      }
   }
   public void stop(){
      solver_thread = null; 
   }

   //  ==========================================================
   //  Iteration of fractional step method
   //  ==========================================================
   public void run(){
      int    NoP=0;  
      while( solver_thread != null ){  
         // ===== First step of FSmethod =========
         invRe=1.0/Data.Re;
		
         calc_contraV(); // calculate contravariant velocity
		 
         calc_auxU();    // calculate auxiliary velocity ua
			
         calc_auxV();    // calculate auxiliary velocity va in direction of ETA

    	 // ===== Second step of FS method  =========
         calc_D();
		 
         //  Solve Poisson's equation of pressure (SOR method)
         NoP=S_Poisson(me.D1P, me.D2P, me.D3P, me.D4P, me.D5P, me.D6P, me.D7P, 
                       me.D8P, me.D0P);

         update_Vel();
		 
         calc_energy();

         //  ===== Output graphics  (once ten iteretions) ==========
         if( (IterNum % 10) == 0)   output_graph(NoP);   
		 
         if(IterNum >= NSmax ){
            solver_thread = null;   
            g.setColor(Color.red);   g.drawString("END", 20,27);
         }
         IterNum++;
      }
   }
   
   //========================================================================
   //   Calculate contravariant velocity on all staggared mesh points
   private void calc_contraV(){
      int i,j,ip,jp,im,ip1,jp1;
      double u0, u1, u2, u3, u4, v0;
      for(i=0; i<mx; i++){      ip=i+1;  
         for(j=0; j<mym; j++){  jp=j+1; jp1=jp+1;
            v0=0.25*(v[i][jp]+v[i][jp1]+v[ip][jp]+v[ip][jp1]); 				   
            U[i][j]=me.JW[i][j]*(me.XIxW[i][j]*u[ip][jp]+me.XIyW[i][j]*v0);
         }
      }

      for(i=0; i<mxm; i++){     ip=i+1; ip1=ip+1; im=i-1;
         for(j=0; j<my; j++){   jp=j+1; 
            //  u0=0.25*(u[ip][j]+u[ip][jp]+u[ip1][j]+u[ip1][jp]);  
            u1=u[ip][j]; u2=u[ip][jp]; u3=u[ip1][j]; u4=u[ip1][jp];
            if(j == 0){
               if(BClow[i]=="periodic" && BClow[ip]=="wall") u3=u[mx1-ip1][1];
               if(BClow[i]=="wall" && BClow[ip]=="periodic") u3=-u4;
            }else if(j == mym){
               if(BCupp[i]=="slip" && BCupp[ip]=="in")  u4=u3;
               if(BCupp[i]=="in" && BCupp[ip]=="slip")  u4=u[i][my];
            }
            u0=0.25*(u1+u2+u3+u4);
            V[i][j]=me.JS[i][j]*(me.ETAxS[i][j]*u0+me.ETAyS[i][j]*v[ip][jp]);
         }
      }  
   }
   
   //=========================================================================
   //  Calculate auxiliary velocity (ua) in direction of XI 
   private void calc_auxU(){
      int ip, jp, ip1, jp1, im, jm, ip2, jp2;
      double Uw, Ue, Vn, Vs, u0, u1, u2, u3, u4, u5, u6, u7, u8;
      double ue, uw, un, us, aa=0.0, a1, a2, bb;
      for(int i = 1; i<mxm; i++){     ip=i+1; im=i-1; ip1=ip+1; ip2=ip1+1;
         for(int j = 0; j<mym; j++){  jp=j+1; jm=j-1; jp1=jp+1; jp2=jp1+1;
            Uw=0.5*(U[im][j]+U[i][j]);  Ue=0.5*(U[i][j]+U[ip][j]); 
            Vs=0.5*(V[im][j]+V[i][j]);  Vn=0.5*(V[im][jp]+V[i][jp]);
            u1=u[i][j];    u2=u[i][jp];    u3=u[i][jp1]; 
            u4=u[ip][j];   u0=u[ip][jp];   u5=u[ip][jp1]; 
            u6=u[ip1][j];  u7=u[ip1][jp];  u8=u[ip1][jp1];
            if(j == 0){
               if(BClow[i]=="periodic" && BClow[ip]=="wall") u6=u[mx1-ip1][1];
               if(BClow[i]=="wall" && BClow[ip]=="periodic") u6=-u7;
               if(BClow[im]=="wall" && BClow[i]=="periodic") u4=-u0;
            }else if(j == mym-1){
               if(BCupp[i]=="slip" && BCupp[ip]=="in")     u8=u[ip1][my];
               if(BCupp[i]=="in" && BCupp[ip]=="slip")     u8=u[ip][my];
               if(BCupp[i]=="slip" && BCupp[im]=="in")     u5=u[i][my];
            }
            if(scheme == "donor-cell"){
               //  ** donor cell sheme (upstream scheme) **
               aa=0.5*(Ue*(u0+u7)+Math.abs(Ue)*(u0-u7)-Uw*(u2+u0)-Math.abs(Uw)*(u2-u0)
                  +Vn*(u0+u5)+Math.abs(Vn)*(u0-u5)-Vs*(u4+u0)-Math.abs(Vs)*(u4-u0));
            }else if(scheme == "first-order"){ // first order upstream scheme
               Uw=U[i][j]; Ue=Math.abs(Uw);
               Vs=0.25*(V[im][j]+V[im][jp]+V[i][j]+V[i][jp]); 
               Vn=Math.abs(Vs);
               aa=0.5*(Uw*(u7-u2)-Ue*(u2-2.0*u0+u7)+Vs*(u5-u4)-Vn*(u4-2.0*u0+u5));
            }else if(scheme == "QUICK"){    // QUICK scheme
               uw=u[im][jp]; ue=u[ip2][jp]; a1=0.0; a2=0.0;			      
               aa=Ue*(-u2+9.0*(u0+u7)-ue)+Math.abs(Ue)*(-u2+3.0*(u0-u7)+ue)
                  -Uw*(-uw+9.0*(u2+u0)-u7)-Math.abs(Uw)*(-uw+3.0*(u2-u0)+u7);
               if(j==0){
                  un=u[ip][jp2]; 
                  if(BClow[i]=="wall") a2=0.0; 
                  else if(BClow[i]=="periodic"){
                     if(BClow[im]=="wall") a2=0.0;
                     else{ 
                        us=u[mx1-ip][2];
                        a2=Vs*(-us+9.0*(u4+u0)-u5)+Math.abs(Vs)*(-us+3.0*(u4-u0)+u5);
                     } 
                  }
                  a1=Vn*(-u4+9.0*(u0+u5)-un)+Math.abs(Vn)*(-u4+3.0*(u0-u5)+un);
               }else if(j==mym-1){
	          us=u[ip][jm];
                  if(BCupp[i]=="wall"){  a1=0.0; }
                  else{
                     un=u5; 
                     a1=Vn*(-u4+9.0*(u0+u5)-un)+Math.abs(Vn)*(-u4+3.0*(u0-u5)+un);
                  }
                  a2=Vs*(-us+9.0*(u4+u0)-u5)+Math.abs(Vs)*(-us+3.0*(u4-u0)+u5);
               }else{ 
                  us=u[ip][jm];  un=u[ip][jp2];
                  a1=Vn*(-u4+9.0*(u0+u5)-un)+Math.abs(Vn)*(-u4+3.0*(u0-u5)+un);
                  a2=Vs*(-us+9.0*(u4+u0)-u5)+Math.abs(Vs)*(-us+3.0*(u4-u0)+u5);
               }
               aa=(aa+a1-a2)/16.0;
            }
            // diffusion term
            Lap=me.D1u[i][j]*u1+me.D2u[i][j]*u2+me.D3u[i][j]*u3
                +me.D4u[i][j]*u4+me.D0u[i][j]*u0+me.D5u[i][j]*u5
                +me.D6u[i][j]*u6+me.D7u[i][j]*u7+me.D8u[i][j]*u8;
            // auxiliary velocity ua
            bb = -(aa-invRe*Lap)/me.JW[i][j];
            ua[ip][jp] = u0+dt*(3.0*bb-Gx[i][j])*0.5;
            Gx[i][j] = bb;
         }
      }
      S_BCfor_ua();
   }

   private void S_BCfor_ua(){
      int    i, j, jp;
      double u0;
      //   right-hand boundary  (outflow) 
      for(j = 0; j<mym; j++){ jp=j+1;
         //u0=u[mxm][jp];  
         ua[mx][jp]=0.0; //u[mx][jp]-dt*Uout*(u[mx][jp]-u0)/me.JW[mxm][j];
         ua[mx+1][jp]=ua[mx-1][jp];
      }
      for(i=0; i<mx; i++){
         ua[i+1][0]=-ua[i+1][1]; ua[i+1][my]=-ua[i][my-1];  
      }
      //   left-hand boundary of computational domain (inflow or outflow ) 
      for(j = 0; j<mym; j++){      jp=j+1;    
      //if(BClef[j] == "in")      ua[1][jp]=u[1][jp]; 
      //else  ua[1][jp]=u[1][jp]-dt*Uout*(u[1][jp]-u[2][jp])/me.JW[0][j];
	 ua[1][jp]=0.0;  ua[0][jp]=ua[2][jp];
      }  	   
   }
   
   //=====================================================================
   //            auxiliary velocity va in direction of ETA
   private void calc_auxV(){
      int im, jm, ip, jp, ip1, jp1, j0, jp2, ip2;
      double Uw, Ue, Vn, Vs, v0, v1, v2, v3, v4, v5, v6, v7, v8;
      double ve, vw, vn, vs, aa=0.0, a1, a2, bb;
      double 	 Gr3=Data.Ra/Data.Pr*invRe*invRe;
      for(int i = 0; i<mxm; i++){      im=i-1; ip=i+1; ip1=ip+1; ip2=ip1+1;
         if(BClow[i]=="periodic" && i<mx/2) j0=0;  else j0=1;
         for(int j = j0; j<mym; j++){  jm=j-1; jp=j+1; jp1=jp+1; jp2=jp1+1;
            if(j == 0){ 
               Uw=0.5*(U[i][0]-U[mxm-i][0]); Ue=0.5*(U[ip][0]-U[mxm-ip][0]); 
               Vs=0.5*(V[i][0]-V[mxm-ip][1]); Vn=0.5*(V[i][0]+V[i][1]); 
               if(BClow[ip]=="wall") Ue=0.0; //0.5*(U[ip][0]-U[mxm-ip][0]);
            }else{  
                  Uw=0.5*(U[i][jm]+U[i][j]);   Ue=0.5*(U[ip][jm]+U[ip][j]); 
                  Vs=0.5*(V[i][jm]+V[i][j]);   Vn=0.5*(V[i][j]+V[i][jp]);
            }
            v1=v[i][j];     v2=v[i][jp];    v3=v[i][jp1]; 
            v4=v[ip][j];    v0=v[ip][jp];   v5=v[ip][jp1]; 
            v6=v[ip1][j];   v7=v[ip1][jp];  v8=v[ip1][jp1];
            if(scheme == "donor-cell"){      // donor cell scheme
               aa=0.5*(Ue*(v0+v7)+Math.abs(Ue)*(v0-v7)-Uw*(v2+v0)-Math.abs(Uw)*(v2-v0)
                  +Vn*(v0+v5)+Math.abs(Vn)*(v0-v5)-Vs*(v4+v0)-Math.abs(Vs)*(v4-v0));
            }else if(scheme == "first-order"){ // first order upstream scheme
               if(j == 0)  Uw=0.25*(U[i][0]+U[ip][0]-U[mxm-i][0]-U[mxm-ip][0]);
               else        Uw=0.25*(U[i][jm]+U[i][j]+U[ip][jm]+U[ip][j]);
               Ue=Math.abs(Uw);
               Vs=V[i][j]; Vn=Math.abs(Vs);
               aa=0.5*(Uw*(v7-v2)-Ue*(v2-2.0*v0+v7)+Vs*(v5-v4)-Vn*(v4-2.0*v0+v5));
            }else if(scheme == "QUICK"){   // QUICK scheme 
               if(i == 0)       { vw=v2;         ve=v[3][jp]; }
               else if(i==mxm-1){ vw=v[im][jp];  ve=v7; }
               else{              vw=v[im][jp];  ve=v[ip2][jp];}
               aa=Ue*(-v2+9.0*(v0+v7)-ve)+Math.abs(Ue)*(-v2+3.0*(v0-v7)+ve)
                  -Uw*(-vw+9.0*(v2+v0)-v7)-Math.abs(Uw)*(-vw+3.0*(v2-v0)+v7);
               vn=v[ip][jp2];
               if(j==0)  vs=v[mx-ip][3]; else  vs=v[ip][jm];  
               a1=Vn*(-v4+9.0*(v0+v5)-vn)+Math.abs(Vn)*(-v4+3.0*(v0-v5)+vn);
               a2=Vs*(-vs+9.0*(v4+v0)-v5)+Math.abs(Vs)*(-vs+3.0*(v4-v0)+v5);
               aa=(aa+a1-a2)/16.0; 
            }
            // diffusion term
            Lap=me.D1v[i][j]*v1+me.D2v[i][j]*v2+me.D3v[i][j]*v3
                +me.D4v[i][j]*v4+me.D0v[i][j]*v0+me.D5v[i][j]*v5
                +me.D6v[i][j]*v6+me.D7v[i][j]*v7+me.D8v[i][j]*v8;
            // auxiliary velocity va
            bb=-(aa-invRe*Lap)/me.JS[i][j]+Gr3*T[ip][jp];
            va[ip][jp]=v0+dt*(3.0*bb-Gy[i][j])*0.5;
            Gy[i][j]=bb;
         }
      } 
      S_BCfor_va();
   }   
   private void S_BCfor_va(){
      int    i, j, ip, jp;
      double v0;
      //   lower boundary (wall or periodic )
      for(i = 0; i<mxm; i++){     ip=i+1;
         if(BClow[i] == "wall"){   va[ip][1]=0.0; va[ip][0]=va[ip][2];}
         // else if(i>mx/2 && BClow[i]=="periodic") va[ip][1]=va[mx-ip][1]; 
      }
      //   upper boundary (wall, in or slip )
      for(i = 0; i<mxm; i++){ ip=i+1;
         // if(BCupp[i]=="in")   va[ip][my]=v[ip][my]; 
         // else                 va[ip][my]=0.0;
        va[ip][my]=0.0;  va[ip][my+1]=va[ip][my-1];
      }
      //   left-hand boundary of computational domain (inflow or outflow ) 
      for(j = 0; j<my; j++){      jp=j+1;    
         //if(BClef[j] == "in")     va[0][jp]=v[0][jp]; 
         //else  va[0][jp]=v[0][jp]-dt*Uout*(v[0][jp]-v[1][jp])/me.JS[0][j];
         va[0][jp]=-v[1][jp];
      }
      //   right-hand boundary  (outflow) 
      for(j = 0; j<my; j++){ jp=j+1;
         // v0=v[mxm][jp];  
         // va[mx][jp]=v[mx][jp]-dt*Uout*(v[mx][jp]-v0)/me.JS[mxm-1][j];
         va[mx][jp]=-va[mx-1][jp];
      }
   }
   
   //========================================================================
   // calculate contravariant auxiliary velocity (Ua, Va)
   private void calc_D(){
      int i, j, ip, jp, ip1, jp1, j0;
      double u0, v0;
      for(i=0; i<mx; i++){      ip=i+1; ip1=ip+1;
         for(j=0; j<mym; j++){  jp=j+1; jp1=jp+1; 
            v0=0.25*(va[i][jp]+va[i][jp1]+va[ip][jp]+va[ip][jp1]);
            Ua[i][j]=me.JW[i][j]*(me.XIxW[i][j]*ua[ip][jp]+me.XIyW[i][j]*v0); 
         }
      }
      for(i=0; i<mxm; i++){      ip=i+1; ip1=ip+1;
         if(BClow[i]=="periodic" && i<mx/2) j0=0;  else j0=1;
         for(j=j0; j<mym; j++){  jp=j+1;
            if(j==0) u0=0.25*(ua[ip][1]+ua[ip1][1]+ua[mx1-ip][1]+ua[mx1-ip1][1]);
            else     u0=0.25*(ua[ip][j]+ua[ip][jp]+ua[ip1][j]+ua[ip1][jp]);
            Va[i][j]=me.JS[i][j]*(me.ETAxS[i][j]*u0+me.ETAyS[i][j]*va[ip][jp]);          }
      }
      S_BCfor_Ua();
      //   Calculate residual of continuity condition D
      for(i = 0; i<mxm; i++){     
         for(j = 0; j<mym; j++){  
            D[i][j]=-Ua[i][j]+Ua[i+1][j]-Va[i][j]+Va[i][j+1];
         }
      }
   }

   //　== Reinstall boundary values for contravariant velocity ==========
   private void S_BCfor_Ua(){
      int    i, j;	   
      //   lower boundary
      for(i = 0; i<mxm; i++){
         if(BClow[i] == "wall")  Va[i][0]=0.0; 
         // else if(BClow[i]=="periodic" && i>mx/2) Va[i][0]=-Va[mxm-1-i][0];  
      }
      //    upper boundary
      for(i = 0; i<mxm; i++){
         if(BCupp[i]=="wall" || BCupp[i]=="slip") Va[i][mym]=0.0;
         // else if(BCupp[i]=="in")                  Va[i][mym]=V[i][mym]; 
      }
      //   left-hand boundary of computational domain 
      for(j = 0; j<mym; j++){          
         //if(BClef[j] == "in")  Ua[0][j]=U[0][j];
         Ua[0][j]=0.0; Ua[mx-1][j]=0.0;
      } 
      //   right-hand boundary 
   }   
   
   //========================================================================
   // ====== SOR method for Poisson equation of pressure =====================
   private int S_Poisson( double D1[][], double D2[][], double D3[][], 
                   double D4[][], double D5[][], double D6[][], double D7[][], 
                   double D8[][], double D0[][]){
      double diff, ERR, ERRmax, ERRabs, error, Pij, aa, dti=1.0/dt;
      int    i, j, ip, jp, im, jm, iERR, jERR, NoP = 0;
      int    imax ;      // Maximum iteration number for Poisson equation
      double eps=1.0e-7, eps0=1.0e-5;   // Maximum allowance for error  
      double omega = 1.6;     // Relaxation coefficient for SOR method
      boolean flag = true; 
      if(IterNum < 19) imax = 200; else imax = 100;
      if(IterNum == 0){
         for(i=0; i<mx1; i++){ for(j=0; j<my1; j++) P[i][j]=10.0; }
      }

      ERRmax = 0.0; 
      while( flag ){  
         ERR = 0.0;  ERRabs = 0.0;  
         for(i=1; i<mx; i++){     im=i-1;  ip=i+1;
            for(j=1; j<my; j++){  jm=j-1;  jp=j+1;
               Pij=P[i][j];
               aa=D1[im][jm]*P[im][jm]+D2[im][jm]*P[im][j]+D3[im][jm]*P[im][jp]
                  +D4[im][jm]*P[i][jm] +D5[im][jm]*P[i][jp]+D6[im][jm]*P[ip][jm]
                  +D7[im][jm]*P[ip][j] +D8[im][jm]*P[ip][jp];
               diff = (dti*D[im][jm]-aa)/D0[im][jm]-Pij;
               ERR = ERR+diff*diff;  
               error = Math.abs(diff); 
               if(error > ERRabs)  ERRabs = error; //{ ERRmax = ERR;  iERR=i;  jERR=j; }
               P[i][j] = Pij+omega*diff;
            }
         }
         BCpres();  // reinstall boundary values of pressure
         if( ERR < eps0 )     flag = false; 
         // if(ERRabs < eps0) flag = false; 
         NoP++;
         if(NoP >= imax)      flag = false;
         // if(Math.abs(ERRmax-ERRabs) < eps)  flag=false;
         ERRmax = ERRabs;
      } 
      return NoP;
   } 
   
   //  Update pressure on boundaries  
   private void BCpres(){
      int    i, j, ip, jp;
      // left sides boundary (inflow or outflow )
      for(j = 0; j<mym; j++){ jp=j+1;  P[0][jp]=P[1][jp];   }
      // right sides boundary ( outflow )
      for(j = 0; j<mym; j++){ jp=j+1;  P[mx][jp]=P[mxm][jp]; }
      // lower side boundary ( wall or periodic ) 
      for(i = 0; i<mxm; i++){   ip=i+1;
         if(BClow[i] == "wall")        P[ip][0]=P[ip][1]; 
         else if(BClow[i]=="periodic") P[ip][0]=P[mx-ip][1];
      }
      // upper side boundary  ( wall, in or slip )
      for(i = 0; i<mxm; i++){    ip=i+1;  P[ip][my]=P[ip][mym];}
      // corner
      if(BClow[0]=="periodic")   P[0][0]=P[mx][1]; else P[0][0]=P[0][1];  
      if(BClow[mxm]=="periodic") P[mx][0]=P[0][1]; else P[mx][0]=P[mx][1];
      P[0][my]=P[0][mym]; 	     P[mx][my]=P[mx][mym];  
   }  // end of BCpres
   
   //========================================================================
   //   Update velocity (U, u) from (Ua, ua)
   private void update_Vel(){
      int i, j, im, jm, ip, jp, ip1, jp1, j0, i0;
      double dpx, dpy;
      if(BClef[0]=="in") i0=1; else i0=0;						   
      for(i=i0; i<mx; i++){      im=i-1; ip=i+1;
         for(j=0; j<mym; j++){   jp=j+1; jp1=jp+1;
            dpx=-P[i][jp]+P[ip][jp];
            dpy=0.25*(-P[i][j]+P[i][jp1]-P[ip][j]+P[ip][jp1]);
            u[ip][jp]=ua[ip][jp]-dt*(me.XIxW[i][j]*dpx+me.ETAxW[i][j]*dpy);
         }
      }
      //                 (V) from (Va)
      for(i=0; i<mxm; i++){      ip=i+1; ip1=ip+1;
         if(BClow[i]=="periodic" && i<mx/2) j0=0;   else j0=1;
         for(j=j0; j<mym; j++){   jp=j+1;  
            dpx=-0.25*(P[i][j]+P[i][jp]-P[ip1][j]-P[ip1][jp]);
            dpy=-P[ip][j]+P[ip][jp];
            if(j==0 && BClow[i]=="periodic" && BClow[ip]=="wall")
               dpx=-0.25*(P[i][j]+P[i][jp]-P[ip1][jp]-P[mx-ip1][1]);
               v[ip][jp]=va[ip][jp]-dt*(me.XIyS[i][j]*dpx+me.ETAyS[i][j]*dpy); 
         }
      }
      S_updateBC();  //   Update boundary velue of (u, v)
   }

   //　 Update velocity on boundaries for next iteration ===========
   private void S_updateBC(){
      int    i, j, ip, jp, ip1;
      //   lower side boundary (wall or periodic ) 
      for(i = 0; i<mxm; i++){  ip=i+1; ip1=ip+1;  
         // if(BClow[i]=="periodic"){ u[ip][0]=u[mx1-ip][1]; v[ip][0]=v[mx-ip][2]; 
         //    if(i>mx/2)         v[ip][1]=v[mx-ip][1];
         // }else if(BClow[i]=="wall"){
         v[ip][0]=v[ip][2]; u[ip][0]=-u[ip][1]; //}
      }
      if(BClow[mxm]=="wall") u[mx][0]=-u[mx][1]; else u[mx][0]=u[1][1];
      //   upper side boundary (wall, inflow or slip )
      for(i = 0; i<mxm; i++){  ip=i+1;  ip1=ip+1;
         if(BCupp[i]=="wall"){
            v[ip][my1]=v[ip][mym]; u[ip][my]=-u[ip][mym]; // v[ip][my]=0.0;
         }// else if(BCupp[i]=="slip"){
          //    v[ip][my]=0.0; v[ip][my1]=-v[ip][mym]; u[ip][my]=u[ip][mym];}
      }
      if(BCupp[mxm]=="wall") u[mx][my]=-u[mx][mym]; else u[mx][my]=u[mx][mym];
      //   left side boundary ( slip )// ( inflow or outflow ) 
      for(j=0; j<mym; j++){   jp=j+1;
         // if(BClef[j]=="out"){ u[0][jp]=u[1][jp]; v[0][jp]=va[0][jp];} 
         u[0][jp]=u[2][jp]; v[0][jp]=-v[1][jp];				              }
      if(BClef[mym]=="out") v[0][my]=v[1][my];
      //   right side boundary( slip ) //( outflow )
      for(j = 0; j<mym; j++){   jp=j+1;  
         //  u[mx1][jp]=u[mx][jp]; v[mx][jp]=v[mxm][jp];
         u[mx1][jp]=u[mx-1][jp]; v[mx][jp]=-v[mxm][jp];
      }
      // cornar
      // v[0][0]=-v[mx][2];  v[mx][my]=v[mxm][my];   
   }
   
   //========================================================================
   // control of off-screen image
   private void output_graph(int NoP){
      // int i, j, im, jm, ip, jp, k, graph_tool=3;
      drw.deleteBuffer();
      oneDimU();  //  pp=oneDimT();  
      drw.drawBenard( up, vp, oneDimT(), vorticity() );

      bg.setColor(Color.red);   // color for vectors
      bg.drawString("   Step No.="+IterNum, 4, 12);
      bg.setColor(Color.magenta);   // color for vectors
      bg.drawString("   /"+NSmax,90,12);
      bg.drawString("  dt="+dt,160,12);
      bg.setColor(Color.blue);	
      bg.drawString("iter="+NoP, 345,12);

      g.drawImage(drw.getBuffer(), 0, 0, drw);  // put the image on screen 
   }

   private void oneDimU(){
      int i, j, ip, jp, im, k; 
      float uin=0.0f, Qrat, Q1=0.0f;
      if(flow != "bend" && graph=="re_vectors")  uin=1.0f;
      for(i = 0; i<mx; i++){     ip=i+1; im=i-1;
         for(j = 0; j<my; j++){  jp=j+1; 
            k=my*i+j;
            up[k]=(float)(0.5*(u[ip][j]+u[ip][jp]))-uin;
            vp[k]=(float)(0.5*(v[i][jp]+v[ip][jp]));
            if(j==0 && im>0 ){
               if(BClow[i]=="wall" && BClow[im]=="periodic") vp[k]=0.0f;
               if(BClow[i]=="periodic" && BClow[im]=="wall"){up[k]=0.0f; vp[k]=0.0f;}
            }
         }
      }
   }
   public float[] getUarray(){ return up; }
   public float[] getVarray(){ return vp; }
   
   private float [] oneDimP(){
      int i, j, ip, jp, k;
      float ppp;  
      float pp [] = new float[kmax]; 
      for(i = 0; i<mx; i++){  ip=i+1;
         for(j=0; j<my; j++){
            k=my*i+j;  jp=j+1;
            pp[k]=(float)(0.25*(P[i][j]+P[i][jp]+P[ip][j]+P[ip][jp]));
         }
      }
      /*if(gr.getFlow()=="cylinder"){
         i=gr.mxL1*my;  ip=gr.mxL2*my; 
         ppp=0.5f*(pp[i]+pp[ip]);
         pp[i]=ppp;     pp[ip]=ppp;
      }*/
      return pp;
   }	
   private float [] oneDimT(){
      int i, j, ip, jp, k;
      // float ppp;  
      float pp [] = new float[kmax]; 
      for(i = 0; i<mx; i++){  ip=i+1;
         for(j=0; j<my; j++){
            k=my*i+j;  jp=j+1;
            pp[k]=(float)(0.25*(T[i][j]+T[i][jp]+T[ip][j]+T[ip][jp]));
         }
      }
      return pp;
   }

   //   calculate vorticity at u-node 
   private float[] vorticity(){  
      int i, j, ip, jp, jp1, k;
      float u0,u1,u2,u3; 
      float omeg[] = new float [kmax]; 
      float vor[][] = new float[mx][mym];
      for(i=0; i<mx; i++){
         for(j=0; j<mym; j++){ ip=i+1; jp=j+1; jp1=jp+1;
            u0=(float)(-me.XIxW[i][j]*(v[i][jp]+v[i][jp1]-v[ip][jp]-v[ip][jp1]));
            u1=(float)(-me.ETAxW[i][j]*(v[i][jp]-v[i][jp1]+v[ip][jp]-v[ip][jp1]));
            u2=(float)(me.XIyW[i][j]*(u[i][jp]-u[ip+1][jp]));
            u3=(float)(me.ETAyW[i][j]*(u[ip][j]-u[ip][jp1]));
            vor[i][j]=0.5f*(u0+u1+u2+u3);
         }
      }
      for(i=0; i<mx; i++){
         for(j=0; j<my; j++){
            k=my*i+j;
            if(j==0){
               if(BClow[i]=="wall") omeg[k]=0.0f; 
               else omeg[k]=0.5f*(vor[i][0]+vor[mxm-i][0]);
            }else if(j==mym){
               if(BCupp[i]=="wall") omeg[k]=0.0f;
               else omeg[k]=vor[i][mym-1];
            }else   omeg[k]=0.5f*(vor[i][j-1]+vor[i][j]);
               if(j==0 && i>1 ){
                  if(BClow[i]=="periodic" && BClow[i-1]=="wall") omeg[k]=0.0f;
               }
            }
      }
      return omeg;
   } // end of vorticity
   
   //========================================================================
   // substitution of initial velocity on staggered node pints
   //　 Install initail velocity on ficticious grid points 
   private void S_initBC(){
      int    i, j, ip, jp;  
      //   lower side boundary (wall or periodic )
      for(i = 0; i<mxm; i++){  
         ip=i+1;              v[ip][0] = v[ip][2]; 
         if(BClow[i]=="wall") u[ip][0]=-u[ip][1]; else u[ip][0]=u[ip][1];
      }
      if(BClow[mxm]=="wall")  u[mx][0]=-u[mx][1]; else u[mx][0]=u[mx][1];
      //   upper side boundary (wall, slip or inflow )
      for(i = 0; i<mxm; i++){  
         ip=i+1;              v[ip][my1]=v[ip][mym]; 
         if(BCupp[i]=="wall") u[ip][my]=-u[ip][mym]; else u[ip][my]=u[ip][mym];   
      }
      if(BCupp[mxm]=="wall")  u[mx][my]=-u[mx][mym]; else u[mx][my]=u[mx][mym];
      //   left side boundary  (inflow or outflow)
      for(j=0; j<mym; j++){   
         jp=j+1;   u[0][jp]=u[1][jp];  v[0][jp]=v[1][jp];   
      }
      v[0][0]=v[0][1];  v[0][my]=v[0][mym];
      //   right side boundary (outflow)
      for(j = 0; j<mym; j++){   
         jp=j+1;   u[mx1][jp]=u[mx][jp]; v[mx][jp]=v[mxm][jp];  
      }
      v[mx][my]=v[mxm][my];
   }
   
   //========================================================================
   // calculation of Energy equation
   private void calc_energy(){
      int    im, jm, ip, jp, ip1, jp1;
      double Uw, Ue, Vs, Vn, t0, t1, t2, t3, t4, t5, t6, t7, t8;
      double tw, te, tn, ts, aa=0.0, a1, a2, bb, Lap, Jp;
      alph=invRe/Data.Pr;

      BCforT();
      for(int i = 1; i<mx; i++){     im=i-1; ip=i+1; ip1=ip+1;
         for(int j = 1; j<my; j++){  jm=j-1; jp=j+1; jp1=jp+1;
            Uw=U[im][jm];  Ue=U[i][jm];  Vs=V[im][jm];  Vn=V[im][j];
            t1=T[im][jm];  t2=T[im][j];  t3=T[im][jp]; 
            t4=T[i][jm];   t0=T[i][j];   t5=T[i][jp]; 
            t6=T[ip][jm];  t7=T[ip][j];  t8=T[ip][jp];
            if(scheme == "donor-cell"){
               //  ** donor cell sheme (upstream scheme) **
               aa=0.5*(Ue*(t0+t7)+Math.abs(Ue)*(t0-t7)-Uw*(t2+t0)-Math.abs(Uw)*(t2-t0)
                  +Vn*(t0+t5)+Math.abs(Vn)*(t0-t5)-Vs*(t4+t0)-Math.abs(Vs)*(t4-t0));
             /* }else if(scheme == "first-order"){ // first order upstream scheme
                Uw=U[i][j]; Ue=Math.abs(Uw);
                Vs=0.25*(V[im][j]+V[im][jp]+V[i][j]+V[i][jp]); 
                Vn=Math.abs(Vs);
                aa=0.5*(Uw*(u7-u2)-Ue*(u2-2.0*u0+u7)+Vs*(u5-u4)-Vn*(u4-2.0*u0+u5)); */
            }else if(scheme == "QUICK"){    // QUICK scheme
               if(i==1)   tw=T[0][j];  else tw=T[im-1][j]; 
               if(i==mxm) te=T[mx][j]; else te=T[ip1][j]; 
               aa= Ue*(-t2+9.0*(t0+t7)-te)+Math.abs(Ue)*(-t2+3.0*(t0-t7)+te)
                  -Uw*(-tw+9.0*(t2+t0)-t7)-Math.abs(Uw)*(-tw+3.0*(t2-t0)+t7);
               if(j==1){
                  tn=T[i][jp1];   a2=0.0; 
                  a1=Vn*(-t4+9.0*(t0+t5)-tn)+Math.abs(Vn)*(-t4+3.0*(t0-t5)+tn);
               }else if(j==mym){
                  ts=T[i][jm-1];  a1=0.0;
                  a2=Vs*(-ts+9.0*(t4+t0)-t5)+Math.abs(Vs)*(-ts+3.0*(t4-t0)+t5); 
               }else{
                  ts=T[i][jm-1];  tn=T[i][jp1];
                  a1=Vn*(-t4+9.0*(t0+t5)-tn)+Math.abs(Vn)*(-t4+3.0*(t0-t5)+tn);
                  a2=Vs*(-ts+9.0*(t4+t0)-t5)+Math.abs(Vs)*(-ts+3.0*(t4-t0)+t5);
               }
               aa=(aa+a1-a2)/16.0;		   
            }
            // diffusion term
            Lap=me.D1P[im][jm]*t1+me.D2P[im][jm]*t2+me.D3P[im][jm]*t3
                +me.D4P[im][jm]*t4+me.D0P[im][jm]*t0+me.D5P[im][jm]*t5
                +me.D6P[im][jm]*t6+me.D7P[im][jm]*t7+me.D8P[im][jm]*t8;
            Jp=0.25*(me.JW[im][jm]+me.JW[i][jm]+me.JS[im][jm]+me.JS[im][j]);
            bb=-(aa-alph*Lap)/Jp;
            T[i][j]=t0+dt*(3.0*bb-Gt[i][j])*0.5;
            Gt[i][j]=bb;
         }
      }  
   }
   //  Update temperature on boundaries  
   private void BCforT(){
      int    i, j, ip, jp;
      // left sides boundary (inflow or outflow )
      for(j = 0; j<mym; j++){  
         jp=j+1;  T[0][jp]=T[1][jp];
         // if(BClef[j] == "in")  P[0][jp]=P[1][jp]; // p0;  
         // else                  P[0][jp]=P[1][jp];
      }
      // right sides boundary ( outflow )
      for(j = 0; j<mym; j++){ 
         jp=j+1;  T[mx][jp]=T[mxm][jp];
      }
      // lower side boundary ( wall or periodic ) 
      for(i = 0; i<mxm; i++){   ip=i+1;
         // if(BClow[i] == "wall") Tnp[ip][0]=Tnp[ip][1]; 
         //else                   Tnp[ip][0]=Tnp[mx-ip][1];
         T[ip][0]=Tsouth;
      }
      // upper side boundary  ( wall, in or slip )
      for(i = 0; i<mxm; i++){  ip=i+1;   
         // if(BCupp[i]=="in")  Tnp[ip][my]=Tnp[ip][mym];  // p0; 
         // else  
         T[ip][my]=0.0; // Tnp[ip][mym];  
      }
      // if(BClef[0]=="in")   Tnp[0][0]=Tnp[0][1];  else Tnp[0][0]=Tnp[mx][1];
      //if(BClef[mym]=="in") Tnp[0][my]=Tnp[0][mym]; else Tnp[0][my]=Tnp[0][mym];
      //if(BClow[mxm]=="wall") Tnp[mx][0]=Tnp[mx][1]; else Tnp[mx][0]=Tnp[0][1];
      //Tnp[mx][my]=Tnp[mx][mym];
      T[0][0]=T[1][0];    T[0][my]=T[1][my];
      T[mx][0]=T[mxm][0]; T[mx][my]=T[mxm][my];
   } 
} 