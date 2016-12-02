/**====================================================================
   Graphics Routines for simulating Benard natural heat convection                        　DrawCanvas.java                          　         　　
      All rights reserved, Copyright (C) 2001-3003, Kiyoshi Minemura
          ver.2.0    Last update: December 25, 2002. 
=====================================================================*/
import java.awt.*;

public class DrawCanvasB extends Canvas { 
   public Graphics g, bg;
   private Image buffer;
   private ViewPort  vp0, vp1, vp2, vp3;
   private int Width, Height;  // drawing range(pixel) given by html 
   private int step = 0;       // Computational step
   private float sx_t0, sx_t1, sy_t0, sy_t1; // 
   private int   mx, my, kmax, nmax;
   private float xp[], yp[], up[], vp[];
   private int E[][]; 
   private String BCupp[],BClow[];
   private float range[], xmin,xmax,ymin,ymax;
   
   // access methods
   public void setStep( int i){ step = i; }
   public void setRange( double ra[] ){
      xmin = (float)ra[0];  xmax = (float)ra[1]; 
      ymin = (float)ra[2];  ymax = (float)ra[3];
      range = new float[]{ xmin, xmax, ymin, ymax };
      // init_tr();
   }
   public void setGrid( int co[] ){
      mx = co[0];  my = co[1];	
      kmax=mx*my; nmax=2*(mx-1)*(my-1);
      numbering();
   }
   public void setXpYpArray( float x[], float y[] ){
      xp = x;   yp = y;
   }
   public void setBCarray( String up[], String low[] ){
      BCupp = up;   BClow = low; 
   }
   public void setVarray( float u[], float v[]){
      up = u;  vp = v;
   }
   
   //  =====  main routine  =====
   public void paint(Graphics g){
      if(step == 0) initBenard(); // When "GRID" button is pressed;
      if(step == 1) initGrid();   // When "INIT" button is pressed;
      if(step == 2) initVeloc();  // When "START" button is pressed;
      if(step == 3) transView();  // When "STOP" button is presseed;
      if(step == 4) g.drawImage(buffer, 0,0, this); 	  
   }
   public void deleteBuffer(){
      bg.clearRect( 0, 0, Width, Height );
   }

   public void initBenard(){
      setBackground( Color.white ); //  setForeground( Color.black );
      Width = getSize().width;    Height = getSize().height;	
      buffer = createImage( Width, Height );
      bg = buffer.getGraphics();
      g = getGraphics();
  
      g.setColor(Color.blue); 
      g.setFont( new Font("TimesRoman",Font.ITALIC+Font.BOLD,24)); //Font.ITALIC,24));
      g.drawString("Welcome to BetaFlow for natural heat convection!",50,60);
      g.setFont( new Font("TimesRoman",Font.PLAIN,12));
      g.drawString("With this BetaFlow you can simulate a natural heat convection, Benard convection.", 60,130);
      // g.drawString("[ duct ] = inlet pipe flow suddenly started in a straight duct,", 110,90);
      // g.drawString("[ bend ] = flow suddenly started in a curved pipe,", 110, 110);
      // g.drawString("[ cylinder ] = flow suddenly started around a cylinder.",110,130);
      g.drawString("By pushing GRID-button, you will find the grid to be calculated.",80,150);
      g.drawString("Push INITIAL V-button, then you will see the initial flow velocity distribution.",80,170);
      g.drawString("After these operations, the simulation will be started by pushing START-button.",80,190);
      g.drawString("Before the start, you can select such calculating conditions as Rayleigh-, Prandtle-, Reynolds-numbers.",80,210);
//      g.drawString("After the start, you can change such conditions as upstream scheme, graph and Scrollbars.",80,230);
	  
      g.drawString("After starting the flow calculation, some numerical conditions will be indicated in the top of the figure;",60,280);
      g.drawString("Step No. means the number of dimensionless time step, dt is the incremental time to calculate, ",80,300);
      // g.drawString("Qout/Qin the ratio of outlet flow rate to inlet one, from which the numerical accuracy can be confirmed,",80,300);
      g.drawString("and iter the iteration number for solving the Poisson equation with SOR method.",80,320);
      g.drawString("The contour lines of temperature and vorticity indicate the levels not of absolute value but of relative one.",80,340);

      g.drawString("BetaFlow calculates the Navier-Stokes equation for incompressible fluid flow",60,370);
      g.drawString("by discretizing it on a staggard arrangement of variables for a general coordinate system",80,390);
      g.drawString("based on a finite difference method and implemented by using a fractional step method.",80,410);
	  
      g.drawString("BetaFlow is intended to provide a simulator for the students on purpose to educate the fluid dynamics.",60,440);

      g.setColor(Color.green);
      g.drawString("All rights reserved,  Copyright (C) 2001-2003,   Kiyoshi MINEMURA (Nagoya University, Japan)",110,470);
   }
   
   /* public void callScale(boolean scale){
       if(scale==true){
         xrange=gr.getXrange(); yrange=gr.getYrange();
       }else{
         xrange=da.getXrange(); yrange=da.getYrange();
       }
       init_tr();	//  Recalculate coordinate transformation coefficients
     }*/
    /* public void setScrollbar(Scrollbar scrH, Scrollbar scrV){
         this.scrH=scrH; this.scrV=scrV;
    }*/

   // display mathod for initial velocity distribution
   public void initVeloc(){
      deleteBuffer(); 
      drawGridV(  );  
      drawVelocI();   // Draw initial velocity distribution
      g.setColor( Color. blue ); 
      g.drawString("Initial velocity distributions (=0)",80,100);
   }
   // === method for drawing the results ===
   
   public void drawBenard(float up[], float vp[], float p[], float vor[]){
      // temperature
      bg.setColor(Color.blue);
      bg.drawString("Contour map for tempature", 40,60);
      drawTemperature( p );  
   
      // velocity vectors
      bg.setColor(Color.blue);
      bg.drawString("Velocity vectors", 50, 210);
      drawVectors( up, vp );
   
      // vorticity
      bg.drawString("Contour map for vorticity",50,350);
      drawVorticity( vor );
   }
   public void transView(){
      vp1 = new ViewPort( Width, Height, range );
      vp1.transView( 50,310, true );
      vp2=new ViewPort( Width, Height, range );
      vp2.transView( 200,160, true );
      vp3=new ViewPort( Width, Height, range );
      vp3.transView( 340,20,true );
   }
   
   // ===== Graphics of pressure distributions ===============================
   public void drawTemperature( float pp[] ){
      float zmin, zmax, dz;

      zmin=pp[0]; zmax=pp[0];
      for(int k=1; k<kmax; k++){
         if(zmax < pp[k]) zmax=pp[k];
         if(zmin > pp[k]) zmin=pp[k];
      }
      dz = (zmax-zmin)/10.99f;
  
      S_contour(pp, dz, zmin);
      // so.bg.setColor(Color.black);
      // so.bg.drawString(" zmax="+zmax,10,40);
      // so.bg.drawString("   zmin="+zmin,200,40);
   }
   
   // numbering method for all the sides of triangular meshes 
   public void numbering(){
      E = new int [nmax][3];
      int k=-1, j;
      for(int i=0; i<mx-1; i++){
         for(j=0; j<my-1; j++){
            k++; E[k][0]=my*i+j; E[k][1]=my*i+j+1; E[k][2]=my*(i+1)+j;
            k++; E[k][0]=my*i+j+1; E[k][1]=my*(i+1)+j; E[k][2]=my*(i+1)+j+1;
         }
      }
   }
   
   //   Method to draw the contour lines:  (11 equally divided)
   public void S_contour(float z[], float dz, float zmin){
      // int   nmax=2*(mx-1)*(my-1);
      int    jmax, jmin, jmid=0, L, nt, k, ip, mxm=mx-1;
      float xs, ys, xe = 0.0f, ye=0.0f, xp1, yp1;
      float zf, x0, y0, z0, x1, y1, z1, x2, y2, z2, zf0, t;
      Color co[]={Color.black, Color.gray, Color.blue, Color.cyan, Color.green,
                  Color.pink, Color.red, Color.magenta, Color.orange, Color.yellow};
  
      for(k=0; k<nmax; k++){
         jmax = 0; z2 = (z[ E[k][0] ]-zmin); jmin = 0;  z0 = z2; 
         for(int m=1; m<=2; m++){
            zf = z[ E[k][m] ]-zmin; 
            if(zf > z2){      jmax = m;  z2 = zf;}
            else if(zf < z0){ jmin = m;  z0 = zf;}
         }
         nt=(int)(z2/dz)-(int)(z0/dz);
         for(int j=0; j<3; j++){ 
            if(j != jmax && j != jmin) jmid=j;
         }
         if(nt <= 0) continue;
         jmax = E[k][jmax]; jmin = E[k][jmin]; jmid = E[k][jmid];
         x2 = xp[jmax]; x0 = xp[jmin];  y2 = yp[jmax];  y0 = yp[jmin];
         z1 = z[jmid]-zmin;  x1 = xp[jmid];  y1 = yp[jmid];
         for(int m=1; m<=nt; m++){
            L = (int)(z0/dz)+m; 
            zf = (float)L*dz;   L = (int)((zf)/dz)-1;
            t = (zf-z0)/(z2-z0); ip=1;
            xs=(1.0f-t)*x0+t*x2;
            ys=(1.0f-t)*y0+t*y2;  
            if(zf > z1){
               t=(zf-z1)/(z2-z1);  if(t>1.0) t=1.0f;
               xe=(1.0f-t)*xp[jmid]+t*x2;
               ye=(1.0f-t)*yp[jmid]+t*y2;  ip++; // =ip+1;
            }else{
               t=(zf-z0)/(z1-z0); if(t>1.0) t=1.0f;
               xe=(1.0f-t)*x0+t*xp[jmid];
               ye=(1.0f-t)*y0+t*yp[jmid];  ip++;  // =ip+1;
            }
            if(ip == 2){
               bg.setColor(co[L]);
               bg.drawLine(vp1.xtr(xs), vp1.ytr(ys), vp1.xtr(xe), vp1.ytr(ye));
            }else ip=0;
         }
      }
      //   draw boundary shape (wall):
      bg.setColor(Color.lightGray);
      for(int i=0; i<mxm; i++){
         if(BCupp[i]=="wall") { 
            k = my*i+my-1;  xs=xp[k]; ys=yp[k];
            k += my;        xe=xp[k]; ye=yp[k];
            bg.drawLine(vp1.xtr(xs), vp1.ytr(ys), vp1.xtr(xe), vp1.ytr(ye));
         }
      }	  
      for(int i=0; i<mxm; i++){
         if(BClow[i]=="wall") { 
            k = my*i;   xs=xp[k]; ys=yp[k];
            k += my;    xe=xp[k]; ye=yp[k];
            bg.drawLine(vp1.xtr(xs), vp1.ytr(ys), vp1.xtr(xe), vp1.ytr(ye));
         }
      }
   }
   
   // ===== output velocity vectors on a buffer image  =====================
   public void drawVectors(float up[], float vp[]){
      float vScale=0.2f, yal=0.25f;    // scales for vector and its arrow
      double rad15=15.0*Math.PI/180.0; // openning angle of arrow
      float ya[][]={{0.8f,0.54f},{0.8f,-0.54f}};  // basic form of arrow
      float th, x0, y0, x1, y1, vleng, costh, sinth;
      float sin15=(float)Math.sin(rad15), cos15=(float)Math.cos(rad15);
      float xs, xe, ys, ye;
      int    j0 , k;

      // draw boundary shape (wall):
      bg.setColor(Color.lightGray);
      for(int i=0; i<mx-1; i++){
         if(BCupp[i]=="wall") { 
            k = my*i+my-1;   xs=xp[k]; ys=yp[k];
            k += my;         xe=xp[k]; ye=yp[k];
            bg.drawLine(vp2.xtr(xs), vp2.ytr(ys), vp2.xtr(xe), vp2.ytr(ye));
         }
      }	  
      for(int i=0; i<mx-1; i++){
         if(BClow[i]=="wall") { 
            k = my*i;    xs=xp[k]; ys=yp[k];
            k += my;     xe=xp[k]; ye=yp[k];
            bg.drawLine(vp2.xtr(xs), vp2.ytr(ys), vp2.xtr(xe), vp2.ytr(ye));
         }
      }
  
      bg.setColor(Color.blue);    // draw vectors with blue color
      for(int i=0; i<mx; i++){// im=i-1; 
         j0=0; 
         for(int j=j0; j<my; j++){     
            k=my*i+j;  x0=xp[k]; y0=yp[k]; // x0=x[i][j]; y0=y[i][j];
            vleng=(float)Math.sqrt(up[k]*up[k]+vp[k]*vp[k]);
            if(vleng >= 0.01f){
               costh=up[k]/vleng;  sinth=vp[k]/vleng;
            }else{  
               sinth=0.0f; costh=0.0f; 
            }
            vleng=vleng*vScale; x1=x0+vleng*costh; y1=y0+vleng*sinth;
            for(k=1; k<=3; k++){
               xs=x0; ys=y0; xe=x1; ye=y1; 
               if(k==2){
                  xs=x1; ys=y1; 
                  xe=x1-vleng*yal*(costh*ya[0][0]-sinth*ya[0][1]);
                  ye=y1-vleng*yal*(sinth*ya[0][0]+costh*ya[0][1]);}
               else if(k==3){ 
                  xs=x1; ys=y1;
                  xe=x1-vleng*yal*(costh*ya[1][0]-sinth*ya[1][1]);
                  ye=y1-vleng*yal*(sinth*ya[1][0]+costh*ya[1][1]);}  
               bg.drawLine(vp2.xtr(xs), vp2.ytr(ys), vp2.xtr(xe), vp2.ytr(ye)); 
            }
         }
      }
   }
   
    public void drawVorticity(float pp[]){
      float zmin, zmax, dz;

      zmin=pp[0]; zmax=pp[0];
      for(int k=1; k<kmax; k++){
         if(zmax < pp[k]) zmax=pp[k];
         if(zmin > pp[k]) zmin=pp[k];
      }
      dz = (zmax-zmin)/10.99f;
 
      S_contourV(pp, dz, zmin);
      // so.bg.setColor(Color.black);
      // so.bg.drawString(" zmax="+zmax,10,40);
      // so.bg.drawString("   zmin="+zmin,200,40);
   }  
   //   Method to draw the contour lines:  (11 equally divided)
   public void S_contourV(float z[], float dz, float zmin){
      // int   nmax=2*(mx-1)*(my-1);
      int    jmax, jmin, jmid=0, L, nt, k, ip, mxm=mx-1;
      float xs, ys, xe = 0.0f, ye=0.0f, xp1, yp1;
      float zf, x0, y0, z0, x1, y1, z1, x2, y2, z2, zf0, t;
      Color co[]={Color.black, Color.gray, Color.blue, Color.cyan, Color.green,
                  Color.pink, Color.red, Color.magenta, Color.orange, Color.yellow};
  
      for(k=0; k<nmax; k++){
         jmax = 0; z2 = (z[ E[k][0] ]-zmin); jmin = 0;  z0 = z2; 
         for(int m=1; m<=2; m++){
            zf=z[ E[k][m] ]-zmin; 
            if(zf > z2){      jmax = m;  z2 = zf;}
            else if(zf < z0){ jmin = m;  z0 = zf;}
         }
         nt=(int)(z2/dz)-(int)(z0/dz);
         for(int j=0; j<3; j++){ 
            if(j != jmax && j != jmin) jmid=j;
         }
         if(nt <= 0) continue;
         jmax = E[k][jmax]; jmin = E[k][jmin]; jmid = E[k][jmid];
         x2 = xp[jmax]; x0 = xp[jmin];  y2 = yp[jmax];  y0 = yp[jmin];
         z1 = z[jmid]-zmin;  x1 = xp[jmid];  y1 = yp[jmid];
         for(int m=1; m<=nt; m++){
            L = (int)(z0/dz)+m; 
            zf = (float)L*dz;   L = (int)((zf)/dz)-1;
            t = (zf-z0)/(z2-z0); ip=1;
            xs=(1.0f-t)*x0+t*x2;
            ys=(1.0f-t)*y0+t*y2;  
            if(zf > z1){
               t=(zf-z1)/(z2-z1);  if(t>1.0) t=1.0f;
               xe=(1.0f-t)*xp[jmid]+t*x2;
               ye=(1.0f-t)*yp[jmid]+t*y2;  ip++; // =ip+1;
            }else{
               t=(zf-z0)/(z1-z0); if(t>1.0) t=1.0f;
               xe=(1.0f-t)*x0+t*xp[jmid];
               ye=(1.0f-t)*y0+t*yp[jmid];  ip++;  // =ip+1;
            }
            if(ip == 2){
               bg.setColor(co[L]);
               bg.drawLine(vp3.xtr(xs), vp3.ytr(ys), vp3.xtr(xe), vp3.ytr(ye));
            }else ip=0;
         }
      }
      //   draw boundary shape (wall):
      bg.setColor(Color.lightGray);
      for(int i=0; i<mxm; i++){
         if(BCupp[i]=="wall") { 
            k = my*i+my-1;  xs=xp[k]; ys=yp[k];
	    k += my;        xe=xp[k]; ye=yp[k];
            bg.drawLine(vp3.xtr(xs), vp3.ytr(ys), vp3.xtr(xe), vp3.ytr(ye));
         }
      }
      for(int i=0; i<mxm; i++){
         if(BClow[i]=="wall") { 
            k = my*i;   xs=xp[k]; ys=yp[k];
            k += my;    xe=xp[k]; ye=yp[k];
            bg.drawLine(vp3.xtr(xs), vp3.ytr(ys), vp3.xtr(xe), vp3.ytr(ye));
         }
      }
      // k=0; xs=xp[k]; ys=yp[k]; k=my*(mx/2); xe=xp[k]; ye=yp[k];
      //  bg.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye));
   } //en  
   // ===== drawing method of initial velocity vectors =================
   private void drawVelocI(){
      float vScale=0.25f, yal=0.25f;    // scales of vector and its arrow
      double rad15=15.0*Math.PI/180.0;  // openning angle of arrow
      float ya[][]={{0.8f,0.54f},{0.8f,-0.54f}};  // basic form of arrow
      float th, x0, y0, x1, y1, vleng, costh, sinth;
      float sin15=(float)Math.sin(rad15), cos15=(float)Math.cos(rad15);
      float xs, xe, ys, ye;

      g.setColor(Color.blue);
      for(int i=0; i<mx; i++){    
         for(int j=0; j<my; j++){
            int k=my*i+j; x0=xp[k]; y0=yp[k];  
            vleng=(float)Math.sqrt(up[k]*up[k]+vp[k]*vp[k]); 
            if(vleng >= 0.05f){
               costh=up[k]/vleng;  sinth=vp[k]/vleng; 
            }else{
               sinth=0.0f; costh=0.0f;
            }
            vleng=vleng*vScale; x1=x0+vleng*costh; y1=y0+vleng*sinth;
            for(k=1; k<=3; k++){
               xs=x0; ys=y0; xe=x1; ye=y1;  
               if(k==2){
                  xs=x1; ys=y1; 
                  xe=x1-vleng*yal*(costh*ya[0][0]-sinth*ya[0][1]);
                  ye=y1-vleng*yal*(sinth*ya[0][0]+costh*ya[0][1]);}
               else if(k==3){
                  xs=x1; ys=y1;
                  xe=x1-vleng*yal*(costh*ya[1][0]-sinth*ya[1][1]);
                  ye=y1-vleng*yal*(sinth*ya[1][0]+costh*ya[1][1]);}  
               g.drawLine(vp0.xtr(xs), vp0.ytr(ys), vp0.xtr(xe), vp0.ytr(ye));
             }
	  }
      }
   }
   // method for drawing the mesh selected before drawing initila velocities
   private void drawGridV(){
      g.setColor(Color.lightGray);
      int i, j, k, k1; // , mx=gr.getMx(), my=gr.getMy();
      float xs, ys, xe, ye; 
      for(j=0; j<my; j++){
         for(i=0; i<mx-1; i++){
            k = my*i+j; k1 = k+my; 
            xs=xp[k]; ys=yp[k]; xe=xp[k1]; ye=yp[k1];  
            g.drawLine(vp0.xtr(xs), vp0.ytr(ys), vp0.xtr(xe), vp0.ytr(ye));  
         }
      }
      for(i=0; i<mx; i++){
         for(j=0; j<my-1; j++){
            k = my*i+j;  k1 = k+1;   
            xs=xp[k]; ys=yp[k]; xe=xp[k1]; ye=yp[k1];
            g.drawLine(vp0.xtr(xs), vp0.ytr(ys), vp0.xtr(xe), vp0.ytr(ye));  
         }
      }
   }
	
   //  ======= Method for numbering and drawing the mesh ================= 
   //  ===== method for drawing grid employed ==========
   public void initGrid(){
      g.clearRect( 0, 0, Width, Height );
      // kmax=mx*my; nmax=2*(mx-1)*(my-1);
      vp0 = new ViewPort( Width, Height, range );
      vp0.transView( 0, 20, true ); // init_tr();
      g.setColor( Color. blue ); 
      g.drawString("Mesh for the computational domain",80,100);
      drawGrid();	  
   }
   // === method for drawing squer grid  
   public void drawGrid(){
      int i, j, k, k1;
      float xs, ys, xe, ye; 
      g.setColor( Color.blue ); // For horizontal lines 
      for(j=0; j<my; j++){
         for(i=0; i<mx-1; i++){
            k = my*i+j; k1 = k+my; 
            xs=xp[k]; ys=yp[k]; xe=xp[k1]; ye=yp[k1];  
            g.drawLine(vp0.xtr(xs), vp0.ytr(ys), vp0.xtr(xe), vp0.ytr(ye));  
         }
      }
      g.setColor(Color.green);  // For vetical lines
      for(i=0; i<mx; i++){
         for(j=0; j<my-1; j++){
            k = my*i+j;  k1 = k+1;   
            xs=xp[k]; ys=yp[k]; xe=xp[k1]; ye=yp[k1];
            g.drawLine(vp0.xtr(xs), vp0.ytr(ys), vp0.xtr(xe), vp0.ytr(ye));  
         }
      }
   }
	
   //  method for calculating coordinate transformation
   public void init_tr() {
      float  sx0, sx1, sy0, sy1, s=1.0f, ss; 
      float  sx_size=0.9f, sy_size=0.9f;  // ratio of drawing range in x,y
      float  sWidth =sx_size*(float)Width; 
      float  sHeight=sy_size*(float)Height;
      s=sWidth/sHeight*(xmax-xmin)/(ymax-ymin);
         if(s >= 1.0f) sx_size=sx_size/s; 
         else          sy_size=sy_size/s;  
       sx0=(float)(Width)*0.5f*(1.0f-sx_size); 
       sx1=Width-sx0 ;
       sy1=(float)(Height)*0.5f*(1.0f-sy_size); 
       sy0=Height-sy1 ;
       sx_t0=(float)((sx0*xmax-sx1*xmin)/(xmax-xmin));
       sx_t1=(float)((sx1-sx0)/(xmax-xmin));
       sy_t0=(float)((sy0*ymax-sy1*ymin)/(ymax-ymin)); 
       sy_t1=(float)((sy1-sy0)/(ymax-ymin));
   }

   private int xtr(float x){    // tansformation of x-coordinate
      return (int)(sx_t1*x+sx_t0);
   }

   private int ytr(float y){    // transformation of y-coordinate
      return(int)(sy_t1*y+sy_t0);
   }

   //  drawing method of coordinates
   private void axis(Graphics g,float x0,float y0,float x1,float y1,
                     int xd, int yd, String y_t, String x_t){ 
      float v;
      int    i, strW, strH, m=5;
      String value;
      Font font = new Font("TimesRoman", Font.BOLD, 10);
      strH = (getFontMetrics(font)).getLeading();
      strH -= (getFontMetrics(font)).getDescent();
      strH += (getFontMetrics(font)).getAscent();
      //	  y_t="y";
      g.setColor(Color.black);
      g.drawLine(xtr(x0),ytr(y0),xtr(x1),ytr(y0));   // x-axis
      g.drawLine(xtr(x0),ytr(y0),xtr(x0),ytr(y1));   // y-axis
      g.drawString(x_t,xtr(x1)+5,ytr(y0));           // name of x-axis
      g.drawString(y_t,xtr(x0),ytr(y1)-5);           // name of y-axis
	  
      for(i=0; i<=xd; i++){    // draw x-axis
          v = x0+(float)i*(x1-x0)/(float)xd;
          v=(float)((int)(v*1000.0f+0.5f))*0.001f;
          value = String.valueOf(v);
          strW = (getFontMetrics(font)).stringWidth(value);
          g.drawLine(xtr(v), ytr(y0)-m, xtr(v), ytr(y0)+m);
          g.drawString(value, xtr(v)-strW/2, ytr(y0)+m+2*strH);
      }
      for(i=0; i<=yd; i++){    // draw y-axis
          float vv;
          v = y0+(float)i*(y1-y0)/(float)yd;
          v=(float)((int)(v*1000.0f+0.5f))*0.001f;
          value = String.valueOf(v);
          strW = (getFontMetrics(font)).stringWidth(value)+5;
          g.drawLine(xtr(x0)-m, ytr(v), xtr(x0)+m, ytr(v));
          g.drawString(value, xtr(x0)-strW, ytr(v)+strH/2);
      }
   }
	
   // drawing method for connecting two points
   private void line(float x0, float y0, float x1, float y1, Color c){
      Color current_color = g.getColor();
      g.setColor(c);
      g.drawLine(xtr(x0), ytr(y0), xtr(x1), ytr(y1));
      g.setColor(current_color);
   }    
 
   public Image getBuffer(){ return buffer;}
}


// ============ View port transformation ===================================
class ViewPort{
   private float sx = 0.9f, sy = 1.0f; // ratio of drawing area to screen area 
   private int   Width, Height, Nxmin, Nymin; 
   private float xmin, xmax, ymin, ymax, rx, ry;
   
   // constructor
   public ViewPort( int width, int height, float range[] ){
      Width = width;  Height = height; 
      xmin = range[0];  xmax = range[1];
      ymin = range[2];  ymax = range[3]; 
   }

   public void transView( int Nup, int Nbottom, boolean aspect ){
      //  Nup,Nbottom = upside and bottom area of no drawing  
      //  aspect = true for keeping the aspect ratio of physical space.  
      float ap = (ymax-ymin)/(xmax-xmin); 
      float aw = (float)(Height-Nup-Nbottom)*sy/(float)Width/sx;
      rx = (float)Width*sx/(xmax-xmin); 
      ry = (float)(Height-Nup-Nbottom)*sy/(ymax-ymin);
      Nxmin = (int)((float)Width*(1.0f-sx))/2;  
      Nymin = Nbottom+(int)((float)(Height-Nup-Nbottom)*(1.0f-sy))/2;
      if(aspect == true){
         if(ap > aw){ 
            Nxmin += (int)((1.0f-ry/rx)*(float)Width*sx)/2; 
            rx = ry;
         }else{
            Nymin += (int)((1.0f-rx/ry)*(float)Height*sy)/2;
            ry = rx;
         }  
      }
   }
   //  === method for transformation of x-coodinate ===
   public int xtr( double x ){
      return (int)(rx*((float)x-xmin))+Nxmin;
   }
   //  === method for transformation of y-coordinate ===
   public int ytr( double y ){
      return Height-Nymin-(int)(ry*((float)y-ymin));
   }
   //  === over-load methods for the above two methods ===
   public int xtr( float x ){
      return (int)(rx*(x-xmin))+Nxmin;
   }
   //   
   public int ytr( float y ){
      return Height-Nymin-(int)(ry*(y-ymin));
   }
}
