/************************************************************************/
/*   Graphics Routines of BetaFlow                                      */
/*          　DrawCanvas class:                          　         　　*/
/*                                 Last update: Novembe 9, 2001.        */
/*          All rights reserved, Copyright (C) 2001, Kiyoshi Minemura   */
/************************************************************************/
import java.awt.*;

//  ====== Graphics Output ============================================
public class DrawCanvas extends Canvas {   
   public int cWidth, cHeight;  // drawing range(picel) given by html 
   public int step=0;           // Computational step
   public int trans_x, trans_y; // ｘ, y translation of Scrollbars
   private Grid gr;   private Data da;   private Metric me; 
   private Solver so; 
   private Scrollbar scrH, scrV;   
   private float sx_t0, sx_t1, sy_t0, sy_t1; // 
   private int   mx, my, kmax, nmax;
   private int   xd, yd;                 // number of scale 　　　　　　　　　　　　
   private float xrange [];
   private float yrange [];   // Range of graph
   private float xp[], yp[];  //  up[], vp[], pp[];
   private int E[][]; 
   private boolean beta_flag = true;
   
   //   constructor of DrawCanvas
   public DrawCanvas(Scrollbar scrH, Scrollbar scrV, Grid gr, Data da, 
                     Metric me, Solver so){
      this.scrH=scrH;  this.scrV=scrV;  this.gr=gr; this.da=da; 
      this.me=me;      this.so=so;
   }
  
   public void init_P(Graphics g){
      setBackground(Color.white); setForeground(Color.black);
      cWidth=getSize().width;     cHeight=getSize().height;	
	  
      g.setColor(Color.blue); 
      g.setFont( new Font("TimesRoman",Font.ITALIC+Font.BOLD,24)); //Font.ITALIC,24));
      g.drawString("Welcome to BetaFlow !",140,50);
      g.setFont( new Font("TimesRoman",Font.PLAIN,12));
      g.drawString("With this BetaFlow you can simulate three kindes of flow by selecting SAMPLE Choice;", 100,70);
      g.drawString("[ duct ] = inlet pipe flow suddenly started in a straight duct,", 150,90);
      g.drawString("[ bend ] = flow suddenly started in a curved pipe,", 150, 110);
      g.drawString("[ cylinder ] = flow suddenly started around a cylinder.",150,130);
      g.drawString("After selecting the choice, push GRID-button, then you will find the grid to be calculated.",120,150);
      g.drawString("Push INITIAL V-button, then you will see the initial flow velocity distribution.",120,170);
      g.drawString("After these operations, the simulation will be started by pushing FLOW-button.",120,190);
      g.drawString("Before and in the middle of the simulation, you can change such calculating conditions",120,210);
      g.drawString("as Reynolds number, upstream scheme, graph, scale and Scrollbars.",150,230);
	  
     g.drawString("After starting the flow calculation, some numerical conditions will be indicated in the top of the figure;",100,260);
     g.drawString("Step No. means the number of dimensionless time step, dt is the incremental time to calculate, ",120,280);
     g.drawString("Qout/Qin the ratio of outlet flow rate to inlet one, from which the numerical accuracy can be confirmed,",120,300);
     g.drawString("and iter the iteration number for solving the Poisson equation with SOR method.",120,320);
     g.drawString("The contour lines of pressure and vorticity indicate the levels not of absolute value but of relative one.",120,340);

     g.drawString("BetaFlow calculates the Navier-Stokes equation for incompressible fluid flow",100,370);
     g.drawString("by discretizing it on a staggard arrangement of variables for a general coordinate system",120,390);
     g.drawString("based on a finite difference method and implemented by using a fractional step method.",120,410);
	  
     g.drawString("BetaFlow is intended to provide a simulator for the students on purpose to educate the fluid dynamics.",100,440);

     g.setColor(Color.green);
     g.drawString("All rights reserved,  Copyright (C) 2001,   Kiyoshi Minemura (Nagoya University, Japan)",150,470);
   }
   
   public void paint(Graphics g){    //  main routine
      if(beta_flag == true) init_P(g);
      if(step == 0){
         scrH.setValue(0);  scrV.setValue(0);  // initialize Scrollbar
      }
      // When "GRID" button is pressed;
      if(step == 1){
         beta_flag = false;
         g.setColor(Color.white);
         g.fillRect(0, 0, cWidth, cHeight);
         grid_plot(g);
      }
      //  When "INIT" button is pressed;
      if(step == 2){
         Veloc_mesh(g);  VelocIni(g);  // Draw initial velocity distribution
      }
      // When "STOP" button is presseed;
      if(step == 4)    so.drawStopImage();	  
   }  // end of paint	  
	  
   public void callScale(boolean scale){
      if(scale==true){  xrange=gr.getXrange(); yrange=gr.getYrange();
      }else{            xrange=da.getXrange(); yrange=da.getYrange();
      }
      init_tr();   //  Recalculate coordinate transformation coefficients
   }

   public void setScrollbar(Scrollbar scrH, Scrollbar scrV){
      this.scrH=scrH; this.scrV=scrV;
   }

   // ===== Graphics of pressure distributions ===============================
   public void DrawPressure(float pp[]){
      float zmin, zmax, dz;

      zmin=pp[0]; zmax=pp[0];
      for(int k=1; k<kmax; k++){
         if(zmax < pp[k]) zmax=pp[k];
         if(zmin > pp[k]) zmin=pp[k];
      }
      dz=(zmax-zmin)/10.0f;
  
      S_contour(pp, dz, zmin);
      // so.bg.setColor(Color.black);
      // so.bg.drawString(" zmax="+zmax,10,40);
      // so.bg.drawString("   zmin="+zmin,200,40);
   }
   
   //   Method to draw the contour lines:  (11 equally divided)
   public void S_contour(float z[], float dz, float zmin){
      // int   nmax=2*(mx-1)*(my-1);
      int    jmax, jmin, jmid=0, L, nt, k, ip, mxm=mx-1;
      float xs, ys, xe = 0.0f, ye=0.0f, xp1, yp1;
      float zf, x0, y0, z0, x1, y1, z1, x2, y2, z2, zf0, t;
      Graphics bg = so.bg;
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
               bg.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye));
            }else ip=0;
         }
      }
      //   draw boundary shape (wall):
      bg.setColor(Color.lightGray);
      for(int i=0; i<mxm; i++){
         if(da.BCupp[i]=="wall") { 
            k = my*i+my-1;  xs=xp[k]; ys=yp[k];
            k += my;        xe=xp[k]; ye=yp[k];
            bg.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye));
         }
      }
      for(int i=0; i<mxm; i++){
         if(da.BClow[i]=="wall") { 
            k = my*i;   xs=xp[k]; ys=yp[k];
            k += my;    xe=xp[k]; ye=yp[k];
            bg.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye));
         }
      }
      // k=0; xs=xp[k]; ys=yp[k]; k=my*(mx/2); xe=xp[k]; ye=yp[k];
      //  bg.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye));
   } //end of S_contour
   
   // ===== Output method of velocity vectors on a buffer image  ============
   public void S_veloc(float up[], float vp[]){
      float vScale=0.2f, yal=0.25f;   // length scales of vector and arrow
      double rad15=15.0*Math.PI/180.0;  // angle of arrow
      float ya[][]={{0.8f,0.54f},{0.8f,-0.54f}};  // elementary scale of arrow
      float th, x0, y0, x1, y1, vleng, costh, sinth;
      float sin15=(float)Math.sin(rad15), cos15=(float)Math.cos(rad15);
      float xs, xe, ys, ye;
      int    j0 , k;
      Graphics bg = so.bg;

      // draw boundary shape (wall):
      bg.setColor(Color.lightGray);
      for(int i=0; i<mx-1; i++){
         if(da.BCupp[i]=="wall") { 
            k = my*i+my-1;   xs=xp[k]; ys=yp[k];
            k += my;         xe=xp[k]; ye=yp[k];
            bg.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye));
         }
      }
      for(int i=0; i<mx-1; i++){
         if(da.BClow[i]=="wall") { 
            k = my*i;    xs=xp[k]; ys=yp[k];
            k += my;     xe=xp[k]; ye=yp[k];
            bg.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye));
         }
      }
  
      bg.setColor(Color.blue);    // draw vector by blue color
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
               bg.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye)); 
            }
         }
      }
   }// end of S_veloc
   
   // ===== Drawing method of initial velocity vectors ======================
   private void VelocIni(Graphics g){
      float vScale=0.25f, yal=0.25f;    // scales of vector and arrow
      double rad15=15.0*Math.PI/180.0;   
      float ya[][]={{0.8f,0.54f},{0.8f,-0.54f}}; 
      float th, x0, y0, x1, y1, vleng, costh, sinth;
      float sin15=(float)Math.sin(rad15), cos15=(float)Math.cos(rad15);
      float xs, xe, ys, ye;
      int k;
      so.up = da.getUarray();  so.vp = da.getVarray();

      g.setColor(Color.blue);
      for(int i=0; i<mx; i++){    
         for(int j=0; j<my; j++){
            k=my*i+j; x0=xp[k]; y0=yp[k];  // xp=gr.x[i][j]; yp=gr.y[i][j];
            vleng=(float)Math.sqrt(so.up[k]*so.up[k]+so.vp[k]*so.vp[k]); 
            if(vleng >= 0.05f){
               costh=so.up[k]/vleng;  sinth=so.vp[k]/vleng; 
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
               g.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye));
            }
         }
      }
   }// end of S_veloc

   // draw the mesh selected 
   private void Veloc_mesh(Graphics g){
      g.setColor(Color.lightGray);
      int i, j, k, k1; // , mx=gr.getMx(), my=gr.getMy();
      float xs, ys, xe, ye; 
      for(j=0; j<my; j++){
         for(i=0; i<mx-1; i++){
            k = my*i+j; k1 = k+my; 
            xs=xp[k]; ys=yp[k]; xe=xp[k1]; ye=yp[k1];  
            g.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye));  
         }
      }
      for(i=0; i<mx; i++){
         for(j=0; j<my-1; j++){
            k = my*i+j;  k1 = k+1;
            xs=xp[k]; ys=yp[k]; xe=xp[k1]; ye=yp[k1];
            g.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye));  
         }
      }
   }
	
   //  ======= Method for numbering and drawing the mesh ================= 
   private void grid_plot(Graphics g){
      float xs,xe,ys,ye;
      mx = gr.getMx(); my = gr.getMy(); kmax=mx*my; nmax=2*(mx-1)*(my-1);
      int i, j, k, k0, k1, k2;
      xp = new float [kmax];   yp = new float [kmax];
      E = new int [nmax][3];
	  
      xrange=gr.getXrange(); yrange=gr.getYrange();
      init_tr();
      xp = gr.getXarray();  yp = gr.getYarray();

      // number all the sides of triangular meshes 
      k=-1;
      for(i=0; i<mx-1; i++){
         for(j=0; j<my-1; j++){
            k++; E[k][0]=my*i+j; E[k][1]=my*i+j+1; E[k][2]=my*(i+1)+j;
            k++; E[k][0]=my*i+j+1; E[k][1]=my*(i+1)+j; E[k][2]=my*(i+1)+j+1;
         }
      }
      // draw horizontal and vertical mesh lines
      for(i=0; i<mx-1; i++){
         for(j=0; j<my-1; j++){
            k=2*((my-1)*i+j);  k0=E[k][0]; k1=E[k][1]; k2=E[k][2]; 
            xs=xp[k0]; ys=yp[k0]; xe=xp[k1]; ye=yp[k1];  
            g.setColor(Color.green);  // vertical line
            g.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye));
            xe=xp[k2]; ye=yp[k2]; 
            g.setColor(Color.blue);   // horizontal line
            g.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye));
            if(j == my-2){
               k++;  k0=E[k][0]; k2=E[k][2];
               xs=xp[k0]; ys=yp[k0]; xe=xp[k2]; ye=yp[k2];  
               g.setColor(Color.blue);
               g.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye));
            }
            if(i == mx-2){
               if(j != my-2) k++; 
               k1=E[k][1]; k2=E[k][2];
               xs=xp[k1]; ys=yp[k1];  xe=xp[k2]; ye=yp[k2];  
               g.setColor(Color.green);
               g.drawLine(xtr(xs), ytr(ys), xtr(xe), ytr(ye));
            }
         }
      }
   }
   

   //  Method for calculating coordinate transformation coeffcients
   public void init_tr() {
      float  sx0, sx1, sy0, sy1, s, ss; 
      float  sx_size=0.9f, sy_size=0.9f;     
      float  sWidth =sx_size*(float)cWidth; 
      float  sHeight=sy_size*(float)cHeight;
      float xmin=xrange[0], xmax=xrange[1], ymin=yrange[0], ymax=yrange[1];

      s=sWidth/sHeight*(xmax-xmin)/(ymax-ymin);
      if(s >= 1.0f) sx_size=sx_size/s;
      else          sy_size=sy_size/s;
      sx0=(float)(cWidth)*0.5f*(1.0f-sx_size); 
      sx1=cWidth-sx0 ;
      sy1=(float)(cHeight)*0.5f*(1.0f-sy_size); 
      sy0=cHeight-sy1 ;
      sx_t0=(float)((sx0*xmax-sx1*xmin)/(xmax-xmin)-trans_x);
      sx_t1=(float)((sx1-sx0)/(xmax-xmin));
      sy_t0=(float)((sy0*ymax-sy1*ymin)/(ymax-ymin)-trans_y); 
      sy_t1=(float)((sy1-sy0)/(ymax-ymin));
   }

   private int xtr(float x){    
      return (int)(sx_t1*x+sx_t0);
   }

   private int ytr(float y){    
      return(int)(sy_t1*y+sy_t0);
   }

   //  Drawing method for corrdinate axes
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
      g.drawLine(xtr(x0),ytr(y0),xtr(x1),ytr(y0));   
      g.drawLine(xtr(x0),ytr(y0),xtr(x0),ytr(y1));   
      g.drawString(x_t,xtr(x1)+5,ytr(y0));           
      g.drawString(y_t,xtr(x0),ytr(y1)-5);           
	  
      for(i=0; i<=xd; i++){   
          v = x0+(float)i*(x1-x0)/(float)xd;
          v=(float)((int)(v*1000.0f+0.5f))*0.001f;
          value = String.valueOf(v);
          strW = (getFontMetrics(font)).stringWidth(value);
          g.drawLine(xtr(v), ytr(y0)-m, xtr(v), ytr(y0)+m);
          g.drawString(value, xtr(v)-strW/2, ytr(y0)+m+2*strH);
      }
      for(i=0; i<=yd; i++){   
          float vv;
          v = y0+(float)i*(y1-y0)/(float)yd;
          v=(float)((int)(v*1000.0f+0.5f))*0.001f;
          value = String.valueOf(v);
          strW = (getFontMetrics(font)).stringWidth(value)+5;
          g.drawLine(xtr(x0)-m, ytr(v), xtr(x0)+m, ytr(v));
          g.drawString(value, xtr(x0)-strW, ytr(v)+strH/2);
      }
   }
	
   // 
   private void line(float x0, float y0, float x1, float y1, Color c){
      Graphics g=getGraphics();
		   
      Color current_color = g.getColor();
      g.setColor(c);
      g.drawLine(xtr(x0), ytr(y0), xtr(x1), ytr(y1));
      g.setColor(current_color);
   }
}