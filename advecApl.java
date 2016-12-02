/**************************************************************************/
/*    Simulation of an Advection Equation　　　                           */
/*       advecApl.java                                                    */
/*                 All Rights Reserved, Copyright (C) 2001, K. MINEMURA   */
/*       　　　　      Last updated by K. Minemura on November 6, 2001.   */
/**************************************************************************/
import java.applet.Applet;
import java.awt.*;
import java.awt.event.*;

//    Java Applet for Advection Equation
public class advecApl extends Applet implements
	         ItemListener,ActionListener,Runnable{
   // declaration of instance variables for graphics
   private Image  buffer;    //  buffer image
   private Graphics  bg;     //  graphic object of buffer image
   private Thread simulator; //  thread for simulation control
   private Button bt_start,bt_stop,bt_reset;  // button objects
   private int    wx, wy;   //  width and height of applet (pixel unit)
   private double xmin,xmax,ymin,ymax; // physical drawing range for applet
   private double sx_t0,sx_t1,sy_t0,sy_t1; // coordinate transformation
   private Choice dt_Choi, sche_Choi;   // choice objects
   private TextField Cour_t;            // text field object

   // declaration of instance variables for computation
   private String scheme;               // upstream scheme
   private int  delay=50;               // sleep period for animation
   // private boolean sml_start=false;     // flag for animation
   private double f[],f0[],fb[],fb0[],fc[],fc0[],g[],g0[],fa[]; // velocities
   private double x[];      // grid coordinate
   private int    Nx=100;   // number of grid
   private double Cour=0.5; // Courant number
   private double dx;       // grid width
   private double dt=0.025; // time step to calculate
   private double xa;       // x-coordinate of wave front
   private double t=0.0;    // elapsed time
   private double u0 = 1.0; // velocity

   /*****  Declaration of methods for applet　*************************/
   //  **　Initializing method for the applet
   public void init(){
      add(new Label("   dt="));   //  paste of label "dt="
      dt_Choi = new Choice();     //  generation of dt-choice
      dt_Choi.addItem("0.01");    //  setting up of dt selection
      dt_Choi.addItem("0.025"); dt_Choi.addItem("0.045");
      dt_Choi.addItem("0.05");  dt_Choi.addItem("0.0505");
      dt_Choi.select(1);      // institution of second initial value of dt
      add(dt_Choi);           // paste of dt-choice

      add(new Label("  Cour No.=")); // paste of Courant number label
      String Cour_No=new Double(Cour).toString();
                    //　convertion of Courant number to character
      Cour_t = new TextField();  //generation of text-field for Courant No.
      Cour_t.setText(Cour_No);   Cour_t.setEditable(false);
      add(Cour_t);                // paste of Courant number as read only
      add(new Label("  Scheme="));  // paste of scheme label
      sche_Choi = new Choice();     // generation of scheme choice
      sche_Choi.addItem("Up-stream");
      sche_Choi.addItem("LaxWendrof"); sche_Choi.addItem("CIP");
      sche_Choi.select(0);   add(sche_Choi); // paste of scheme choice

      setBackground(Color.white ); // setting up of background color(defort)
      setForeground (Color.black );    // setting up of drawing color as black

      bt_start= new Button("Start"); add(bt_start);
      bt_stop = new Button("Stop");  add(bt_stop);
      bt_reset= new Button("Reset"); add(bt_reset);// paste of reset button

      // setting up of event interface
      dt_Choi.addItemListener(this);
      bt_start.addActionListener(this);
      bt_stop.addActionListener(this);
      bt_reset.addActionListener(this);
      sche_Choi.addItemListener(this);

      // simulator = new Thread(this);  // generation of thread object
      wx = getSize().width;          // acquisition of applet width
      wy = getSize().height;         // acquisition of applet height
      buffer = createImage(wx, wy);  // generation of buffer image
      bg     = buffer.getGraphics(); // graphics object for buffer image
            //  setting up of display range for physical coordinate
      getRange(0.0, -0.5, 5.0, 1.5);
      x = new double [Nx+1];         // arrangement for x coordinate
      dx=(xmax-xmin)/(double)Nx;     // incremental distance in x-direction
      for(int i=0; i<=Nx; i++){ x[i]=dx*(double)i; }
      trans();         // calculate coordiante transformation for display
      initFlow();      // set up initial value
      drawFlow();      // draw initial value
   }
   // **  Manegement method for the buttons.
   public void actionPerformed(ActionEvent ev){
      if(ev.getSource() == bt_start){  // when start-button is pressed,
         // sml_start = true;
         if(simulator == null ) simulator = new Thread(this);     // generation of thread.
         simulator.start();
      }else if(ev.getSource() == bt_stop){ // when stop-button is pressed,
         // sml_start = false;
         simulator = null;
      }else if(ev.getSource() == bt_reset){ // when reset-button is pressed,
         // sml_start = false;
         simulator = null;
         initFlow();          // recalculate the initial value
      }
   }
   // **  Management method for the choices.
   public void itemStateChanged(ItemEvent ev){
      String str;
      if(ev.getSource() == dt_Choi){ // When dt-choice is selected,
         str = dt_Choi.getSelectedItem();        // get character
         dt = Double.valueOf(str).doubleValue(); // convert to number.
      }else if(ev.getSource() == sche_Choi){ // When scheme choice
         str = sche_Choi.getSelectedItem(); //                is selected,
         scheme="Up-strem";                      //   defort scheme,
         if(str == "LaxWendrof") scheme="LaxWendrof";
         else if(str=="CIP")     scheme="CIP";
      }
      initFlow();                            //  initialize flow field,
      Cour_t.setText(new Double(Cour).toString()); // update Courant No.
   }

   //  ** Execution method for the thread,
   //        practice run() when thread.start() is called.
   public void run(){
      double maxT=(xmax-xmin-0.2)/u0;  // maximum iteration number
      while( simulator != null){ // sml_start ){   // When flug is true, run
         if(t >= maxT) initFlow(); //  initialize every iteration
         drawFlow();               //  draw on the buffer memory
         //  Calculate for next step values
         solvAdvec();              //  compute by first upstream scheme
         if(scheme == "LaxWendrof")  LaxWendrof(); // Lax-Wendroff
         else if(scheme=="CIP")      solvCip();    // CIP method
         repaint();                // call up paint method for renewal
         t+=dt;                    // update time
         try{ simulator.sleep(delay);     // take a pause
            }catch(InterruptedException e){} // (it is necessary to sleep)
      }
      // simulator = null;
   }

   // **　Starting method for the thread start.
   /* public void start(){  //
      if(simulator == null){  sml_start = true;	 }
   } */

   // **  Stopping method for the thread killed.
   // public void stop(){   sml_start= false;       }

   //  **  Paint method --this method is called by update() --
   public void paint(Graphics g){
	  g.drawImage(buffer, 0, 0, this);  // copy buffer image
   }

   //  **  Update mathod -- this method is called by repaint() --
   public void update(Graphics g){
	  paint(g);
   }

   /***** Calculation of advection equation by finite difference method ***/
   //  **   setting up of initial values
   public void initFlow(){
      f  = new double [Nx+1];   // arrangement for calculated results
      f0 = new double [Nx+1];   // values of one step ahead
      fb = new double [Nx+1];
      fb0= new double [Nx+1];
      fc = new double [Nx+1];
      fc0= new double [Nx+1];
      g  = new double [Nx+1];
      g0 = new double [Nx+1];
      fa = new double [Nx+1];
      Cour=u0*dt/dx;              // update of Courant number
      t=0.0;                      // initialization of elapsed time
      for(int i=0; i<=Nx; i++){   // install initial values
         if(i < (int)(0.2/dx)) f[i] = 1.0; else f[i] = 0.0;
      }
      for(int i=0; i<=Nx; i++) f0[i]=f[i];
      for(int i=0; i<=Nx; i++) fb0[i]=f[i];
      for(int i=0; i<=Nx; i++) fc0[i]=f[i];
      for(int i=0; i<=Nx; i++) fa[i]=f[i];
      resetImage();            // delete the buffer image
      drawFlow();
      repaint();
   }
   // **  Upstream scheme
   public void solvAdvec(){
      for(int i=1; i<=Nx; i++)  f[i]=f0[i]-Cour*(f0[i]-f0[i-1]);
      setBoundary();   //  install baundary condition
      f0=nextStep( f );
   }
   // **  Lax-Wendroff scheme
   public void LaxWendrof(){
      for(int i=1; i<Nx; i++){
         fb[i]=fb0[i]-0.5*Cour*((fb0[i+1]-fb0[i-1])
               -Cour*(fb0[i+1]-2.0*fb0[i]+fb0[i-1]));
         fb[0]=1.0;   //  install baundary condition
      }
      fb0 = nextStep(fb);
   }
   // **  CIP method
   public void solvCip(){
      double xx,fdif,xam1,xbm1;
      for(int i=1; i<Nx; i++){
         xx=-u0*dt;
         fdif=(fc0[i]-fc0[i-1])/dx;
         xam1=(g0[i]+g0[i-1]-2.0*fdif)/(dx*dx);
         xbm1=(-3.0*fdif+2.0*g0[i]+g0[i-1])/dx;
         fc[i]=((xam1*xx+xbm1)*xx+g0[i])*xx+fc0[i];
         g[i]=(3.0*xam1*xx+2.0*xbm1)*xx+g0[i];
      }
      fc[0]=1.0;  g[0]=0.0; fc[Nx]=fc[Nx-1]; g[Nx]=g[Nx-1];
      fc0 = nextStep(fc);  g0 = nextStep(g);
   }

   // **  Install method of baundary value
   public void setBoundary(){
	  f[0]=1.0;
   }
   // **  Permutation method for next step
   public double[] nextStep(double fun[]){
      double fd []=new double [Nx+1];
      for(int i=0; i<=Nx; i++)  fd[i]=fun[i];
      return fd;
   }

   /*****  Methods for various graphics   */
   //　** 　Drawing method for the results obtained
   public void drawFlow(){
      int  i, x0, x1, y0, y1;
      resetImage();      //  delete the display
      axis();            //  draw axes
      for(i=1; i<=Nx; i++)
         //  draw the strict solution by black color
         if(i< (int)((0.2+u0*t)/dx+0.5)) fa[i]=1.0;
         else fa[i]=0.0;
         drawCurve(fa, Color.black);
         bg.drawLine(50,wy-30,100,wy-30);
         if(t==0.0)  bg.drawString(" initial value",105,wy-25);
		 else        bg.drawString(" Strict solution,",105,wy-25);
         bg.drawRect(0,0,wx-1,wy-1);

         //  draw the results by Lax-Wendroff method
         if(scheme == "LaxWendrof"){
            drawCurve(fb0,Color.blue);
            bg.drawLine(400,wy-30,470,wy-30);
            bg.drawString(" Lax-Wendroff scheme",475,wy-25);
         }else if(scheme=="CIP"){ // draw the results by CIP method
            drawCurve(fc0,Color.blue);
            bg.drawLine(400,wy-30,470,wy-30);
            bg.drawString(" CIP scheme",475,wy-25);
        }
        //  draw the results by the first upstream scheme
        drawCurve(f0,Color.red);
        bg.drawLine(200,wy-30,270,wy-30);
        bg.drawString(" 1st upstream scheme,",275,wy-25);
   }
   //  **  Drawing method of curves by indicated color
   public void drawCurve(double y[], Color col){
      int x0,x1,y0,y1;
      x0=xtr(x[0]);      y0=ytr(f[0]);
      bg.setColor(col);
      for(int i=1; i<Nx; i++){
         x1=xtr(x[i]);  y1=ytr(y[i]);
         bg.drawLine(x0,y0,x1,y1);
         x0=x1;         y0=y1;
      }
   }
   //  **　Setting up method of drawing range in the physical coordinate
   public void getRange(double x0, double y0, double x1, double y1){
      this.xmin = x0;       this.ymin = y0;
      this.xmax = x1;       this.ymax = y1;
   }
   //  **  Delete method of display by background color
   public void resetImage(){
      bg.setColor(Color.white);
      bg.fillRect(0, 0, wx, wy);
   }
   //  **  Method for calculating coordinate transformation coefficients
   public void trans(){
      double sx_size = 0.8, sy_size = 0.7;// drawing rates in x- & y-direction
      double sx0, sx1, sy0, sy1;
      sx0 = wx*0.5*(1.0-sx_size);   sx1 = wx-sx0;
      sy1 = wy*0.5*(1.0-sy_size);   sy0 = wy-sy1;
      sx_t0 =(sx0*xmax-sx1*xmin)/(xmax-xmin);
	  sx_t1 =(sx1-sx0)/(xmax-xmin);
      sy_t0 =(sy0*ymax-sy1*ymin)/(ymax-ymin);
	  sy_t1 =(sy1-sy0)/(ymax-ymin);
   }
   //  **  Transformation method of x-coordinate
   int xtr(double x){
      return (int)(sx_t1*x+sx_t0);
   }

   //  **  Transformation method of y-coordinate
   int ytr(double y){
      return (int)(sy_t1*y+sy_t0);
   }

   //  **  Drawing method for coordinate axes
   public void axis(){
      double v;
      int  i, charaW, charaH;
      int  xd = 5,  yd = 4;     // number of scales in  x- & y-axis
      int  m = 4;               // length of scale
      bg.setColor(Color.black);
      String st;
      Font font = new Font("TimesRoman", Font.BOLD, 10);
      charaH  = (getFontMetrics(font)).getLeading();
      charaH -= (getFontMetrics(font)).getDescent();
      charaH -= (getFontMetrics(font)).getAscent();

      bg.drawLine(xtr(xmin), ytr(0.0), xtr(xmax), ytr(0.0)); // x-axis
      bg.drawLine(xtr(0.0), ytr(ymin), xtr(0.0), ytr(ymax)); // y-axis
      bg.drawString("x", xtr(xmax)+10, ytr(0.0));            // name of x-axis
      bg.drawString("f", xtr(0.0), ytr(ymax)-10);            // name of y-axis

      for(i=0; i<=xd; i++){             // draw x-axis
         v = xmin+(double)i*(xmax-xmin)/(double)xd;
         st = String.valueOf(v);
         charaW = (getFontMetrics(font)).stringWidth(st);
	 bg.drawLine(xtr(v), ytr(0.0)-m, xtr(v), ytr(0.0)+m);
         bg.drawString(st,xtr(v)-charaW/2,ytr(0.0)-2*charaH);
      }
      for(i=0; i<=yd; i++){             // draw y-axis
         v = ymin+(double)i*(ymax - ymin)/(double)yd;
         st = String.valueOf(v);
         charaW = (getFontMetrics(font)).stringWidth(st);
         bg.drawLine(xtr(0.0)-m,ytr(v),xtr(0.0)+m,ytr(v));
         bg.drawString(st,xtr(0.0)-charaW,ytr(v)+charaH/2);
      }
   }
}
