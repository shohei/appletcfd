/**===================================================================
    Simulation of a Convection Equation                             
    　 ConvecApl.java  ( main program )                                             
             All Rights Reserved, Copyright (C)       K. MINEMURA   
       　　　　　　 Last updated by K. Minemura on Feb. 28, 2007.  
=====================================================================*/
import java.applet.Applet;
import java.awt.*;
import java.awt.event.*; 

//     Main Applet for Convection Equation
public class ConvecApl extends Applet implements 
                       ItemListener, ActionListener, Runnable{
	
   private ConvSolver sol = new ConvSolver();
   private DrawCanvas drw;
   private Thread thread;   
   
   //  Declaration of instance variables related GUI
   private Button bt_start, bt_stop, bt_reset; // button objects
   private Choice dt_Choi, sche_Choi;   // choice objects
   private TextField Cour_t, Diff_t;    // text-field objects
   
   private String scheme = "FTCS";
   private double x[], ys[];
   private Image buffer;
   private Graphics bg;
   private int Width, Height;
   
   //  ===　Initializing method of applet ===
   public void init(){
      add(new Label("dt=", Label.RIGHT));     //  time step
      dt_Choi = new Choice();    
      dt_Choi.addItem("0.01");    dt_Choi.addItem("0.02"); 
      dt_Choi.addItem("0.025");   dt_Choi.addItem("0.0255");
      dt_Choi.addItem("0.026");   dt_Choi.addItem("0.027");
      dt_Choi.addItem("0.028");
      dt_Choi.select(1);    sol.setDt( 0.02 ); // initial condition     
      add(dt_Choi);              
	  
      add(new Label("diff.=", Label.RIGHT));  // diffusion number
      Diff_t = new TextField(); 
      String str = "" + sol.getDiff();
      Diff_t.setText( str.substring(0,6) );   Diff_t.setEditable( false );
      add(Diff_t);
	  
      add(new Label("  Courant="));            // Courant Number
      Cour_t = new TextField();   
      str = ""+ sol.getCourant(); 
      Cour_t.setText( str.substring(0,5) );   Cour_t.setEditable( false );
      add( Cour_t );                
	  
      add(new Label("Scheme=", Label.RIGHT ));  //  scheme
      sche_Choi = new Choice();
      sche_Choi.addItem( "FTCS" );  sche_Choi.addItem( "P-C" );
      sche_Choi.addItem( "QUICK" ); sche_Choi.addItem( "K-K" );
      sche_Choi.select( 0 );   add( sche_Choi );
	  
      setBackground( Color.white );   
      // setForeground (Color.black );    
	  
      bt_start = new Button( "Start" );  add( bt_start ); 	  
      bt_stop  = new Button( "Stop" );   add( bt_stop );  	  
      bt_reset = new Button( "Reset" );  add( bt_reset ); 	  
	  
      dt_Choi.addItemListener( this );    
      bt_start.addActionListener( this );
      bt_stop.addActionListener( this );
      bt_reset.addActionListener( this );
      sche_Choi.addItemListener( this );
	  
      Width = getSize().width;       
      Height = getSize().height;      
      buffer = createImage( Width, Height );  
      bg  = buffer.getGraphics(); 
	  
      drw = new DrawCanvas( buffer, bg, Width, Height );
      drw.winPort( sol.getRange() ); 
	  
      sol.initFlow();  // set up initial values for calculation
      drawIntDistr();  // draw initial distribution
   }

   // ===  Managing method for the button action === 
   public void actionPerformed( ActionEvent ev ){
      if(ev.getSource() == bt_start){          // For start-button;
         if(thread == null) thread = new Thread( this );
         thread.start();
      }
      if(ev.getSource() == bt_stop){           // For stop-button;
         stop();
      }
      if(ev.getSource() == bt_reset){          // For reset-button; 
         stop();
         sol.initFlow();    // re-install the initial values  
         drawIntDistr();    // draw the initial condition
      }
      repaint();
   }

   // ==  Managing method for the choice selection ===
   public void itemStateChanged( ItemEvent ev ){
      String str;
      if(ev.getSource() == dt_Choi){             // For dt-choice;
         str = dt_Choi.getSelectedItem();            
         sol.setDt( Double.valueOf(str).doubleValue()); 
         Diff_t.setText( "" + sol.getDiff() );
         Cour_t.setText( "" + sol.getCourant() );
      }
      if(ev.getSource() == sche_Choi){           // For scheme-choice; 
         scheme = sche_Choi.getSelectedItem(); //     
      }
      stop();  
      sol.initFlow();                //  initialize flow field
   }
   
   //  === Execution method of the thread;
   //           practice run() when simulator.start() is called. ===
   public void run(){   
      int  delay = 35;             // sleep period( msec )
      if(sol.dt >= 0.02) delay = 50; 
      while( thread != null){      // When the thread lives, run on

         // calculate next step using the following schemes
         sol.solvFTCS();                   // FTCS scheme 
         sol.solvCip();                    // CIP scheme
         drawCurve( x, sol.getCipArray() );    

         if(scheme == "P-C"){       // Predictor-Corector scheme 
            sol.P_C();    
            drw.drawCurve( x, sol.getPCarray(), x.length, Color.green );
         } 
         if(scheme == "QUICK"){     // QUICK scheme
            sol.QUICK(); 
            drw.drawCurve( x, sol.getQuickArray(), x.length, Color.green );
         }
         if(scheme == "K-K"){       // Kuwahara-Kuwahara scheme
            sol.K_K();   
            drw.drawCurve( x, sol.getKKarray(), x.length, Color.green );
         }
         repaint(); 
         try { thread.sleep( delay ); }  // take a pause before next step 
         catch(InterruptedException e){} // 
      }
   }
   //  === stop method for kill the thread === 
   public void stop(){   thread = null; }
   //  === method for drawing the initial distribution ===
   public void drawIntDistr(){
      drw.resetImage();             //  delete the buffer memory
      x = sol.getXarray();          //  intake x-array
      ys = sol.getYsArray();        //  intake the initial distribution
      drw.axis( "x", 5, "z", 4 );   //  draw axes
      drw.drawCurve( x, ys, x.length, Color.darkGray ); //  initial distribution
      notes();                      //  draw notes
   }
   // === method for drawing the elementary results  ===
   public void drawCurve( double x[], double y[] ){
      drw.resetImage();
      drw.axis( "x", 5, "z", 4 );
      drw.drawCurve(  x, ys, x.length, Color.darkGray ); 
      drw.drawCurve( x, sol.getFtcsArray(), x.length, Color.red ); // FTCS
      drw.drawCurve( x, sol.getCipArray(), x.length, Color.blue ); // CIP
      notes(); addNotes();
  }
   // === method for adding notes === 
   public void notes(){
      bg.setColor( Color.darkGray );
      bg.drawLine(30, Height-30, 90, Height-30);
      bg.drawString(" initial value", 90, Height-25);
      bg.setColor( Color.red );
      bg.drawLine( 165, Height-30, 210, Height-30);
      bg.drawString(" FTCS scheme", 215, Height-25);
      bg.setColor( Color.blue );
      bg.drawLine(315, Height-30, 375, Height-30);
      bg.drawString(" CIP method", 380, Height-25);
   }
   // === method for adding additional notes ===
   public void addNotes(){
      bg.setColor( Color.green ); 
      bg.drawLine( 460, Height-30, 505, Height-30);	
      if(scheme == "P-C")
         bg.drawString(" Predictor-Corrector scheme", 510, Height-25);
      if(scheme == "QUICK")
         bg.drawString(" QUICK scheme", 510, Height-25);
      if(scheme == "K-K")
         bg.drawString(" Kawahara-Kuwahara scheme", 510, Height-25);
   }
   
   
   //  =====  paint method -- this method is call out by update() =====
   public void paint( Graphics g ){
      g.drawImage( drw.getBuffer(), 0, 0, this ); // copy the buffered image
   }

   //  ===  update method -- this method is call out by repaint() === 
   public void update( Graphics g ){   paint(g);   }
}