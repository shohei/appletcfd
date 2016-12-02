/**=====================================================================
   Main frame of flow simulator  " BetaFlow " for heat convection 
           Benard.java   Preparation of control panel             
                         All rights reserved, Copyright (C) 2001-2003,
       Ver.2.0.     Last update December 25,2002, Kiyoshi Minemura.
=======================================================================*/
import java.applet.Applet;
import java.awt.*;
import java.awt.event.*;

// =====  main program of Benard  =====
public class Benard extends Applet{
   //         implements AdjustmentListener{
   // private ControlPanel panel;  
   // private DrawCanvasB  canvas; 
   private Scrollbar scrH, scrV;  
	
   //     Initializing method
   public void init(){
      //  classes used
      Grid gr = new Grid(); 
      Data da = new Data();
      Metric me = new Metric();   
      DrawCanvasB canvas = new DrawCanvasB( ); 
      SolverB so = new SolverB( canvas, me );
      ControlPanel panel = new ControlPanel( canvas, gr, da, me, so ); 

      // generation of GUI layout using BorderLayout
      setLayout(new BorderLayout()); 
      scrH = new Scrollbar(Scrollbar.HORIZONTAL, 0,10,-200,200);
      scrV = new Scrollbar(Scrollbar.VERTICAL,   0,10,-200,200);
      add( canvas, "Center" );    add( panel, "North" );
      add( scrH, "South" );       // scrH.addAdjustmentListener(this);  
      add( scrV, "East" );        // scrV.addAdjustmentListener(this);  
      Button ku = new Button(""); add( ku,"West" ); // left-side frame(no use)
   }

   /*//   Control method for responding to scrollbar operation
      public void adjustmentValueChanged(AdjustmentEvent ev){
        if(ev.getAdjustable() == scrH){
          canvas.trans_x = -scrH.getValue();  // get horizontal displacement
        }
        else if(ev.getAdjustable() == scrV){
          canvas.trans_y = -scrV.getValue();  // get vartical displacement
        }
        canvas.init_tr(); // recalculate coordinate transformation coefficients
    } */
}

//  ======== Control panel for the computation ===================
class ControlPanel extends Panel 
	           implements ActionListener, ItemListener {
   private DrawCanvasB canvas; 
   private Grid gr;        private Data da;  
   private SolverB so;     private Metric me;
   private Button btn_GRID, btn_INIT, btn_START, btn_STOP, btn_ReStart; 
   private Choice c_Re, c_Ra, c_Pr;   //   c_upst, c_graph; 
   private String flow = "duct";
   
   //  Constructor 
   public ControlPanel( DrawCanvasB canvas, Grid gr, Data da, 
                        Metric me, SolverB so ){
      this.canvas = canvas;    this.gr = gr;      this.da = da; 
      this.me = me;            this.so = so;

      initPanel();
   }

   //   Method for setting up GUI components
   public void initPanel(){    
      setBackground( SystemColor.activeCaptionBorder );
      Label la_Re, la_Ra, la_Pr, la_v; // la_upstream,la_graph,la_sample,scale;

      // install choice for selecting calculation object
      // c_obj=new Choice();  
      // c_obj.addItem("duct");c_obj.addItem("bend");c_obj.addItem("cylinder");
      // install choice for selecting Rayleigh number
      c_Ra=new Choice();   
      c_Ra.addItem(" 2500"); c_Ra.addItem(" 5000"); c_Ra.addItem("10000");
      c_Ra.addItem("20000"); c_Ra.addItem("40000");
      c_Ra.select(2);
      // install choice for selecting Reynolds number
      c_Re=new Choice();   
      c_Re.addItem("  50"); c_Re.addItem("  75");
      c_Re.addItem(" 100"); c_Re.addItem(" 250"); c_Re.addItem(" 500");
      c_Re.select(1);
      // install choice for selecting Rrandtle number
      c_Pr=new Choice();  
      c_Pr.addItem(" 0.5"); c_Pr.addItem(" 1.0"); c_Pr.addItem(" 2.5");
      c_Pr.select(1);
      // install choice for selecting graph
      /* c_graph = new Choice();
         c_graph.addItem("vectors");  c_graph.addItem("temper");
         c_graph.addItem("vorticity");  //  c_graph.addItem("pressure"); 
         c_graph.select(1); */
      /* set up of choice for selecting upstream scheme
         c_upst = new Choice();
         c_upst.addItem("QUICK"); c_upst.addItem("donor-cell"); 
         c_upst.addItem("first-order");
         c_upst.select(0);  */
      /* c_scale = new Choice();
         c_scale.addItem("original"); c_scale.addItem("enlarge");
      */

      // put of GUI components using LayoutManager ( GridBagLayout )
      GridBagLayout gb = new GridBagLayout(); 
      setLayout(gb); 
      GridBagConstraints gbc = new GridBagConstraints();
      //   
      gbc.ipadx=0;
      gbc.fill = GridBagConstraints.BOTH; // expand the panel to both directions
      gbc.anchor = GridBagConstraints.EAST;
      gbc.insets = new Insets(2,2,2,2);   // margin(picel) to up,down,left,right 
      /*
        gbc.gridy=0; gbc.gridx=0; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(la_sample = new Label("  SAMPLE=>"),gbc);
        la_sample.setForeground(Color.blue);    add(la_sample);
        gbc.gridy=0; gbc.gridx=1; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(c_obj, gbc);         add(c_obj);
      */
      gbc.gridy=0; gbc.gridx=1; gbc.gridwidth=1; gbc.gridheight=1;
      gb.setConstraints(btn_GRID =new Button(" GRID "),gbc);
      add(btn_GRID);

      gbc.gridy=0; gbc.gridx=2; gbc.gridwidth=1; gbc.gridheight=1;
      gb.setConstraints(btn_INIT = new Button("INITIAL V"), gbc);
      add(btn_INIT);

      gbc.gridy=0; gbc.gridx=3; gbc.gridwidth=1; gbc.gridheight=1;
      gb.setConstraints(btn_START = new Button("START"),gbc);
      add(btn_START);

      gbc.gridy=0; gbc.gridx=4; gbc.gridwidth=1; gbc.gridheight=1;
      gb.setConstraints(la_v = new Label("    "),gbc);
      la_v.setForeground(Color.blue);    add(la_v);

      gbc.gridy=0; gbc.gridx=5; gbc.gridwidth=1; gbc.gridheight=1;
      gb.setConstraints(la_Ra = new Label("Ra. No.=",Label.RIGHT),gbc);
      la_Ra.setForeground(Color.blue);    add(la_Ra);

      gbc.gridy=0; gbc.gridx=6; gbc.gridwidth=1; gbc.gridheight=1;
      gb.setConstraints(c_Ra, gbc);       add(c_Ra);	
		
      gbc.gridy=0; gbc.gridx=7; gbc.gridwidth=1; gbc.gridheight=1;
      gb.setConstraints(la_Pr = new Label("Pr. No.=",Label.RIGHT),gbc);
      la_Pr.setForeground(Color.blue);    add(la_Pr);
		
      gbc.gridy=0; gbc.gridx=8; gbc.gridwidth=1; gbc.gridheight=1;
      gb.setConstraints(c_Pr, gbc);       add(c_Pr);	

      gbc.gridy=0; gbc.gridx=9; gbc.gridwidth=1; gbc.gridheight=1;
      gb.setConstraints(la_Re=new Label("Re. No.=",Label.RIGHT),gbc);
      la_Re.setForeground(Color.blue);       add(la_Re);
		
      gbc.gridy=0; gbc.gridx=10; gbc.gridwidth=1; gbc.gridheight=1;
      gb.setConstraints(c_Re,gbc);           add(c_Re);

      /*
        gbc.gridy=0; gbc.gridx=5; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(la_scale = new Label("  scale="),gbc);
        la_scale.setForeground(Color.blue);    add(la_scale);

        gbc.gridy=0; gbc.gridx=6; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(c_scale, gbc);       add(c_scale);	
      */
      gbc.gridy=1; gbc.gridx=3; gbc.gridwidth=1; gbc.gridheight=1;
      gb.setConstraints(btn_STOP=new Button("STOP"),gbc);
      add(btn_STOP);

      gbc.gridy=1; gbc.gridx=2; gbc.gridwidth=1; gbc.gridheight=1;
      gb.setConstraints(btn_ReStart=new Button(" ReStart "),gbc);
      add(btn_ReStart);

      /* gbc.gridy=1; gbc.gridx=4; gbc.gridwidth=1;
        gb.setConstraints(la_graph = new Label("  graph="),gbc);
        la_graph.setForeground(Color.blue);    add(la_graph);
		
        gbc.gridy=1; gbc.gridx=5; gbc.gridwidth=1;
        gb.setConstraints(c_graph, gbc);       add(c_graph);

        gbc.gridy=1; gbc.gridx=6; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(la_upstream = new Label("Scheme"),gbc);
        la_upstream.setForeground(Color.blue);    add(la_upstream);

        gbc.gridy=1; gbc.gridx=7; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(c_upst,gbc);        add(c_upst);		
      */
      btn_GRID.addActionListener( this );
      btn_INIT.addActionListener( this );
      btn_START.addActionListener( this );
      btn_STOP.addActionListener( this );
      btn_ReStart.addActionListener( this );

      // c_obj.addItemListener(this); c_graph.addItemListener(this);
      // c_upst.addItemListener(this); c_scale.addItemListener(this);
      c_Re.addItemListener( this );
      c_Ra.addItemListener( this );
      c_Pr.addItemListener( this );

      // set up of response condition of button
      btn_GRID.setEnabled(true);     btn_INIT.setEnabled(false);
      btn_START.setEnabled(false);   btn_STOP.setEnabled(false);
      btn_ReStart.setEnabled(false); c_Re.setEnabled(false);
      c_Ra.setEnabled(false);        c_Pr.setEnabled(false);
      // c_obj.setEnabled(true);     c_graph.setEnabled(false);     
      // c_upst.setEnabled(false);   c_scale.setEnabled(false);
   }

   //   Control method for responding to choice events
   public void itemStateChanged( ItemEvent ev ){
      /*if(ev.getSource() == c_obj){
          flow = c_obj.getSelectedItem();
          canvas.repaint();
          btn_GRID.setEnabled( true );
      }*/
      if( ev.getSource() == c_Ra ){
         String st = c_Ra.getSelectedItem();
         da.Ra = Double.valueOf(st.trim()).doubleValue();
         btn_START.setEnabled( true );
      }
      if( ev.getSource() == c_Pr ){
         String st = c_Pr.getSelectedItem();
         da.Pr = Double.valueOf(st.trim()).doubleValue();
      }
      if( ev.getSource() == c_Re ){
         String st = c_Re.getSelectedItem();
         da.Re = Double.valueOf(st.trim()).doubleValue();
      }/*else if(ev.getSource() == c_graph){
         so.graph=c_graph.getSelectedItem(); 
        }else if(ev.getSource() == c_upst){	
         if(c_upst.getSelectedItem().equals("donor-cell"))
            so.scheme = "donor-cell";
         else if(c_upst.getSelectedItem().equals("first-order"))
            so.scheme = "first-order";
         else if(c_upst.getSelectedItem().equals("QUICK"))
            so.scheme = "QUICK";
      }*//*else if(ev.getSource() == c_scale){
            if(c_scale.getSelectedItem().equals("original"))
               scale=true;  
            else if(c_scale.getSelectedItem().equals("enlarge"))
               scale= false;
            canvas.callScale(scale);
      }*/
   }

   //   Control method for responding to button events
   public void actionPerformed( ActionEvent ev ){ 
      // When GRID button is pressed, 
      if( ev.getSource() == btn_GRID ){ 
         canvas.setStep( 1 );  
         so.stopSolver();       // stop the calculation
         gr.selectGrid( flow ); // select Grid for flow
         canvas.setRange( gr.getRange() );
         canvas.setGrid( gr.getGrid() );  
         canvas.setXpYpArray( gr.getXpArray(), gr.getYpArray() );
         // so.setNSmax( da.getNSmax() );
         // canvas.trans_x=0;   canvas.trans_y=0;

         btn_INIT.setEnabled( true );     btn_GRID.setEnabled( false );
         btn_START.setEnabled( false );   btn_STOP.setEnabled( false );
         btn_ReStart.setEnabled( false ); c_Re.setEnabled( false );
         c_Ra.setEnabled( false );        c_Pr.setEnabled( false );
         // c_graph.setEnabled(false);    c_upst.setEnabled(false); 
         // c_obj.setEnabled(true);   
      }
      //  When INIT button is pressed,
      if( ev.getSource() == btn_INIT ){
         canvas.setStep( 2 );          da.selectData( flow ); 
         da.setGrid( gr.getGrid() );   so.setGrid( gr.getGrid() );
         canvas.setVarray( da.getUParray(), da.getVParray() );

         btn_GRID.setEnabled(true);    btn_START.setEnabled(true); 
         btn_INIT.setEnabled(false);   c_Re.setEnabled(true);     
         c_Ra.setEnabled(true);        c_Pr.setEnabled(true);
         // c_graph.setEnabled(true);  c_upst.setEnabled(true);
      }
      // When START button is pressed, calculate the flow
      if( ev.getSource() == btn_START ){ 
         canvas.setStep( 3 );   
         me.setXYarray( gr.getXarray(), gr.getYarray() );
         me.setGrid( gr.getGrid() );      me.metric();  //  calculate metrics,
         canvas.setBCarray( da.getBCupp(), da.getBClow() );
         so.setBC( da.getBCupp(), da.getBClow(), da.getBCrig(), da.getBClef());
         so.setParameters( da.getParameters() ); so.setNSmax( da.getNSmax() );

         so.initSolver( flow );    //  initialize variables,
         so.startSolver();         //  start calculation.
         btn_STOP.setEnabled( true );  btn_ReStart.setEnabled( false );
         btn_INIT.setEnabled( false ); btn_GRID.setEnabled( true );
         c_Ra.setEnabled( false );     c_Pr.setEnabled( false );  
         c_Re.setEnabled( false );   
         // if(flow=="cylinder")  c_scale.setEnabled(true);
         // else c_scale.setEnabled(false);
      }
      // When STOP button is pressed,
      if( ev.getSource() == btn_STOP ){ 
         canvas.setStep( 4 );
         so.stopSolver();                   //        kill thread.
         btn_STOP.setEnabled( false );   btn_ReStart.setEnabled( true );
      }
      // When ReStart button is pressed,
      if( ev.getSource() == btn_ReStart ){ 
         so.startSolver();
         btn_STOP.setEnabled( true );   btn_ReStart.setEnabled( false );
      }
      canvas.repaint();
   } 
} 
