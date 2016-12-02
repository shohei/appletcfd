/***********************************************************************/
/*   Main frame of flow simulator  " BetaFlow "                        */
/*           BetaFlow class:  Preparation of control panel             */
/*                            Last update November 9,2001.             */
/*           All rights reserved, Copyright (C) 2001, Kiyoshi Minemura */ 
/***********************************************************************/
import java.applet.Applet;
import java.awt.*;
import java.awt.event.*;

//  ====== Main program of BetaFlow ==========
public class BetaFlow extends Applet
            implements AdjustmentListener{
   private ControlPanel panel;  
   private DrawCanvas   canvas; 
   private Scrollbar scrH, scrV;  
	
   //  Initializing method 
   public void init(){
      Grid gr = new Grid();     Data da = new Data();
      Metric me = new Metric(); Solver so = new Solver();
      // setting up of GUI layout using BorderLayout()
      setLayout(new BorderLayout());  //  
      scrH = new Scrollbar(Scrollbar.HORIZONTAL, 0,10,-200,200);
      scrV = new Scrollbar(Scrollbar.VERTICAL,   0,10,-200,200);
      canvas = new DrawCanvas(scrH, scrV, gr, da, me, so); 
      panel  = new ControlPanel(canvas, gr, da, me, so);  
      add(canvas,"Center"); add(panel,"North");
      scrH.addAdjustmentListener(this);  add(scrH, "South");
      scrV.addAdjustmentListener(this);  add(scrV, "East");
      Button ku = new Button("");    add(ku,"West"); // left-side frame(no use)
   }

   //   Control method for responding to scrollbar operation
   public void adjustmentValueChanged(AdjustmentEvent ev){
     if(ev.getAdjustable() == scrH){
        canvas.trans_x = -scrH.getValue();  // get horizontal displacement
     }
     else if(ev.getAdjustable() == scrV){
        canvas.trans_y = -scrV.getValue();  // get vartical displacement
     }
     canvas.init_tr();    // recalculate coordinate transformation coefficients
   } 
}

//  ======== Instalation class of control panael ===================
class ControlPanel extends Panel   
	   implements ActionListener,ItemListener{
   private DrawCanvas canvas; private Metric me;
   private Data da;  private Solver so;  private Grid gr; 
					   
   private Button btn_GRID, btn_INIT, btn_FLOW, btn_STOP, btn_ReStart; 
   private Choice c_obj, c_graph, c_Re, c_upst, c_scale; 
   private String flow;
   
   //  constructor 
   public ControlPanel(DrawCanvas canvas, Grid gr, Data da, 
                       Metric me, Solver so){
      this.canvas=canvas;      this.gr=gr;      this.da=da; 
      this.me= me;  this.so=so; 
      init();
   }

   //   Initializing method for setting up GUI components
   public void init(){     
      setBackground(SystemColor.activeCaptionBorder);
      Label la_sample, la_Re, la_void, la_graph, la_upstream, la_scale; 
   
      c_obj=new Choice();  // install choice for selecting calculation object
         c_obj.addItem("duct"); c_obj.addItem("bend"); c_obj.addItem("cylinder");
      c_Re=new Choice();   // install choice for selecting Reynolds number
         c_Re.addItem("  100"); c_Re.addItem("  300");
         c_Re.addItem("  500"); c_Re.addItem("  800"); c_Re.addItem(" 1500");
         c_Re.select(1);
      c_graph = new Choice(); // install choice for selecting graph
         c_graph.addItem("vectors");  c_graph.addItem("re_vectors"); 
         c_graph.addItem("pressure"); c_graph.addItem("vorticity");
      c_upst = new Choice();  // set up of choice for selecting upstream scheme
         c_upst.addItem("QUICK"); c_upst.addItem("donor-cell"); 
         c_upst.addItem("first-order");
      c_scale = new Choice(); // install choice for selecting display scale
         c_scale.addItem("original"); c_scale.addItem("enlarge");
   
      // put of GUI components using LayoutManager ( GridBagLayout)
      GridBagLayout gb = new GridBagLayout(); 
      setLayout(gb); 
      GridBagConstraints gbc = new GridBagConstraints();
      //   
      gbc.ipadx=0;
      gbc.fill = GridBagConstraints.BOTH;  // expand this panel to both directions
      gbc.anchor = GridBagConstraints.EAST;
      gbc.insets = new Insets(2,2,2,2);  // margin(pixel) to up,down,left,right

      gbc.gridy=0; gbc.gridx=0; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(la_sample = new Label("  SAMPLE=>"),gbc);
        la_sample.setForeground(Color.blue);    add(la_sample);
		
      gbc.gridy=0; gbc.gridx=1; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(c_obj, gbc);         add(c_obj);

      gbc.gridy=0; gbc.gridx=2; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(btn_GRID =new Button("GRID"),gbc);
	add(btn_GRID);

      gbc.gridy=0; gbc.gridx=4; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(btn_INIT = new Button("INITIAL V"), gbc);
        add(btn_INIT);
		
      gbc.gridy=0; gbc.gridx=5; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(la_scale = new Label("  scale="),gbc);
        la_scale.setForeground(Color.blue);    add(la_scale);
		
      gbc.gridy=0; gbc.gridx=6; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(c_scale, gbc);       add(c_scale);	

      gbc.gridy=1; gbc.gridx=0; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(btn_FLOW = new Button("FLOW"),gbc);
        add(btn_FLOW);

      gbc.gridy=1; gbc.gridx=1; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(btn_STOP=new Button("STOP"),gbc);
        add(btn_STOP);

      gbc.gridy=1; gbc.gridx=2; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(btn_ReStart=new Button(" ReStart "),gbc);
        add(btn_ReStart);

      gbc.gridy=1; gbc.gridx=3; gbc.gridwidth=1;
        gb.setConstraints(la_graph = new Label("  graph="),gbc);
        la_graph.setForeground(Color.blue);    add(la_graph);
		
      gbc.gridy=1; gbc.gridx=4; gbc.gridwidth=1;
        gb.setConstraints(c_graph, gbc);       add(c_graph);

      gbc.gridy=1; gbc.gridx=5; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(la_Re=new Label("     RE="),gbc);
        la_Re.setForeground(Color.blue);       add(la_Re);
	
      gbc.gridy=1; gbc.gridx=6; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(c_Re,gbc);           add(c_Re);

      gbc.gridy=1; gbc.gridx=7; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(la_upstream = new Label("Scheme"),gbc);
        la_upstream.setForeground(Color.blue);    add(la_upstream);

      gbc.gridy=1; gbc.gridx=8; gbc.gridwidth=1; gbc.gridheight=1;
        gb.setConstraints(c_upst,gbc);        add(c_upst);		

      btn_GRID.addActionListener(this);
      btn_INIT.addActionListener(this);
      btn_FLOW.addActionListener(this);
      btn_STOP.addActionListener(this);
      btn_ReStart.addActionListener(this);

      c_obj.addItemListener(this);
      c_graph.addItemListener(this);
      c_upst.addItemListener(this);
      c_Re.addItemListener(this);
      c_scale.addItemListener(this);

      // set up of response condition of button
      c_obj.setEnabled(true);        btn_GRID.setEnabled(false);     
      btn_INIT.setEnabled(false);  
      btn_FLOW.setEnabled(false);    btn_STOP.setEnabled(false);    
      btn_ReStart.setEnabled(false);
      c_graph.setEnabled(false);     c_Re.setEnabled(false);
      c_upst.setEnabled(false);      c_scale.setEnabled(false);
   }

   //   Control method for responding to choice
   public void itemStateChanged(ItemEvent ev){
      boolean scale=true;
      if(ev.getSource() == c_obj){
         flow = c_obj.getSelectedItem();
         canvas.repaint();
         btn_GRID.setEnabled(true);
      }else if(ev.getSource() == c_Re){
         String st = c_Re.getSelectedItem();
         da.Re = Double.valueOf(st.trim()).doubleValue();
      }else if(ev.getSource() == c_graph){
         so.graph=c_graph.getSelectedItem(); 
      }else if(ev.getSource() == c_upst){	
         if(c_upst.getSelectedItem().equals("donor-cell"))
            so.scheme = "donor-cell";
         else if(c_upst.getSelectedItem().equals("first-order"))
            so.scheme = "first-order";
         else if(c_upst.getSelectedItem().equals("QUICK"))
            so.scheme = "QUICK";
      }else if(ev.getSource() == c_scale){
         if(c_scale.getSelectedItem().equals("original"))
            scale=true;  
         else if(c_scale.getSelectedItem().equals("enlarge"))
            scale= false;
         canvas.callScale(scale);
      }
   }

   //   Control method for responding to button events
   public void actionPerformed(ActionEvent ev){
      if(ev.getSource() == btn_GRID){      // When GRID button is pressed, 
         canvas.step=1;  canvas.repaint(); //    generate grid.
         btn_INIT.setEnabled(true);    c_obj.setEnabled(true);	   
         btn_FLOW.setEnabled(false);      //     set button no response
         btn_STOP.setEnabled(false);   btn_ReStart.setEnabled(false);
         c_graph.setEnabled(false);    c_Re.setEnabled(false);
         c_upst.setEnabled(false);     btn_GRID.setEnabled(false);
		 
         so.solver_stopwater();        //  stop FLOW-execution
         gr.grid_select(flow);         //  select other grid for flow
         canvas.trans_x=0;   canvas.trans_y=0;
      }
      if(ev.getSource()==btn_INIT){    //  When INIT button is pressed,
         canvas.step=2;  
         btn_GRID.setEnabled(true);    btn_FLOW.setEnabled(true); 
         btn_INIT.setEnabled(false);   c_graph.setEnabled(true);     
         c_Re.setEnabled(true);        c_upst.setEnabled(true);
         da.Data_select(gr, flow);    //    install initila condition
      }

      //  ****** Flow calculation *****************************************
      if(ev.getSource()==btn_FLOW)  { // When FLOW button is pressed, 
         canvas.step=3;
         btn_STOP.setEnabled(true);  btn_ReStart.setEnabled(false);
         btn_INIT.setEnabled(false); btn_GRID.setEnabled(true);
         c_Re.setEnabled(true); 
         if(flow=="cylinder")  c_scale.setEnabled(true);
         else c_scale.setEnabled(false);

         me.metric(gr);                      //  calculate metrics,
         so.solver_init(canvas, gr, me, da); //  initialize variables,
         so.solver_start();                  //  start calculation.
      }
      if(ev.getSource() == btn_STOP){ // When STOP button is pressed, 
         btn_STOP.setEnabled(false);   btn_ReStart.setEnabled(true);  
         canvas.step=4;

         so.solver_stop();            //        kill thread.
      }	  
      if(ev.getSource() == btn_ReStart) { // When ReStart button is pressed,
         btn_STOP.setEnabled(true);   btn_ReStart.setEnabled(false);
		  
         so.solver_start();               //      restart calculation.
      }
      canvas.repaint();               // repaint
   }  // end of ActionPerformed
}  // end of UpsidePanel
