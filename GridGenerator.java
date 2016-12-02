/**====================================================================
   Java Applet for generating Mesh                                
        GridGenerator.java  ( main program )                                       
             All Rights Reserved, Copyright (C) 2001-2002, K. MINEMURA  
       　　　　　　 Last updated by K. Minemura on November 5, 2002.  
======================================================================*/
import java.applet.Applet;
import java.awt.*;
import java.awt.event.*;

public class GridGenerator extends Applet
                           implements ItemListener, ActionListener{
	
   private Grid gr = new Grid();  // generate Grid object
   private DrawCanvas drw;
   private Graphics g;
   private Choice choice;
   private int Width, Height;
   private Button btnA, btnB; 
   private boolean step = false;
        
   public void init(){
      add(new Label("Select=> ")); 
      choice = new Choice();
      choice.addItem( "Duct" );     choice.addItem( "Elbow" );
      choice.addItem( "Cylinder" ); choice.addItem( "Cylinder_C" );
      choice.select( 0 );           
      add( choice );                choice.setEnabled( true );
      choice.addItemListener( this ); 
	  
      btnA = new Button( "Adjust-1" );  add( btnA );
      btnB = new Button( "Adjust-2" );  add( btnB );
      btnA.setEnabled( false );         btnB.setEnabled( false );
      btnA.addActionListener( this );
      btnB.addActionListener( this );
	  
      setBackground( Color.white ); 
      setForeground( Color.black );
      Width = getSize().width;      
      Height = getSize().height;
	  
      g = getGraphics();
      drw  = new DrawCanvas( g, Width, Height );
   }

   public void itemStateChanged(ItemEvent ev){
      String prob = "Duct"; 
      if(ev.getSource() == choice){ 
         prob = choice.getSelectedItem();  step=false;
      }
      if(prob =="Duct"){          
         gr.gridDuct();  
         btnA.setEnabled( false );  btnB.setEnabled( false );
      }
      if(prob == "Elbow"){   
         gr.gridElbow();  
         btnA.setEnabled( false );  btnB.setEnabled( false );
      }
      if(prob == "Cylinder"){  
         gr.gridO();    
         btnA.setEnabled( false );  btnB.setEnabled( false );
      }
      if(prob == "Cylinder_C"){    
         gr.gridCylinder(); 
         btnA.setEnabled( true );   btnB.setEnabled( false );
      }
      drw.winPort( gr.getRange() );
      repaint();
   }

   public void actionPerformed( ActionEvent ev ){	   
      if(ev.getSource() == btnA){
         btnA.setEnabled( false );  btnB.setEnabled( true );
         gr.gridE();    step = true;
      }
      if(ev.getSource() == btnB){
         btnB.setEnabled( false );
         gr.gridEP();   step = true;
      }
      repaint();
   }

   //  === paint method ===  
   public void paint(Graphics g){ 
      gr.linkNode(); 
      drw.drawGrid( gr.getMx(), gr.getMy(), gr.getLink(), gr.getXarray(), gr.getYarray() );
      if(step == true){
         g.setColor(Color.black);
         g.drawString("  k=" + gr.getIter(), 10, 420); 
         g.drawString("  residual = " + gr.getResid(), 10, 435);
      }
   }
}