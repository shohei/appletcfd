/*   Preliminary program for numerical simulation  No. 1
　　  ===　Particle rebound on walls ====
                       Last update: October 10, 2001 bt K. Minemura
   <applet code="Anime.class" width="350" height=""50"></applet> */
import java.applet.Applet;
import java.awt.*;
import java.awt.event.*;   // management for events

public class Anime extends Applet  // name this applet "Anime.java"
       implements ActionListener, Runnable{  // interfaces for events
   int x = 250, y = 140, r = 10; // particle coordinate and its radius
   int dx = 8, dy = 5;           // unit length of particle movement
   int width, height;            // size of applet window
   Image buff;               // off-screen image for double buffering
   Graphics bg;              // Graphics object for off-screen
   Button bt_stop, bt_start; // buttons to start and stop for animation
   Thread th = null;         // initial condition of thread object is dead

   // Initializing method, read when applet starts
   public void init(){
      width  = getSize().width;            // get the size of applet width
      height = getSize().height;
      setBackground(Color.lightGray);      // color the background light gray
      buff = this.createImage(width, height); // create the buffer memory
      bg = buff.getGraphics();      // create the off-screen Graphics object
      bt_start= new Button("START");       // create START button
      add(bt_start);                       // put the button on applet
      bt_stop = new Button("STOP");  add(bt_stop);
      bt_start.addActionListener(this);    // register the event-listener
      bt_stop.addActionListener(this);
   }

   // Managing method when the buttons are pressed
   public void actionPerformed(ActionEvent ev){
      if(ev.getSource() == bt_start){  // When START button is pressed,
          if(th == null) th = new Thread(this);       //   create the thread object
          th.start();                  //   start run method for animation
      }else if(ev.getSource() == bt_stop){// When STOP button is pressed,
          th = null;                   //   nullify the thread object
      }
   }

   // Paint method for displaying the off-screen image on the applet
   public void paint(Graphics g){
      bg.clearRect(0, 0, width, height); // delete the applet display
      bg.setColor(Color.red);            // designate red color for painting
      bg.fillOval(x-r, y-r, r*2, r*2);   // draw particle with radius r on (x,y)
      g.drawImage(buff, 0, 0, this);     // copy buffer image on applet window
   }

   // Update method for dealing with paint method without delete of applet,
   //    which is called by repaint method.
   // It is necessary to over-write update method for double buffering.
   public void update( Graphics g ){
      paint(g);
  }

   //　Run method for animation
   public void run(){
      while( th != null ){      // execute when the thread lives.
         if((x-r+dx<0) || (x+r+dx> width))  dx=-dx;//　moving unit length
         if((y-r+dy<0) || (y+r+dy> height)) dy=-dy;
         x += dx;  y += dy;         // update particle coordinate
         repaint();                 // call update method for refreshing
         try{  Thread.sleep(100);  }           // take a rest for 0.01 seconds
         catch( InterruptedException e){ } // sense exceptions
      }
   }
}