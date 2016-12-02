/*   Preliminary program for numerical simulation No. 2
　　===　Particles rebound on walls  ====
                     Last update: October 10,2001 by K. Minemura  */
import java.applet.*;
import java.awt.*;
import java.awt.event.*;

public class Anime2 extends Applet   // name this applet "anime2.java"
	   implements ActionListener, MouseListener, Runnable{
   int  r = 5;               //  particle radius
   int  k=0, ico;            //  k=particle number, ico=color number
   int width, height;        // size of applet window
   Image buff;               // off-screen image for double buffering
   Graphics bg;              // Graphics object for off-screen
   Thread th = null;         // dead thread for animation
   Button  bt_clear;         // start button for animation
   int px[]=new int [25]; int py[]=new int[25]; // particle coordinates
   int pdx[]=new int[25]; int pdy[]=new int[25];// unit length of particle movement
   Color pco[]=new Color[25]; // color of each particle
   Color co[]={Color.red, Color.blue, Color.green}; // variable for colors

   // Initializing method, called when applet starts
   public void init(){
      width = getSize().width;          // get the size of applet width
      height = getSize().height;
      setBackground(Color.lightGray);   // color the background light gray

      buff = this.createImage(width, height); // create the buffer memory
      bg = buff.getGraphics();        // create the off-screen Graphics object
      bt_clear = new Button("CLEAR"); add(bt_clear);// create and put button

      bt_clear.addActionListener(this); // register the event listener
      addMouseListener(this);           // register the mouse listener
   }

   // Managing method when the button is pressed
   public void actionPerformed(ActionEvent ev){
      if(ev.getSource() == bt_clear){// When clear button is pressed,
         k=0;  ico=0;          //  re-initialize numbers of particle and color
         th = null;            //  nullify the thread object
          //bg.clearRect(0,0,width,height);
          repaint();
      }
   }

   // Managing method for the mouse which is pressed
   public void mousePressed(MouseEvent ev){
      if(th == null){             // When the thread object is not alive,
         th = new Thread(this);   //   create thread object
         th.start();              //   initiate start action of the thread
         repaint();               //   call update method through repaint method
      }
      if(k<25){                   // When released particle numbers are less than 25,
         px[k]=ev.getX();  py[k]=ev.getY(); // get the mouse position
         pdx[k]=8;   pdy[k]=5;  //  set initial unit length for movement
         pco[k]=co[ico];        //  designate particle color
         k++;
      }else {th = null;} // nullify the thread when k>25.
   }
   // Following three methods are necessary though they are not used.
   public void mouseReleased(MouseEvent ev){}
   public void mouseClicked(MouseEvent ev){}
   public void mouseEntered(MouseEvent ev){}
   public void mouseExited(MouseEvent ev){
      ico++; if(ico>2) ico=0;  // change color within the range of three.
   }

   // Paint method for displaying the off-screen image on the applet
   public void paint(Graphics g){
      bg.clearRect(0,0,width,height);  // delete painted image
      for(int i=0; i<k; i++){          // repeat k-times
         bg.setColor(pco[i]);          //   change particle color
         bg.fillOval(px[i]-r,py[i]-r, r*2, r*2); // draw the particle
      }
      g.drawImage(buff, 0, 0, this);   // copy buffer image on applet
   }

   // Update method for painting without deleting
   public void update( Graphics g ){
      paint(g);
   }

   //　Run method for animation
   public void run(){
      while( th != null){ // execute when the thread object lives
	 for(int i=0; i<k; i++){  // decide unit length of particle movement
            if((px[i]-r+pdx[i]<0) || (px[i]+r+pdx[i]> width))  pdx[i]=-pdx[i];
            if((py[i]-r+pdy[i]<0) || (py[i]+r+pdy[i]> height)) pdy[i]=-pdy[i];
               px[i]+=pdx[i];	py[i]+=pdy[i]; // update particle coordinates
        }
        repaint();
        try{ Thread.sleep(100);  }        // take a rest for 0.01 seconds
        catch( InterruptedException e){ } // sense exceptions
        }
   }
}