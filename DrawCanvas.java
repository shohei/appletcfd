/**===================================================================
    Simulation of a Convection Equation
    　 DrawCanvas.java  ( Graphics Parts )
             All Rights Reserved, Copyright (C) 2001-2002, K. MINEMURA
       　　　　　　 Last updated by K. Minemura on November 4, 2002.
=====================================================================*/
import java.awt.*;

//    Graphics Class for Convection Equation
public class DrawCanvas extends Canvas{
   //  Declaration of instance variables related GUI and graphics
   private Graphics  bg;    //  graphics object for buffer image
   private Image  buffer;   //  object of buffer image
   private int    wx, wy;   //  width and height of applet(pixel unit)
   private double xmin,xmax,ymin,ymax; // drawing range of physical space
   private double sx_t0,sx_t1,sy_t0,sy_t1; // transformation coefficients

   // Constructor
   public DrawCanvas(Image im, Graphics g, int width, int height){
      buffer = im;  bg = g;  wx = width; wy = height;
   }

   // === transformation method for coordinates between physical and screen spaces ===
   public void winPort( double range[]){
      xmin=range[0]; xmax=range[1]; ymin=range[2]; ymax=range[3];
      trans();
   }

   // === method for drawing a curve with designated color ===
   public void drawCurve(double x[], double y[], int N, Color col){
      int xp[] = new int [ N ],  yp[] = new int [ N ];
      for(int i=0; i<N; i++)  xp[i]=xtr( x[i] );
      for(int i=0; i<N; i++)  yp[i]=ytr( y[i] );
      bg.setColor( col );
      bg.drawPolyline( xp, yp, N );
   }

   // === delete method of buffered image using background color ===
   public void resetImage(){
      bg.clearRect( 0, 0, wx, wy );
   }

   //  === access method for buffer memory ===
   public Image getBuffer(){
      return buffer;
   }

   // === calculating method for coordinate transformation coefficients ===
   public void trans(){
      double sx_size = 0.75, sy_size = 0.65; // drawing rates in x- & y-direction
      double sx0, sx1, sy0, sy1;
      sx0 = wx*0.5*(1.0-sx_size);   sx1 = wx-sx0;
      sy1 = wy*0.5*(1.0-sy_size);   sy0 = wy-sy1;
      sx_t0 =(sx0*xmax-sx1*xmin)/(xmax-xmin);
      sx_t1 =(sx1-sx0)/(xmax-xmin);
      sy_t0 =(sy0*ymax-sy1*ymin)/(ymax-ymin);
      sy_t1 =(sy1-sy0)/(ymax-ymin);
   }

   // === transformation method for x-coordinate ===
   public int xtr( double x ){
      return (int)(sx_t1*x+sx_t0);
   }

   // === transformation method for y-coordinate ===
   public int ytr( double y ){
      return (int)(sy_t1*y+sy_t0);
   }

   // === drawing method for coordinate axes ===
   public void axis( String strX, int divX, String strY, int divY ){
      double v;
      int  i, charaW, charaH;
      int  m=4;              // length of scale
      bg.setColor(Color.black);
      String st;
      Font font = new Font("TimesRoman", Font.PLAIN, 12);
      bg.setFont( font );
      FontMetrics fm=bg.getFontMetrics( font );
      charaH = fm.getHeight();

      bg.drawLine(xtr(xmin), ytr(0.0), xtr(xmax), ytr(0.0)); // x-axis
      bg.drawLine(xtr(0.0), ytr(ymin), xtr(0.0), ytr(ymax)); // y-axis
      bg.drawString( strX, xtr(xmax)+10, ytr(0.0));          // name of x-axis
      bg.drawString( strY , xtr(0.0), ytr(ymax)-10);         // name of y-axis

      for(i=0; i<=divX; i++){             // draw x-axis
         v = xmin+(double)i*(xmax-xmin)/(double)divX;
         st = String.valueOf(v);
         charaW = fm.stringWidth( st );
         bg.drawLine(xtr(v), ytr(0.0)-m, xtr(v), ytr(0.0)+m);
         bg.drawString( st, xtr(v)-2*charaW/3, ytr(0.0)+charaH+2 );
      }
      for(i=0; i<=divY; i++){              // draw y-axis
         v = ymin+(double)i*(ymax - ymin)/(double)divY;
         st = String.valueOf(v);
         charaW = fm.stringWidth( st );
         bg.drawLine(xtr(0.0)-m, ytr(v), xtr(0.0)+m,ytr(v) );
         bg.drawString(st, xtr(0.0)-3*charaW/2, ytr(v)+charaH/2 );
      }
   }
}
