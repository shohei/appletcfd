/****************************************************************/
/*   array.java :  Example of Java applet program               */
/*                                 to calculate a determinant   */
/*                        Last update: November 21,2001         */
/*                                     K. Minemura              */
/************************************************************** */
import java.awt.*;
import java.applet.*;

public class array extends Applet{
   int ii=3, jj=3;		
   double x[][] = new double [ii][jj];
   double result=0.0, a=0.0;
		
   public void init(){		
      for(int i=0; i<ii; i++){
         for(int j=0; j<jj; j++){
            x[i][j]=a;
            a++;
         }
      }
   }

   double array(double x[][], double result){
      result=x[0][0]*x[1][1]*x[2][2]+x[1][0]*x[2][1]*x[0][2]
             +x[2][0]*x[0][1]*x[1][2]-x[0][2]*x[1][1]*x[2][0]
             -x[1][0]*x[2][2]*x[0][1]-x[0][0]*x[1][2]*x[2][1];
      return result;
   }
	
   public void paint(Graphics g){
      for(int i=0; i<ii; i++){
         for(int j=0; j<jj; j++){
            g.drawString(String.valueOf(x[i][j]),50*(i+1),50*(j+1));
         }
      }
      g.drawString(String.valueOf("=   " + result), 200, 100);
   }
}
