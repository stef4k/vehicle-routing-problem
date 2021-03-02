package vrp;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;


public class SolutionDrawer {
    
    static void drawRoutes(ArrayList <Node> allNodes, Solution s, String fileName) 
    {
        
        int VRP_Y = 800;
        int VRP_INFO = 200;
        int X_START = 200;
        int X_GAP = 600;
        int Y_START = 200;
        int Y_GAP = 800;
        int margin = 30, VEH_INDEX;
        int marginNode = 1;
        int drawMulti = 12;
        
        int XXX =  VRP_INFO + X_GAP;
        int YYY =  VRP_Y;
        
              
        BufferedImage output = new BufferedImage(XXX, YYY, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = output.createGraphics();
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, XXX, YYY);
        g.setColor(Color.BLACK);

        
        double minX = Double.MAX_VALUE;
        double maxX = Double.MIN_VALUE;
        double minY = Double.MAX_VALUE;
        double maxY = Double.MIN_VALUE;
        for (int i = 0; i < allNodes.size(); i++)
        {
            Node n = allNodes.get(i);
            if (n.x > maxX) maxX = n.x;
            if (n.x < minX) minX = n.x;
            if (n.y > maxY) maxY = n.y;
            if (n.y < minY) minY = n.y;
        }

        int mX = XXX - 2 * margin;
        int mY = VRP_Y - 2 * margin;

        int A, B;
        if ((maxX - minX) > (maxY - minY))
        {
            A = mX;
            B = (int)((double)(A) * (maxY - minY) / (maxX - minX));
            if (B > mY)
            {
                B = mY;
                A = (int)((double)(B) * (maxX - minX) / (maxY - minY));
            }
        }
        else
        {
            B = mY;
            A = (int)((double)(B) * (maxX - minX) / (maxY - minY));
            if (A > mX)
            {
                A = mX;
                B = (int)((double)(A) * (maxY - minY) / (maxX - minX));
            }
        }
        
        for (int r = 0; r < s.routes.size(); r++) {
            Route rt = s.routes.get(r);
            for (int i = 1; i < rt.nodes.size(); i++) {
                Node n;
                n = rt.nodes.get(i - 1);
                int ii1 = (int) ((double) (A) * ((n.x - minX) / (maxX - minX) - 0.5) + (double) mX / 2) + margin;
                int jj1 = (int) ((double) (B) * (0.5 - (n.y - minY) / (maxY - minY)) + (double) mY / 2) + margin;
                n = rt.nodes.get(i);
                int ii2 = (int) ((double) (A) * ((n.x - minX) / (maxX - minX) - 0.5) + (double) mX / 2) + margin;
                int jj2 = (int) ((double) (B) * (0.5 - (n.y - minY) / (maxY - minY)) + (double) mY / 2) + margin;

                g.drawLine(ii1, jj1, ii2, jj2);
            }
        }
        
        for (int i = 0; i < allNodes.size(); i++)
        {
            Node n = allNodes.get(i);

            int ii = (int)((double)(A) * ((n.x - minX) / (maxX - minX) - 0.5) + (double)mX / 2) + margin;
            int jj = (int)((double)(B) * (0.5 - (n.y - minY) / (maxY - minY)) + (double)mY / 2) + margin;
            if (i != 0)
            {
                g.fillOval(ii - 2 * marginNode, jj - 2 * marginNode, 4 * marginNode, 4 * marginNode);
                String id = Integer.toString(n.ID);
                g.drawString(id, ii + 8 * marginNode, jj+ 8 * marginNode);
            }
            else
            {
                g.fillRect(ii - 4 * marginNode, jj - 4 * marginNode, 8 * marginNode, 8 * marginNode);
                String id = Integer.toString(n.ID);
                g.drawString(id, ii + 8 * marginNode, jj + 8 * marginNode);
            }
        }
        
        String cst = "Cost: " + s.distance;
        g.drawString(cst, 10, 10);
        
        fileName = fileName + ".png";
        File f = new File(fileName);
        try 
        {
            ImageIO.write(output, "PNG", f);
        } catch (IOException ex) {
            Logger.getLogger(VRP.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
    
}