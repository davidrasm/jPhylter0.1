import cern.colt.matrix.*;
import cern.colt.list.*;
//import java.util.ArrayList;

/**
 * For two-population model used in GENETICS paper with coupling scalar M
 * @author David
 */
public class StructCoalModel2Pop extends StructCoalModel {
    
    
    @Override
    public void make() 
    {   
        states = 2; 
        
        DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
        F = factory2D.make(states, states);
        F.assign(0);
        G = factory2D.make(states, states);
        G.assign(0);
        
        DoubleFactory1D factory1D;
        factory1D = DoubleFactory1D.dense;
        Y = factory1D.make(states);
        Y.assign(0);
    }
    
    @Override
    public void updateF(DoubleMatrix1D xCurr, DoubleArrayList theta, double currTime) 
    {
        
        //Full spatial SIR model with seasonality in focal and global pops
        
        //double totalN = xCurr.get(3);
        double beta = theta.get(0);
        double alphaFocal = theta.get(1);
        double alphaGlobal = theta.get(2);
        double deltaFocal = theta.get(3) * 365.25;
        double deltaGlobal = theta.get(4) * 365.25;
        double M = theta.get(5);
        
        //For sinuisoidal forcing
        double seasNow = ((currTime + deltaFocal)/365.25) - Math.floor(currTime/365.25);
        double seasHeight = Math.cos(2*Math.PI * seasNow);
        //if (seasHeight < 0) {
            //seasHeight = 0;
        //}
        double betaNowFocal = beta * (1 + alphaFocal*seasHeight);
        
        seasNow = ((currTime + deltaGlobal)/365.25) - Math.floor(currTime/365.25);
        seasHeight = Math.cos(2*Math.PI * seasNow);
        //if (seasHeight < 0) {
            //seasHeight = 0;
        //}
        double betaNowGlobal = beta * (1 + alphaGlobal*seasHeight);
        
        double fFF = betaNowFocal * (1-M) * xCurr.get(0) * xCurr.get(2) / xCurr.get(4); //rate of transmission within focal pop
        double fFG = betaNowGlobal * M * xCurr.get(1) * xCurr.get(2) / xCurr.get(5);
        double fGF = betaNowFocal * M * xCurr.get(0) * xCurr.get(3) / xCurr.get(4); //rate of transmission from J -> A
        double fGG = betaNowGlobal * (1-M) * xCurr.get(1) * xCurr.get(3) / xCurr.get(5);
        
        F.setQuick(0,0,fFF);
        F.setQuick(0,1,fFG);
        F.setQuick(1,0,fGF);
        F.setQuick(1,1,fGG);
    }
    
    @Override
    public void updateG(DoubleMatrix1D xCurr, DoubleArrayList theta)
    {
        //G Matrix is always zero
    }
    
    @Override
    public void updateY(DoubleMatrix1D xCurr) 
    {
        Y.setQuick(0,xCurr.get(2)); //infections in focal pop
        Y.setQuick(1,xCurr.get(3)); //infected in global pop    
    }
    
}
