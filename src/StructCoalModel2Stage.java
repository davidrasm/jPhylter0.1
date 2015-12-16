import cern.colt.matrix.*;
import cern.colt.list.*;

/**
 * Subclass for 2 stage acute/chronic model
 * @author David
 */
public class StructCoalModel2Stage extends StructCoalModel {
    
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
    public void updateF(DoubleMatrix1D xCurr, DoubleArrayList theta, double timeCurr) 
    {
        
        double totalN = theta.get(4); //change this if pop size can change
        double fA = theta.get(0) * xCurr.get(0) * xCurr.get(1) / totalN; //rate of transmission from J -> J
        double fC = theta.get(1) * xCurr.get(0) * xCurr.get(2) / totalN; //rate of transmission from J -> A
        
        F.setQuick(0,0,fA);
        F.setQuick(0,1,0);
        F.setQuick(1,0,fC);
        F.setQuick(1,1,0);
    }
    
    @Override
    public void updateG(DoubleMatrix1D xCurr, DoubleArrayList theta)
    {
        double gAC = theta.get(2) * xCurr.get(1);
        
        //G.assign(0.0);
        G.setQuick(0,1,gAC);
    }
    
    @Override
    public void updateY(DoubleMatrix1D xCurr) 
    {
        Y.setQuick(0,xCurr.get(1)); //acutes
        Y.setQuick(1,xCurr.get(2)); //chronics    
    }
    
}
