import cern.colt.matrix.*;
import cern.colt.list.*;

/**
 * Subclass for 3 stage early/chronic/AIDS model
 * @author David
 */
public class StructCoalModelSECA extends StructCoalModel {
    
    @Override
    public void make() 
    {   
        states = 3; 
        
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
        currTime = timeCurr;
        double tauScaler = 1.0;
        //if (timeCurr >= 729756) { //if time is >= 1998
            //tauScaler = 1 - theta.get(9);
        //}
        
        double N = xCurr.get(4);
        double transScaler = Math.exp(-theta.get(3) * (xCurr.getQuick(1) + xCurr.getQuick(2) + xCurr.getQuick(3))); //alpha is param 3
        double fEE = tauScaler * transScaler * theta.get(1) * xCurr.getQuick(0) * xCurr.getQuick(1) / N; //betaChronic is param 0
        double fCE = tauScaler * transScaler * theta.get(0) * xCurr.getQuick(0) * xCurr.getQuick(2) / N;
        double fAE = tauScaler * transScaler * theta.get(2) * xCurr.getQuick(0) * xCurr.getQuick(3) / N;
        
        F.setQuick(0,0,fEE);
        F.setQuick(1,0,fCE);
        F.setQuick(2,0,fAE);
        //all other elements are zero in F
    }
    
    @Override
    public void updateG(DoubleMatrix1D xCurr, DoubleArrayList theta)
    {
        double tauScaler = 1.0;
        //if (currTime >= 729756) { //if time is >= 1998
            //tauScaler = 1 - theta.get(9);
        //}
        double gEC = tauScaler * theta.get(4) * xCurr.getQuick(1);
        double gCA = tauScaler * theta.get(5) * xCurr.getQuick(2);
        
        G.assign(0.0);
        G.setQuick(0,1,gEC);
        G.setQuick(1,2,gCA);
        
    }
    
    @Override
    public void updateY(DoubleMatrix1D xCurr) 
    {
        Y.setQuick(0,xCurr.getQuick(1)); //early infections
        Y.setQuick(1,xCurr.getQuick(2)); //chronic
        Y.setQuick(2,xCurr.getQuick(3)); //AIDS
    }
    
}
