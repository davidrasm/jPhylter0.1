import cern.colt.list.*;
import cern.colt.matrix.*;

/**
 * Superclass for the coalescent model
 * Modify CoalModel according to process model
 * @author David
 */
public class CoalModel {
    
    public double computeRate(DoubleArrayList theta, DoubleMatrix1D states, double currTime) 
    {
        //This is Eric's lambda = 2f(t) / (I^2)
	double rate = 0.0; //compute coal rate
        return rate;
    }
    
}
