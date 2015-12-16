import cern.colt.list.*;
import cern.colt.matrix.*;

/**
 * Coalescent model for influenza with super spreading
 * @author David
 */
public class CoalModelExpGrowth extends CoalModel {
    
    double beta; double I;
    
    public double computeRate(DoubleArrayList theta, DoubleMatrix1D states, double currTime) 
    {
        //This is the pairwise coalescent rate!
        
        beta = theta.getQuick(0);
        //delta = theta.getQuick(1);
        I = states.getQuick(0);
        
        //Pairwise coal rate
        double rate = 2.0 * (beta) / (I);
	//double rate = 2.0 * beta * (S/N) * (1.0/I) * (1.0 + ((rho*rho)/(beta*beta))); //is beta^2 corrrect?
        
        return rate;
    }
    
}
