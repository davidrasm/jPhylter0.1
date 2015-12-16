import cern.colt.list.*;
import cern.colt.matrix.*;

/**
 * Coalescent model for influenza with super spreading
 * @author David
 */
public class CoalModelSuperFlu extends CoalModel {
    
    double beta; double S; double I; double N; double rho;
    
    public double computeRate(DoubleArrayList theta, DoubleMatrix1D states, double currTime) 
    {
        //This is the pairwise coalescent rate!
        
        beta = theta.getQuick(0);
        rho = theta.getQuick(1);
        N = theta.getQuick(3);
        
        S = states.getQuick(0);
        I = states.getQuick(1);
        
        //Pairwise coal rate
	double rate = 2.0 * beta * (S/N) * (1.0/I) * (1.0 + ((rho*rho)/(beta*beta))); //is beta^2 corrrect?
        
        return rate;
    }
    
}
