import cern.colt.list.*;
import java.util.ArrayList;
import cern.jet.random.*;
import cern.colt.matrix.*;
//import cern.jet.random.AbstractDistribution;
//import cern.colt.matrix.DoubleMatrix1D;
//import cern.colt.matrix.DoubleMatrix3D;

/**
 * Specifies the epidemiological model
 * @author David
 */
public class EpiModel {
    
    //Superclass instance variables
    DoubleArrayList params = new DoubleArrayList();
    ArrayList<Integer> estParams = new ArrayList<Integer>(); //indexes of the estimated params
    ArrayList<Integer> fixedParams = new ArrayList<Integer>(); //indexes of the fixed params
    DoubleFactory1D factory1D = DoubleFactory1D.dense;
    DoubleFactory2D factory2D = DoubleFactory2D.dense;
    boolean alive;
    int startTimeIndex;
    double rootTime;
    
    
    public void setInitParams()
    {
        //Set params to initial values and determine their order in params list                
    }
    
    public void setRandomGenerator() 
    {
        
    }
    
    public DoubleArrayList getEstParams() 
    {
        DoubleArrayList currEstParams = new DoubleArrayList();
        int p = estParams.size();
        for (int i = 0; i < p; i++) {
            int paramListIndex = estParams.get(i);
            currEstParams.add(params.getQuick(paramListIndex));
        }
        return currEstParams;
        
    }
    
    public DoubleArrayList getFixedParams() 
    {
        DoubleArrayList currFixedParams = new DoubleArrayList();
        int p = fixedParams.size();
        int paramListIndex;
        for (int i = 0; i < p; i++) {
            paramListIndex = fixedParams.get(i);
            currFixedParams.add(params.getQuick(paramListIndex));
        }
        return currFixedParams;
        
    }
    
    public void updateEstParams(DoubleArrayList newParams)
    {
        int p = estParams.size();
        int paramListIndex;
        for (int i = 0; i < p; i++) {
            paramListIndex = estParams.get(i);
            params.setQuick(paramListIndex, newParams.get(i));
        }
        
    }
    
    public double getPriorProb()
    {
        double prior = 0.0;
        return prior;
    }
    
    public DoubleArrayList getInitialConditions() 
    {
        
        DoubleArrayList currInits = new DoubleArrayList();
        return currInits;
        
    }
    
    public boolean checkParamConstraints()
    {
        boolean constraintCheckFail = false;
        return constraintCheckFail;
    }
    
    public double getObservationProb(DoubleMatrix1D currStates, DoubleMatrix1D prevStates, double obsvVal)
    {
        double prob = 1.0;
        return prob;
    }
    
    public DoubleMatrix3D updateStatesTauLeap(int jParticles, DoubleMatrix3D matrix, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
    {
        return matrix;
    }
    
    public DoubleMatrix3D updateStatesDeterministic(int jParticles, DoubleMatrix3D matrix, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
    {
        return matrix;
    }
    
    public DoubleMatrix3D updateStates(int jParticles, DoubleMatrix3D matrix, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
    {
        return matrix;
    }
    
    public double computeTransProb(DoubleMatrix1D xTimeNow, DoubleMatrix1D xTimePlusOne, double dt, double timeNow) 
    {   
        double transProb = 1.0;
        return transProb;
    }
    
       
}
