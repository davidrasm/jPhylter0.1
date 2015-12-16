import cern.colt.list.*;
import java.util.ArrayList;
import cern.jet.random.*;
import cern.colt.matrix.*;
import cern.jet.stat.Gamma;

/**
 * Specifies the epidemiological model
 * Model for influenza with super spreaders (i.e. transmission heterogeniety)
 * @author David
 */
public class EpiModelExpGrowth extends EpiModel {
    
    boolean alive = true;
    Poisson poissDist; 
    //Normal stndNorm;
    
    double beta; 
    double delta;
    double N;
    double initTime;
    
    //double fNoise; double stndNormRand; double fTerm; double eta;
    //double tau;
    
    double dSI; double dIR;
    double newI;
    double initI;
    
    @Override
    public void setInitParams()
    {
        //Set params to initial values and determine their order in params list
        beta = 0.55; params.add(beta); estParams.add(0);
        delta = 0.5; params.add(delta); estParams.add(1);
        initTime = 100.0; params.add(initTime); estParams.add(2); 
        initI = 1.0; params.add(initI); fixedParams.add(3);
        N = 2e07; params.add(N); fixedParams.add(4);
                        
    }
    
    @Override
    public void setRandomGenerator() 
    {
      
        cern.jet.random.engine.RandomEngine engine = new cern.jet.random.engine.MersenneTwister(new java.util.Date());
        poissDist = new Poisson(0.0, engine);
        //stndNorm = new Normal(0.0, 1.0, engine);
        
    }
    
    @Override
    public void updateEstParams(DoubleArrayList newParams)
    {
        int p = estParams.size();
        int paramListIndex;
        for (int i = 0; i < p; i++) {
            paramListIndex = estParams.get(i);
            params.setQuick(paramListIndex, newParams.get(i));
        }
        
        //Update the estimated params
        beta = params.getQuick(0);
        delta = params.getQuick(1);
        initTime = params.getQuick(2);
        
        
    }
    
    @Override
    public double getPriorProb()
    {
        //Compute log prior
        double mu = 1.0;
        double prior = 1.0; 
        if (rootTime > initTime) {
            double deviation = Math.log(rootTime - initTime) - mu; 
            double standDev = 1.25;
            double exponent = -Math.pow(deviation,2.0) / (2.0 * Math.pow(standDev,2.0));
            prior = Math.exp(exponent)/((rootTime - initTime) * standDev * Math.sqrt(2*Math.PI));
        } else {
            prior = 0.0;
        }
        
        
        double logPrior = Math.log(prior);
        return logPrior;
    }
    
    @Override
    public DoubleArrayList getInitialConditions() 
    {
        
        DoubleArrayList currInits = new DoubleArrayList();
        currInits.add(initI);

        return currInits;
        
    }
    
    @Override
    public boolean checkParamConstraints()
    {
        boolean constraintCheckFail = false;
        
        int p = estParams.size();
        int paramListIndex;
        //Check if any params are negative
        for (int i = 0; i < p; ++i) {
            paramListIndex = estParams.get(i);
            if (params.get(paramListIndex) <= 0.0) {
                constraintCheckFail = true;
            }
        }
            
        //Impose futher constraints
        //initS + initI has to be less than one
        if ((params.getQuick(0) > 2.0) || (params.getQuick(1) > 2.0)) {
            constraintCheckFail = true;
        }
        
        if (params.get(2) < 0.0 || params.get(2) > rootTime) {
            constraintCheckFail = true;
            //System.out.println("InitTime proposal out of bounds!");
        }
        
        return constraintCheckFail;
    }
    
    @Override
    public double getObservationProb(DoubleMatrix1D currStates, DoubleMatrix1D prevStates, double obsvVal)
    {
        
        double prob = 1.0;
        //double newIncidence = currStates.getQuick(2) - prevStates.getQuick(2);
        //double deviation = obsvVal - (newIncidence * kappa);
        //double standDev = tau * newIncidence * kappa; //Math.sqrt(newIncidence * kappa);
        //double exponent = -Math.pow(deviation,2.0) / (2.0 * Math.pow(standDev,2.0));
        //prob *= Math.exp(exponent)/(standDev * Math.sqrt(2*Math.PI));
        
        return prob;
    }
    
    @Override
    //This simulates the current model!!
    public DoubleMatrix3D updateStates(int jParticles, DoubleMatrix3D matrix, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
    {
        
        //Exponential growth with demographic stoch
        //Order of state variables is (0) I
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            double dt = xDtTimes.getQuick(n);
            int nextIndex = xIndex + 1;
            //double timeNow = xSubTimes.get(n);
            
            for (int j = 0; j < jParticles; ++j) 
            {
                
                if (xIndex >= startTimeIndex) { //for delayed start with set initial conditions
                
                    //INFECTIONS
                    //stndNormRand = stndNorm.nextDouble();
                    //fTerm = (1/beta) * Math.sqrt((rho*rho)/matrix.getQuick(1, j, xIndex));
                    dSI = poissDist.nextInt(beta * matrix.getQuick(0,j,xIndex) * dt);
                
                    //RECOVERIES
                    dIR = poissDist.nextInt(delta * matrix.getQuick(0,j,xIndex) * dt);

                    newI = matrix.getQuick(0,j,xIndex) + dSI - dIR;
                    //if (newI < 0) {
                        //newI = 0.0;
                        //System.out.println("Hit negative I");
                    //}
                    matrix.setQuick(0,j,nextIndex,newI);
                
                }

            }	
        }
        return matrix;
    }
    
    /**
    @Override
    public DoubleMatrix3D updateStatesDeterministic(int jParticles, DoubleMatrix3D matrix, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
    {
        
        //SIR with parameter noise
        //Order of state variables is (0) S, (1) I
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            double dt = xDtTimes.getQuick(n);
            int nextIndex = xIndex + 1;
            //double timeNow = xSubTimes.get(n);
            
            for (int j = 0; j < jParticles; ++j) 
            {
                //INFECTIONS
                //stndNormRand = stndNorm.nextDouble();
                //fTerm = (1/beta) * Math.sqrt((rho*rho)/matrix.getQuick(1, j, xIndex)); 
                dSI = beta * matrix.getQuick(0, j, xIndex) * matrix.getQuick(1, j, xIndex) * dt / N;
                if (dSI < 0.0) {
                    dSI = 0.0;
                }
                
                //RECOVERIES
                dIR = (nu * matrix.getQuick(1,j,xIndex) * dt);

                newS = matrix.getQuick(0,j,xIndex) - dSI;
                newI = matrix.getQuick(1,j,xIndex) + dSI - dIR;
                newC = matrix.getQuick(2,j,xIndex) + dSI;

                matrix.setQuick(0,j,nextIndex,newS);
                matrix.setQuick(1,j,nextIndex,newI);
                matrix.setQuick(2,j,nextIndex,newC);

            }	
        }
        return matrix;
    }
    */
       
}