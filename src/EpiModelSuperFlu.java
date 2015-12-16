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
public class EpiModelSuperFlu extends EpiModel {
    
    boolean alive = true;
    //Poisson poissDist; 
    Normal stndNorm;
    
    double beta; double rho;
    double N; double nu; 
    double fNoise; double stndNormRand; double fTerm; double eta;
    double tau;
    
    double dSI; double dIR;
    double newS; double newI; double newC;
    double initS; double initI; double kappa;
    
    @Override
    public void setInitParams()
    {
        //Set params to initial values and determine their order in params list
        beta = 2.0 * (1/3.8); params.add(beta); estParams.add(0);
        rho = 1.0; params.add(rho); estParams.add(1);
        nu = 1.0/3.8; params.add(nu); fixedParams.add(2);
        N = 2e07; params.add(N); fixedParams.add(3);
        kappa = 0.05; params.add(kappa); estParams.add(4); //reporting rate
        initS = 0.65; params.add(initS); estParams.add(5);
        initI = 0.0005; params.add(initI); estParams.add(6);
        
        tau = 0.20; params.add(tau); fixedParams.add(7);
                        
    }
    
    @Override
    public void setRandomGenerator() 
    {
      
        cern.jet.random.engine.RandomEngine engine = new cern.jet.random.engine.MersenneTwister(new java.util.Date());
        //poissDist = new Poisson(0.0, engine);
        stndNorm = new Normal(0.0, 1.0, engine);
        
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
        rho = params.getQuick(1);
        kappa = params.getQuick(4);
        initS = params.getQuick(5);
        initI = params.getQuick(6);
        //tau = params.getQuick(7);
        
        
    }
    
    @Override
    public DoubleArrayList getInitialConditions() 
    {
        
        DoubleArrayList currInits = new DoubleArrayList();
        //double initS = 130000.0; currInits.add(initS);
        currInits.add(N * initS);
        //double initI = 100.0; currInits.add(initI);
        currInits.add(N * initI);
        double initInc = 0.0; currInits.add(initInc); //for cumulative incidence
        
        //lastIncidence = 0.0;
        
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
        if ((params.get(5) + params.getQuick(6)) > 1.0) {
            constraintCheckFail = true;
        }
        
        return constraintCheckFail;
    }
    
    @Override
    public double getObservationProb(DoubleMatrix1D currStates, DoubleMatrix1D prevStates, double obsvVal)
    {
        
        double prob = 1.0;
        double newIncidence = currStates.getQuick(2) - prevStates.getQuick(2);
        double deviation = obsvVal - (newIncidence * kappa);
        double standDev = tau * newIncidence * kappa; //Math.sqrt(newIncidence * kappa);
        double exponent = -Math.pow(deviation,2.0) / (2.0 * Math.pow(standDev,2.0));
        prob *= Math.exp(exponent)/(standDev * Math.sqrt(2*Math.PI));
        
        return prob;
    }
    
//    @Override
//    public double getObservationProb(DoubleMatrix1D currStates, DoubleMatrix1D prevStates, double obsvVal)
//    {
//        //With binomial probability density
//        double prob = 1.0;
//        double n = currStates.getQuick(2) - prevStates.getQuick(2); // this is N
//        double k = obsvVal;
//        
//        double logBinomial = Gamma.logGamma(n+1) - Gamma.logGamma(n-k+1) - Gamma.logGamma(k+1);
//        double binomial = Math.exp(logBinomial);
//        if (Double.isInfinite(binomial)) {
//            
//            double logKappaToK = k * Math.log(kappa);
//            double logOneMinusKappa = (n-k) * Math.log(1-kappa);
//            double logProb = logBinomial + logKappaToK + logOneMinusKappa;
//            //double logProb = logBinomial + Math.log(Math.pow(kappa, k)) + Math.log(Math.pow((1-kappa), (n-k)));
//            prob *= Math.exp(logProb);
//        } else {
//            prob *= binomial * Math.pow(kappa, k) * Math.pow((1-kappa), (n-k));
//        }
//        
//        return prob;
//    }
    
    @Override
    //This simulates the current model!!
    public DoubleMatrix3D updateStates(int jParticles, DoubleMatrix3D matrix, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
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
                stndNormRand = stndNorm.nextDouble();
                fTerm = (1/beta) * Math.sqrt((rho*rho)/matrix.getQuick(1, j, xIndex)); 
                dSI = beta * (1.0 + (fTerm * stndNormRand)) * matrix.getQuick(0, j, xIndex) * matrix.getQuick(1, j, xIndex) * dt / N;
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
    
       
}
