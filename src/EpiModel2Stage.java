import cern.colt.list.*;
import java.util.ArrayList;
import cern.jet.random.*;
import cern.colt.matrix.*;

/**
 * Specifies the epidemiological model
 * @author David
 */
public class EpiModel2Stage extends EpiModel {
    
    //Acute/chronic 2-stage SIR model used in original GENETICS manuscript

    Poisson poissDist; Normal stndNorm;
    double betaAcute; double betaChronic; double gamma; double nu; double N; //double fNoise;
    double dSIa; double dIaIc; double dIcR; double R;
    double newS; double newIa; double newIc; double newR; double newP; double pTotal;
    double transSIa; double transSIc;
    
    @Override
    public void setInitParams()
    {
        //Set params to initial values and determine their order in params list
        double beta = 0.00052;
        betaAcute = 100.0 * beta; params.add(betaAcute); estParams.add(0); //0.0304;
        betaChronic = beta; params.add(betaChronic); estParams.add(1); //0.0061;
        gamma = (1/(365.25/4.0)); params.add(gamma); fixedParams.add(2);
        nu = (1/(365.25/1.0)); params.add(nu); fixedParams.add(3);
        N = 2.0e05; params.add(N); fixedParams.add(4);
        //fNoise = 0.0012; params.add(fNoise); fixedParams.add(5);
                        
    }
    
    @Override
    public void setRandomGenerator() 
    {
      
        cern.jet.random.engine.RandomEngine engine = new cern.jet.random.engine.MersenneTwister(new java.util.Date());
        poissDist = new Poisson(0.0, engine);
        stndNorm = new Normal(0.0, 1.0, engine);
        
    }
    
    @Override
    public void updateEstParams(DoubleArrayList newParams)
    {
        int p = estParams.size();
        int paramListIndex;
        for (int i = 0; i < p; i++) {
            paramListIndex = estParams.get(i);
            params.set(paramListIndex, newParams.get(i));
        }
        
        //Update estimate params
        betaAcute = params.getQuick(0);
        betaChronic = params.getQuick(1);
        //gamma = params.getQuick(2);
        
    }
    
    @Override
    public DoubleArrayList getInitialConditions() 
    {
        
        //Order init conditions 1) S 2) Ia 3) Ic 4) N
        DoubleArrayList currInits = new DoubleArrayList();

        //At t = 100.0 days
        double initS = N - 1; currInits.add(initS);
        double initIa = 1.0; currInits.add(initIa);
        double initIc = 0.0; currInits.add(initIc);
        double initR = N - initS - initIa - initIc; currInits.add(initR);
        double initN = N; currInits.add(initN);
        
        currInits.add(0.0); //transSIa
        currInits.add(0.0); //transSIc
        currInits.add(0.0); //dIaIc
        currInits.add(0.0); //dIcR
        
//        x.setQuick(0,j,nextIndex,newS);
//        x.setQuick(1,j,nextIndex,newIa);
//        x.setQuick(2,j,nextIndex,newIc);
//        x.setQuick(3,j,nextIndex,newR);
//        x.setQuick(4,j,nextIndex,newP);
//        x.setQuick(5,j,nextIndex,transSIa);
//        x.setQuick(6,j,nextIndex,transSIc);
//        x.setQuick(7,j,nextIndex,dIaIc);
//        x.setQuick(8,j,nextIndex,dIcR);

        
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
        
        return constraintCheckFail;
    }
    
    @Override
    public DoubleMatrix3D updateStates(int jParticles, DoubleMatrix3D x, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
    {
        
        // DONT USE THIS METHOD BEFORE UPDATING
        
        //Two-stage SIR model with acute/chronic stages
        //Order of state variables is 1) S 2) Ia 3) Ic 4) N
      
        double eta; double dSIa; double dIaIc; double dIcR; double R;
        double newS; double newIa; double newIc; double newR; double newP; double pTotal;
        double fNoiseSIa; double fNoiseSIc;
        double transSIa; double transSIc;
		
        double betaAcute = params.getQuick(0);
	double betaChronic = params.getQuick(1);
        double gamma = params.getQuick(2);
        double nu = params.getQuick(3);
        double fNoise = params.getQuick(5);

	double stndNormRand;
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            double dt = xDtTimes.getQuick(n);
            
            double timeNow = xSubTimes.get(n);
            
            for (int j = 0; j < jParticles; ++j) {
                
                stndNormRand = stndNorm.nextDouble();
                eta = stndNormRand / Math.sqrt(Math.abs(dt)); //used for this particle and time only

                pTotal = x.getQuick(4,j,xIndex);

                //INFECTIONS FROM ACUTE CASES (INDEX 1)
                fNoiseSIa = fNoise * betaAcute * x.getQuick(0,j,xIndex) * x.getQuick(1,j,xIndex) * eta * dt / pTotal;
                transSIa = betaAcute * x.getQuick(0,j,xIndex) * x.getQuick(1,j,xIndex) * dt / pTotal;

                //INFECTIONS CHRONIC CASES
                fNoiseSIc = fNoise * betaChronic * x.getQuick(0,j,xIndex) * x.getQuick(2,j,xIndex) * eta * dt / pTotal;
                transSIc = betaChronic * x.getQuick(0,j,xIndex) * x.getQuick(2,j,xIndex) * dt / pTotal;
                
//                double checkSIa = fNoiseSIa + transSIa; //can this be positive?
//                if (checkSIa > 0.0) {
//                    System.out.println("Check SI flow");
//                }
//                double checkSIc = fNoiseSIc + transSIc; //can this be positive?
//                if (checkSIc > 0.0) {
//                    System.out.println("Check SI flow");
//                }

                dSIa = fNoiseSIa + transSIa + fNoiseSIc + transSIc;

                //TRANSITIONS FROM Ia TO Ic
                dIaIc = gamma * x.getQuick(1,j,xIndex) * dt;
                if (x.getQuick(2,j,xIndex) == 0.0) {
                    dIaIc = 0.0;
                }

                //RECOVERIES
                dIcR = nu * x.getQuick(2,j,xIndex) * dt;
                if (x.getQuick(3,j,xIndex) == 0.0) {
                    dIcR = 0.0;
                }

                newS = x.getQuick(0,j,xIndex) - dSIa;
                newIa = x.getQuick(1,j,xIndex) + dSIa - dIaIc;
                if (newIa < 0.0) {
                    newS -= Math.abs(newIa);
                    newIa = 0.0;
                }
                newIc = x.getQuick(2,j,xIndex) + dIaIc - dIcR;
                if (newIc < 0.0) {
                    newIa -= Math.abs(newIc);
                    newIc = 0.0;
                }
                newR = x.getQuick(3,j,xIndex) + dIcR;
                if (newR < 0.0) {
                    newIc -= Math.abs(newR); //subtract deficit in R from Ic
                    newR = 0.0;
                }
                
                newP = newS + newIa + newIc + newR;
                //newP = pTotal;
                
                //if (newP > 2.0e05) {
                    //System.out.println();
                //}

                int nextIndex = xIndex - 1;
                x.setQuick(0,j,nextIndex,newS);
                x.setQuick(1,j,nextIndex,newIa);
                x.setQuick(2,j,nextIndex,newIc);
                x.setQuick(3,j,nextIndex,newR);
                x.setQuick(4,j,nextIndex,newP);
                //System.out.println();

            }	
        }
        return x;
    }
    
    @Override
    public DoubleMatrix3D updateStatesDeterministic(int jParticles, DoubleMatrix3D x, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
    {
        
        //Two-stage SIR model with acute/chronic stages
        //Order of state variables is 1) S 2) Ia 3) Ic 4) N
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            int nextIndex = xIndex + 1;
            double dt = xDtTimes.getQuick(n);
            
            //double timeNow = xSubTimes.get(n);
            
            for (int j = 0; j < jParticles; ++j) {

                pTotal = x.getQuick(4,j,xIndex);

                //INFECTIONS FROM ACUTE CASES (INDEX 1)
                transSIa = betaAcute * x.getQuick(0,j,xIndex) * x.getQuick(1,j,xIndex) * dt / pTotal;

                //INFECTIONS CHRONIC CASES
                transSIc = betaChronic * x.getQuick(0,j,xIndex) * x.getQuick(2,j,xIndex) * dt / pTotal;

                dSIa = transSIa + transSIc;

                //TRANSITIONS FROM Ia TO Ic
                dIaIc = gamma * x.getQuick(1,j,xIndex) * dt;

                //RECOVERIES
                dIcR = nu * x.getQuick(2,j,xIndex) * dt;

                newS = x.getQuick(0,j,xIndex) - dSIa;
                newIa = x.getQuick(1,j,xIndex) + dSIa - dIaIc;
                newIc = x.getQuick(2,j,xIndex) + dIaIc - dIcR;
                newR = x.getQuick(3,j,xIndex) + dIcR;
                newP = newS + newIa + newIc + newR;

                x.setQuick(0,j,nextIndex,newS);
                x.setQuick(1,j,nextIndex,newIa);
                x.setQuick(2,j,nextIndex,newIc);
                x.setQuick(3,j,nextIndex,newR);
                x.setQuick(4,j,nextIndex,newP);
                x.setQuick(5,j,nextIndex,transSIa);
                x.setQuick(6,j,nextIndex,transSIc);
                x.setQuick(7,j,nextIndex,dIaIc);
                x.setQuick(8,j,nextIndex,dIcR);

            }	
        }
        return x;
    }
    
    @Override
    public DoubleMatrix3D updateStatesTauLeap(int jParticles, DoubleMatrix3D x, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
    {
        
        //Two-stage SIR model with acute/chronic stages
        //Order of state variables is 1) S 2) Ia 3) Ic 4) N
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            int nextIndex = xIndex + 1;
            double dt = xDtTimes.getQuick(n);
            
            //double timeNow = xSubTimes.get(n);
            
            for (int j = 0; j < jParticles; ++j) {

                pTotal = x.getQuick(4,j,xIndex);

                //INFECTIONS FROM ACUTE CASES (INDEX 1)
                transSIa = poissDist.nextInt(betaAcute * x.getQuick(0,j,xIndex) * x.getQuick(1,j,xIndex) * dt / pTotal);

                //INFECTIONS CHRONIC CASES
                transSIc = poissDist.nextInt(betaChronic * x.getQuick(0,j,xIndex) * x.getQuick(2,j,xIndex) * dt / pTotal);

                dSIa = transSIa + transSIc;

                //TRANSITIONS FROM Ia TO Ic
                dIaIc = poissDist.nextInt(gamma * x.getQuick(1,j,xIndex) * dt);

                //RECOVERIES
                dIcR = poissDist.nextInt(nu * x.getQuick(2,j,xIndex) * dt);

                newS = x.getQuick(0,j,xIndex) - dSIa;
                newIa = x.getQuick(1,j,xIndex) + dSIa - dIaIc;
                if (newIa < 0.0) {
                    newIa = 0.0;
                }
                newIc = x.getQuick(2,j,xIndex) + dIaIc - dIcR;
                if (newIc < 0.0) {
                    newIc = 0.0;
                }
                newR = x.getQuick(3,j,xIndex) + dIcR;
                newP = newS + newIa + newIc + newR;

                x.setQuick(0,j,nextIndex,newS);
                x.setQuick(1,j,nextIndex,newIa);
                x.setQuick(2,j,nextIndex,newIc);
                x.setQuick(3,j,nextIndex,newR);
                x.setQuick(4,j,nextIndex,newP);
                x.setQuick(5,j,nextIndex,transSIa);
                x.setQuick(6,j,nextIndex,transSIc);
                x.setQuick(7,j,nextIndex,dIaIc);
                x.setQuick(8,j,nextIndex,dIcR);

            }	
        }
        return x;
    }
    
    public DoubleArrayList fastSim(DoubleArrayList inits, double dt, double endTime)
    {
        DoubleArrayList newInits = new DoubleArrayList();
        return newInits;
    }
    
    @Override
    public double computeTransProb(DoubleMatrix1D xTimeNow, DoubleMatrix1D xTimePlusOne, double dt, double timeNow) 
    {
        double transProb = 1.0;
        
        transSIa = betaAcute * xTimeNow.getQuick(0) * xTimeNow.getQuick(1) * dt / pTotal;
        transSIc = betaChronic * xTimeNow.getQuick(0) * xTimeNow.getQuick(2) * dt / pTotal;
        dIaIc = gamma * xTimeNow.getQuick(1) * dt;
        dIcR = nu * xTimeNow.getQuick(2) * dt;
        
        //Implicitly sets transProb to 1.0 if mean for poisson density function is 0.0, which is technically undefined
        if (transSIa > 0.0) {
            poissDist.setMean(transSIa);
            transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(5));
        }
        if (transSIc > 0.0) {
            poissDist.setMean(transSIc);
            transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(6));
        }
        if (dIaIc > 0.0) {
            poissDist.setMean(dIaIc);
            transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(7));
        }
        if (dIcR > 0.0) {
            poissDist.setMean(dIcR);
            transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(8));
        }
        
        if (transProb <= 0.0) {
            transProb = Double.MIN_VALUE;
        }
        if (Double.isNaN(transProb)) {
            transProb = Double.MIN_VALUE;
            System.out.println("Pop state trans prob is NaN");
        }
        
        return transProb;
    }
    
       
}