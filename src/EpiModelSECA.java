import cern.colt.list.*;
import java.util.ArrayList;
import cern.jet.random.*;
import cern.colt.matrix.*;

/**
 * Specifies the epidemiological model
 * @author David
 */
public class EpiModelSECA extends EpiModel {
    
    //Early/Chronic/AIDS model for HIV
    Poisson poissDist; Normal stndNorm;
    DoubleArrayList currInits;
    
    double betaEarly; double betaChronic; double betaAIDS;
    double gammaEarly; double gammaChronic; double gammaAIDS; 
    double alpha; double mu; double N; double tau; double initE;
    double initTime;
    double cumulativeI;
    
    double transScaler; double transEE; double transCE; double transAE; 
    double dSE; double dEC; double dCA;
    double dBS; double dSD; double dED; double dCD; double dAD;
    double newS; double newE; double newC; double newA; double newP;
    double pTotal;
    
    @Override
    public void setInitParams()
    {
        //Set params to initial values and determine their order in params list
        betaChronic = 1.2943e-04; params.add(betaChronic); estParams.add(0);
        betaEarly = 20 * betaChronic; params.add(betaEarly); estParams.add(1);
        betaAIDS = 5 * betaChronic; params.add(betaAIDS); estParams.add(2);
        alpha = 0.0; params.add(alpha); fixedParams.add(3); //was 0.00015
        //betaChronic = 1.98e-04; params.add(betaChronic); estParams.add(0);
        //betaEarly = 0.0045; params.add(betaEarly); estParams.add(1);
        //betaAIDS = 0.001; params.add(betaAIDS); estParams.add(2);
        //alpha = 0.00017883; params.add(alpha); estParams.add(3);
        
        //Fixed params
        //double timeScaler = 26.0;
        gammaEarly = (1.0/(2.55 * 365.25)); params.add(gammaEarly); fixedParams.add(4);
        gammaChronic = (1.0/(6.31 * 365.25)); params.add(gammaChronic); fixedParams.add(5); //timeScaler * (1.0/(6.31 * 365.25));
        gammaAIDS = (1.0/(2.55 * 365.25)); params.add(gammaAIDS); fixedParams.add(6);
        mu = (1.0/(40.28*365.25)); params.add(mu); fixedParams.add(7);
        N = 1.0e04; params.add(N); fixedParams.add(8);
        tau = 0.0; params.add(tau); fixedParams.add(9);
        
        //Init conditions
        //initE = 1.0; params.add(initE); estParams.add(10);
        
        //Init introduction time
        //initTime = 720990.0 + 365.0; params.add(initTime); estParams.add(10); //Jan-01-1974
        initTime = 0.0; params.add(initTime); fixedParams.add(10); //Jan-01-1974
                        
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
        betaEarly = params.getQuick(1);
        betaChronic = params.getQuick(0);
        betaAIDS = params.getQuick(2);
        //alpha = params.getQuick(3);
        //initE = params.getQuick(10);
        //initTime = params.getQuick(10);
        
    }
    
    @Override
    public DoubleArrayList getInitialConditions() 
    {
        
        //Order init conditions 1) S 2) Ia 3) Ic 4) N
        currInits = new DoubleArrayList();

        //Compute equilibrium conditons without migration
        double initS = 10000.0 - 1; currInits.add(initS);
        initE = 1.0; currInits.add(initE);
        double initC = 0.0; currInits.add(initC);
        double initA = 0.0; currInits.add(initA);
        double initN = 10000.0; currInits.add(initN);
        
        //Flows
        //currInits.add(0.0); //transEE
        //currInits.add(0.0); //transCE
        //currInits.add(0.0); //transAE
        //currInits.add(0.0); //dEC
        //currInits.add(0.0); //dCA
        //currInits.add(0.0); //dBS
        //currInits.add(0.0); //dSD
        //currInits.add(0.0); //dED
        //currInits.add(0.0); //dCD
        //currInits.add(0.0); //dAD
        
        currInits.add(0.0); //cumulative incidence (index 15)
        currInits.add(0.0); //cumulative incidence in early class
        
        return currInits;
        
    }
    
    public DoubleArrayList getRandomInitialConditions() 
    {
        
        return currInits;
        
    }
    
    public double computeInitProbs(DoubleMatrix1D inits)
    {
        return 1.0;
        
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
        
        //720625 is Jan-01-1973
        //722267.0 is the rootTime for tree1
        //722892.0 is the rootTime for tree2
        //722866.0 is the rootTime for tree3
        //721980.0 is the rootTime for tree4
        //721203.0 is the rootTime for tree5
        //721596.0 is the rootTime for tree6
        //722212.0 is the rootTime for tree7;
        //723821.0 is the rootTime for tree8
        //723995.0 is the rootTime for tree9;
        //721819.0 is the rootTime for tree10;
        //if (params.get(10) < 720625 || params.get(10) > 722892.0) {
            //constraintCheckFail = true;
            //System.out.println("InitTime proposal out of bounds!");
        //}
        
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
        
        //Order of state variables is 1) S 2) E 3) C 4) A
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            int nextIndex = xIndex + 1;
            double dt = xDtTimes.getQuick(n);
            double tauScaler = 1.0;
            
            //double timeNow = xSubTimes.get(n);
            if (xSubTimes.get(n) >= 729756) { //if time is >= 1998
                tauScaler = 1 - tau;
            }
            
            for (int j = 0; j < jParticles; ++j) {
                
                if (xIndex >= startTimeIndex) { //for delayed start with set initial conditions

                pTotal = x.getQuick(4,j,xIndex);
                
                transScaler = Math.exp(-alpha * (x.getQuick(1,j,xIndex) + x.getQuick(2,j,xIndex) + x.getQuick(3,j,xIndex)));

                //INFECTIONS FROM EALRY CASES (INDEX 1)
                transEE = tauScaler * transScaler * betaEarly * x.getQuick(0,j,xIndex) * x.getQuick(1,j,xIndex) * dt / pTotal;

                //INFECTIONS FROM CHRONIC CASES (INDEX 2)
                transCE = tauScaler * transScaler * betaChronic * x.getQuick(0,j,xIndex) * x.getQuick(2,j,xIndex) * dt / pTotal;

                //INFECTIONS FROM AIDS CASES (INDEX 3)
                transAE = tauScaler * transScaler * betaAIDS * x.getQuick(0,j,xIndex) * x.getQuick(3,j,xIndex) * dt / pTotal;

                //Total incidence
                dSE = transEE + transCE + transAE;
                if (Double.isNaN(dSE)) {
                    dSE = 0.0;
                }
                
                //DISEASE PROGRESSION
                dEC = tauScaler * gammaEarly * x.getQuick(1,j,xIndex) * dt;
                dCA = tauScaler * gammaChronic * x.getQuick(2,j,xIndex) * dt;
                
                //BIRTHS
                dBS = mu * x.getQuick(4,j,xIndex) * dt;

                //Deaths
                dSD = mu * x.getQuick(0,j,xIndex) * dt;
                dED = mu * x.getQuick(1,j,xIndex) * dt;
                dCD = mu * x.getQuick(2,j,xIndex) * dt;
                dAD = tauScaler * gammaAIDS * x.getQuick(3,j,xIndex) * dt;

                newS = x.getQuick(0,j,xIndex) + dBS - dSE - dSD;
                newE = x.getQuick(1,j,xIndex) + dSE - dED - dEC;
                newC = x.getQuick(2,j,xIndex) + dEC - dCD - dCA;
                newA = x.getQuick(3,j,xIndex) + dCA - dAD;
                newP = newS + newE + newC + newA;

                x.setQuick(0,j,nextIndex,newS);
                x.setQuick(1,j,nextIndex,newE);
                x.setQuick(2,j,nextIndex,newC);
                x.setQuick(3,j,nextIndex,newA);
                x.setQuick(4,j,nextIndex,newP);
                //x.setQuick(5,j,nextIndex,transEE);
                //x.setQuick(6,j,nextIndex,transCE);
                //x.setQuick(7,j,nextIndex,transAE);
                //x.setQuick(8,j,nextIndex,dEC);
                //x.setQuick(9,j,nextIndex,dCA);
                //x.setQuick(10,j,nextIndex,dBS);
                //x.setQuick(11,j,nextIndex,dSD);
                //x.setQuick(12,j,nextIndex,dED);
                //x.setQuick(13,j,nextIndex,dCD);
                //x.setQuick(14,j,nextIndex,dAD);
                                
                //Update incidence
                cumulativeI = x.getQuick(5,j,xIndex) + (transEE + transCE + transAE);
                x.setQuick(5,j,nextIndex,cumulativeI);
                double cumulativeIE = x.getQuick(6, j, xIndex) + transEE;
                x.setQuick(6,j,nextIndex,cumulativeIE);
                
                
                }//if >= startTimeIndex

            }	
        }
        return x;
    }
    
    @Override
    public DoubleMatrix3D updateStatesTauLeap(int jParticles, DoubleMatrix3D x, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
    {
        
        //Order of state variables is 1) S 2) E 3) C 4) A
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            int nextIndex = xIndex + 1;
            double dt = xDtTimes.getQuick(n);
            
            //double timeNow = xSubTimes.get(n);
            
            for (int j = 0; j < jParticles; ++j) {
                
                if (xIndex >= startTimeIndex) { //for delayed start with set initial conditions
                
                pTotal = x.getQuick(4,j,xIndex);
                
                transScaler = Math.exp(-alpha * (x.getQuick(1,j,xIndex) + x.getQuick(2,j,xIndex) + x.getQuick(3,j,xIndex)));

                //INFECTIONS FROM EALRY CASES (INDEX 1)
                transEE = poissDist.nextInt(transScaler * betaEarly * x.getQuick(0,j,xIndex) * x.getQuick(1,j,xIndex) * dt / pTotal);

                //INFECTIONS FROM CHRONIC CASES (INDEX 2)
                transCE = poissDist.nextInt(transScaler * betaChronic * x.getQuick(0,j,xIndex) * x.getQuick(2,j,xIndex) * dt / pTotal);

                //INFECTIONS FROM AIDS CASES (INDEX 3)
                transAE = poissDist.nextInt(transScaler * betaAIDS * x.getQuick(0,j,xIndex) * x.getQuick(3,j,xIndex) * dt / pTotal);

                //Total incidence
                dSE = transEE + transCE + transAE;
                
                //DISEASE PROGRESSION
                dEC = poissDist.nextInt(gammaEarly * x.getQuick(1,j,xIndex) * dt);
                dCA = poissDist.nextInt(gammaChronic * x.getQuick(2,j,xIndex) * dt);
                
                //BIRTHS
                dBS = poissDist.nextInt(mu * x.getQuick(4,j,xIndex) * dt);

                //Deaths
                dSD = poissDist.nextInt(mu * x.getQuick(0,j,xIndex) * dt);
                dED = poissDist.nextInt(mu * x.getQuick(1,j,xIndex) * dt);
                dCD = poissDist.nextInt(mu * x.getQuick(2,j,xIndex) * dt);
                dAD = poissDist.nextInt(gammaAIDS * x.getQuick(3,j,xIndex) * dt);

                newS = x.getQuick(0,j,xIndex) + dBS - dSE - dSD;
                if (newS < 0.0) {
                    newS = 0.0;
                }
                newE = x.getQuick(1,j,xIndex) + dSE - dED - dEC;
                if (newE < 0.0) {
                    newE = 0.0;
                }
                newC = x.getQuick(2,j,xIndex) + dEC - dCD - dCA;
                if (newC < 0.0) {
                    newC = 0.0;
                }
                newA = x.getQuick(3,j,xIndex) + dCA - dAD;
                if (newA < 0.0) {
                    newA = 0.0;
                }
                newP = newS + newE + newC + newA;
                
                //if (newE == 0.0 & newC == 0.0 & newA == 0.0) {
                    //System.out.println("Particle" + j + " is dead");
                //}
                
                x.setQuick(0,j,nextIndex,newS);
                x.setQuick(1,j,nextIndex,newE);
                x.setQuick(2,j,nextIndex,newC);
                x.setQuick(3,j,nextIndex,newA);
                x.setQuick(4,j,nextIndex,newP);
                //x.setQuick(5,j,nextIndex,transEE);
                //x.setQuick(6,j,nextIndex,transCE);
                //x.setQuick(7,j,nextIndex,transAE);
                //x.setQuick(8,j,nextIndex,dEC);
                //x.setQuick(9,j,nextIndex,dCA);
                //x.setQuick(10,j,nextIndex,dBS);
                //x.setQuick(11,j,nextIndex,dSD);
                //x.setQuick(12,j,nextIndex,dED);
                //x.setQuick(13,j,nextIndex,dCD);
                //x.setQuick(14,j,nextIndex,dAD);
                
                //Update incidence
                cumulativeI = x.getQuick(5,j,xIndex) + (transEE + transCE + transAE);
                x.setQuick(5,j,nextIndex,cumulativeI);
                double cumulativeIE = x.getQuick(6, j, xIndex) + transEE;
                x.setQuick(6,j,nextIndex,cumulativeIE);
                
                }
                
            }	
        }
        return x;
    }
    
    public DoubleMatrix3D updateStatesTauLeapSingleParticle(int particle, DoubleMatrix3D x, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
    {
        
        //Order of state variables is 1) S 2) E 3) C 4) A
        int j = particle;
        alive = true;
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            int nextIndex = xIndex + 1;
            double dt = xDtTimes.getQuick(n);
            
            pTotal = x.getQuick(4,j,xIndex);

            transScaler = Math.exp(-alpha * (x.getQuick(1,j,xIndex) + x.getQuick(2,j,xIndex) + x.getQuick(3,j,xIndex)));

            //INFECTIONS FROM EALRY CASES (INDEX 1)
            transEE = poissDist.nextInt(transScaler * betaEarly * x.getQuick(0,j,xIndex) * x.getQuick(1,j,xIndex) * dt / pTotal);

            //INFECTIONS FROM CHRONIC CASES (INDEX 2)
            transCE = poissDist.nextInt(transScaler * betaChronic * x.getQuick(0,j,xIndex) * x.getQuick(2,j,xIndex) * dt / pTotal);

            //INFECTIONS FROM AIDS CASES (INDEX 3)
            transAE = poissDist.nextInt(transScaler * betaAIDS * x.getQuick(0,j,xIndex) * x.getQuick(3,j,xIndex) * dt / pTotal);

            //Total incidence
            dSE = transEE + transCE + transAE;

            //DISEASE PROGRESSION
            dEC = poissDist.nextInt(gammaEarly * x.getQuick(1,j,xIndex) * dt);
            dCA = poissDist.nextInt(gammaChronic * x.getQuick(2,j,xIndex) * dt);

            //BIRTHS
            dBS = poissDist.nextInt(mu * x.getQuick(4,j,xIndex) * dt);

            //Deaths
            dSD = poissDist.nextInt(mu * x.getQuick(0,j,xIndex) * dt);
            dED = poissDist.nextInt(mu * x.getQuick(1,j,xIndex) * dt);
            dCD = poissDist.nextInt(mu * x.getQuick(2,j,xIndex) * dt);
            dAD = poissDist.nextInt(gammaAIDS * x.getQuick(3,j,xIndex) * dt);

            newS = x.getQuick(0,j,xIndex) + dBS - dSE - dSD;
            //if (newS < 0.0) {
              //  newS = 0.0;
            //}
            newE = x.getQuick(1,j,xIndex) + dSE - dED - dEC;
            //if (newE < 0.0) {
              //  newE = 0.0;
            //}
            newC = x.getQuick(2,j,xIndex) + dEC - dCD - dCA;
            //if (newC < 0.0) {
              //  newC = 0.0;
            //}
            newA = x.getQuick(3,j,xIndex) + dCA - dAD;
            //if (newA < 0.0) {
              //  newA = 0.0;
            //}
            newP = newS + newE + newC + newA;

            if (newE < 1.0 & newC < 1.0 & newA < 1.0) {
                alive = false;
                System.out.println("Particle" + j + " is dead");
                break;
            }

            x.setQuick(0,j,nextIndex,newS);
            x.setQuick(1,j,nextIndex,newE);
            x.setQuick(2,j,nextIndex,newC);
            x.setQuick(3,j,nextIndex,newA);
            x.setQuick(4,j,nextIndex,newP);
            //x.setQuick(5,j,nextIndex,transSIa);
            //x.setQuick(6,j,nextIndex,transSIc);
            //x.setQuick(7,j,nextIndex,dIaIc);
            //x.setQuick(8,j,nextIndex,dIcR);

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
        
        pTotal = xTimeNow.getQuick(4);
        transScaler = Math.exp(-alpha * (xTimeNow.getQuick(1) + xTimeNow.getQuick(2) + xTimeNow.getQuick(3)));
        transEE = transScaler * betaEarly * xTimeNow.getQuick(0) * xTimeNow.getQuick(1) * dt / pTotal;
        transCE = transScaler * betaChronic * xTimeNow.getQuick(0) * xTimeNow.getQuick(2) * dt / pTotal;
        transAE = transScaler * betaAIDS * xTimeNow.getQuick(0) * xTimeNow.getQuick(3) * dt / pTotal;
        dEC = gammaEarly * xTimeNow.getQuick(1) * dt;
        dCA = gammaChronic * xTimeNow.getQuick(2) * dt;
        dBS = mu * xTimeNow.getQuick(4) * dt;
        dSD = mu * xTimeNow.getQuick(0) * dt;
        dED = mu * xTimeNow.getQuick(1) * dt;
        dCD = mu * xTimeNow.getQuick(2) * dt;
        dAD = gammaAIDS * xTimeNow.getQuick(3) * dt;
        
        //Implicitly sets transProb to 1.0 if mean for poisson density function is 0.0, which is technically undefined
        if (transEE > 0.0) {
            poissDist.setMean(transEE);
            transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(5));
        } else {
            if (xTimePlusOne.getQuick(5) > 0.0) {
                transProb *= 0.0;
            }
        }
        if (transCE > 0.0) {
            poissDist.setMean(transCE);
            transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(6));
        } else {
            if (xTimePlusOne.getQuick(6) > 0.0) {
                transProb *= 0.0;
            }
        }
        if (transAE > 0.0) {
            poissDist.setMean(transAE);
            transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(7));
        } else {
            if (xTimePlusOne.getQuick(7) > 0.0) {
                transProb *= 0.0;
            }
        }
        if (dEC > 0.0) {
            poissDist.setMean(dEC);
            transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(8));
        } else {
            if (xTimePlusOne.getQuick(8) > 0.0) {
                transProb *= 0.0;
            }
        }
        if (dCA > 0.0) {
            poissDist.setMean(dCA);
            transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(9));
        } else {
            if (xTimePlusOne.getQuick(9) > 0.0) {
                transProb *= 0.0;
            }
        }
        if (dBS > 0.0) {
            poissDist.setMean(dBS);
            transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(10));
        } else {
            if (xTimePlusOne.getQuick(10) > 0.0) {
                transProb *= 0.0;
            }
        }
        if (dSD > 0.0) {
            poissDist.setMean(dSD);
            transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(11));
        } else {
            if (xTimePlusOne.getQuick(11) > 0.0) {
                transProb *= 0.0;
            }
        }
        if (dED > 0.0) {
            poissDist.setMean(dED);
            transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(12));
        } else {
            if (xTimePlusOne.getQuick(12) > 0.0) {
                transProb *= 0.0;
            }
        }
        if (dCD > 0.0) {
            poissDist.setMean(dCD);
            transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(13));
        } else {
            if (xTimePlusOne.getQuick(13) > 0.0) {
                transProb *= 0.0;
            }
        }
        if (dAD > 0.0) {
            poissDist.setMean(dAD);
            transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(14));
        } else {
            if (xTimePlusOne.getQuick(14) > 0.0) {
                transProb *= 0.0;
            }
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