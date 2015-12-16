import cern.colt.list.*;
import java.util.ArrayList;
import cern.jet.random.AbstractDistribution;
import cern.colt.matrix.DoubleMatrix3D;

/**
 * Specifies the epidemiological model
 * @author David
 */
public class EpiModelFullSpatial {

    DoubleArrayList params = new DoubleArrayList();
    ArrayList<Integer> estParams = new ArrayList<Integer>(); //indexes of the estimated params
    ArrayList<Integer> fixedParams = new ArrayList<Integer>(); //indexes of the fixed params
    
    public void setInitParams()
    {
        //Set params to initial values and determine their order in params list
        double beta = 1.000; params.add(beta); estParams.add(0);
        double alphaFocal = 0.05; params.add(alphaFocal); estParams.add(1);
        double alphaGlobal = 0.14; params.add(alphaGlobal); fixedParams.add(2);
        double deltaFocal = 6.0/12.0; params.add(deltaFocal); estParams.add(3);
        double deltaGlobal = 0.0/12.0; params.add(deltaGlobal); fixedParams.add(4);
        double M = 1.0/1000.0; params.add(M); fixedParams.add(5); // assumes symmetric import and export rates
        double NFocal = 5.0e06; params.add(NFocal); fixedParams.add(6);
        double NGlobal = 5.0e06; params.add(NGlobal); fixedParams.add(7);
        double mu = 1/(70*365.2); params.add(mu); fixedParams.add(8);
        double nu = 0.2; params.add(nu); fixedParams.add(9);
        double fNoise = 0.0012; params.add(fNoise); fixedParams.add(10);
        double MSwitch = 0.0; params.add(MSwitch); fixedParams.add(11);
                        
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
            params.set(paramListIndex, newParams.get(i));
        }   
        
    }
    
    public DoubleArrayList getInitialConditions() 
    {
        
        //Order init conditions 1) Sf 2) Sg 3) If 4) Ig 5) Nf 6) Ng
        
        DoubleArrayList currInits = new DoubleArrayList();
        
        //Compute betaAvg and R0 in order to find equilibrium conditions without migration
        double beta = params.get(0); 
        double alpha1 = params.get(1);
        double alpha2 = params.get(2);
        double N1 = params.get(6);
        double N2 = params.get(7);
        double mu = params.get(8);
        double nu = params.get(9);
        double averageBeta1 = beta * (1 + alpha1 * 0.3183); //0.3183 is sin(pi/2) / pi;
        double averageBeta2 = beta * (1 + alpha2 * 0.3183);
        double R01 = averageBeta1 / (mu + nu); //assuming no migration
        double R02 = averageBeta2 / (mu + nu); //assuming no migration
        
        //Compute equilibrium conditons without migration
        double initS1 = N1/R01; currInits.add(initS1);
        double initS2 = N2/R02; currInits.add(initS2);
        double initI1 = mu * N1 * (R01 - 1) / averageBeta1; currInits.add(initI1);
        double initI2 = mu * N2 * (R02 - 1) / averageBeta2; currInits.add(initI2);
        currInits.add(N1);
        currInits.add(N2);
        
        double dt = 0.50;
        double length = 365.25*100;
        currInits = this.fastSim(currInits, dt, length);
        
        return currInits;
        
    }
    
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
        
        //R01 has to be greater than one and less than 30
        double averageBeta1 = params.get(0) * (1 + params.get(1) * 0.3183);
        double R01 = averageBeta1 / (0.2 + (1/(70*365.25)));
        if (R01 < 1 | R01 > 30) {
            constraintCheckFail = true;
        }
        
        //M has to be less than one
        if (params.get(5) > 1.0) {
            constraintCheckFail = true;
        }
        
        //deltaFocal has to be less than 365.25
        if (params.get(3) > (1.0)) { //this means delta is less than one year
            constraintCheckFail = true;
        }
        
        //R01 has to be greater than one and less than 30
        double averageBeta2 = params.get(0) * (1 + params.get(2) * 0.3183);
        double R02 = averageBeta2 / (0.2 + (1/(70*365.25)));
        if (R02 < 1 | R02 > 30) {
            constraintCheckFail = true;
        }
        
        //deltaGlobal has to be less than 365.25
        if (params.get(4) > (1.0)) { //this means delta is less than one year
            constraintCheckFail = true;
        }
        
        
        return constraintCheckFail;
    }
    
    
    public DoubleMatrix3D updateStates(int jParticles, DoubleMatrix3D matrix, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes, AbstractDistribution stndNorm) 
    {
        
        //Spatial SIR with focal and global population
        //Order of state variables is 1) Sf 2) Sg 3) If 4) Ig 5) Nf 6) Ng
      
        double eta; double dSfIf; double dSfIg; double dSgIg; double dSgIf; double dIfRf; double dIgRg;
        double dSfD; double dIfD; double dRfD; double dSgD; double dIgD; double dRgD; double dNfB; double dNgB; double Rf; double Rg;
        double newSf; double newSg; double newIf; double newIg; double newRf;  double newRg; double newNf; double newNg;
        double fTerm;
		
        double beta = params.getQuick(0);
	double alphaFocal = params.getQuick(1);
        double alphaGlobal = params.getQuick(2);
        double deltaFocal = params.getQuick(3) * 365.25;
        double deltaGlobal = params.getQuick(4) * 365.25;
        double M = params.getQuick(5);
        double mu = params.getQuick(8);
        double nu = params.getQuick(9);
        double fNoise = params.getQuick(10);
        double MSwitchPoint = params.getQuick(11);
	double stndNormRand;
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            double dt = xDtTimes.getQuick(n);
            
            double timeNow = xSubTimes.get(n);
            
            //Compute betaSeas for focal pop
            double seasNow = ((timeNow + deltaFocal)/365.25) - Math.floor(timeNow/365.25);
            double seasHeight = Math.cos(2*Math.PI * seasNow);
            if (seasHeight < 0) {
                seasHeight = 0;
            }
            double betaSeasFocal = beta * (1 + alphaFocal*seasHeight);
            
            //Compute betaSeas for global pop
            seasNow = ((timeNow + deltaGlobal)/365.25) - Math.floor(timeNow/365.25);
            seasHeight = Math.cos(2*Math.PI * seasNow);
            if (seasHeight < 0) {
                seasHeight = 0;
            }
            double betaSeasGlobal = beta * (1 + alphaGlobal*seasHeight);
            
            for (int j = 0; j < jParticles; ++j) {
                
                stndNormRand = stndNorm.nextDouble();
                eta = stndNormRand / Math.sqrt(Math.abs(dt)); //used for this particle and time only

                //INFECTIONS IN FOCAL POP FROM FOCAL POP
                if (timeNow >= MSwitchPoint) {
                    fTerm = fNoise * betaSeasFocal * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * eta * dt / matrix.getQuick(4, j, xIndex);
                    dSfIf = (betaSeasFocal * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(4,j,xIndex)) + fTerm;
                } else {
                    dSfIf = 0.0;
                }
                
                //INFECTIONS IN FOCAL POP FROM GLOBAL POP
                if (timeNow >= MSwitchPoint) {
                    dSfIg = betaSeasFocal * M * matrix.getQuick(0,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(4,j,xIndex);
                } else {
                    dSfIg = 0.0;
                } 

                //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
                if (timeNow >= MSwitchPoint) {
                    fTerm = fNoise * betaSeasGlobal * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * eta * dt / matrix.getQuick(5, j, xIndex);
                    dSgIg = (betaSeasGlobal * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(5,j,xIndex)) + fTerm;
                } else {
                    dSgIg = 0.0;
                }
                
                //INFECTIONS IN GLOBAL POP FROM FOCAL POP
                if (timeNow >= MSwitchPoint) {
                    dSgIf = betaSeasGlobal * M * matrix.getQuick(1,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(5,j,xIndex);
                } else {
                    dSgIf = 0.0;
                } 
                
                //RECOVERIES
                if (timeNow >= MSwitchPoint) {
                    dIfRf = nu * matrix.getQuick(2,j,xIndex) * dt;
                    dIgRg = nu * matrix.getQuick(3,j,xIndex) * dt;
                } else {
                    dIfRf = 0.0;
                    dIgRg = 0.0;
                }
                
                //BIRTHS
                if (timeNow >= MSwitchPoint) {
                    dNfB = mu * matrix.getQuick(4,j,xIndex) * dt;
                    dNgB = mu * matrix.getQuick(5,j,xIndex) * dt;
                } else {
                    dNfB = 0.0;
                    dNgB = 0.0;
                }
                
                //DEATHS IN FOCAL POP
                Rf = matrix.getQuick(4,j,xIndex) - matrix.getQuick(0,j,xIndex) - matrix.getQuick(2,j,xIndex);
                Rg = matrix.getQuick(5,j,xIndex) - matrix.getQuick(1,j,xIndex) - matrix.getQuick(3,j,xIndex);
                if (timeNow >= MSwitchPoint) {
                    dSfD = mu * matrix.getQuick(0,j,xIndex) * dt;
                    dIfD = mu * matrix.getQuick(2,j,xIndex) * dt;
                    dRfD = mu * Rf * dt;
                    
                    dSgD = mu * matrix.getQuick(1,j,xIndex) * dt;
                    dIgD = mu * matrix.getQuick(3,j,xIndex) * dt;
                    dRgD = mu * Rg * dt;   
                    
                } else {
                    dSfD = 0.0;
                    dIfD = 0.0;
                    dRfD = 0.0;
                    dSgD = 0.0;
                    dIgD = 0.0;
                    dRgD = 0.0;
                }

                newSf = matrix.getQuick(0,j,xIndex) + dNfB - dSfIf - dSfIg - dSfD;
                newSg = matrix.getQuick(1,j,xIndex) + dNgB - dSgIg - dSgIf - dSgD;
                newIf = matrix.getQuick(2,j,xIndex) + dSfIf + dSfIg - dIfRf - dIfD;
                newIg = matrix.getQuick(3,j,xIndex) + dSgIg + dSgIf - dIgRg - dIgD;
                newRf = Rf + dIfRf - dRfD;
                newRg = Rg + dIgRg - dRgD;
                newNf = newSf + newIf + newRf;
                newNg = newSg + newIg + newRg;

                int nextIndex = xIndex - 1;
                matrix.setQuick(0,j,nextIndex,newSf);
                matrix.setQuick(1,j,nextIndex,newSg);
                matrix.setQuick(2,j,nextIndex,newIf);
                matrix.setQuick(3,j,nextIndex,newIg);
                matrix.setQuick(4,j,nextIndex,newNf);
                matrix.setQuick(5,j,nextIndex,newNg);
            }	
        }
        return matrix;
    }
    
    public DoubleArrayList fastSim(DoubleArrayList inits, double dt, double endTime)
    {
        double dSfIf; double dSfIg; double dSgIg; double dSgIf; double dIfRf; double dIgRg;
        double dSfD; double dIfD; double dRfD; double dSgD; double dIgD; double dRgD; double dNfB; double dNgB; double Rf; double Rg;
        double newSf; double newSg; double newIf; double newIg; double newRf;  double newRg; double newNf; double newNg;
        //double fTerm;
		
        double beta = params.getQuick(0);
	double alphaFocal = params.getQuick(1);
        double alphaGlobal = params.getQuick(2);
        double deltaFocal = params.getQuick(3) * 365.25;
        double deltaGlobal = params.getQuick(4) * 365.25;
        double M = params.getQuick(5);
        double mu = params.getQuick(8);
        double nu = params.getQuick(9);
        double fNoise = params.getQuick(10);
        //double MSwitchPoint = params.getQuick(11);
	//double stndNormRand;
        
        //Get sim times
        DoubleArrayList times = new DoubleArrayList();
        double timeNow = 0.0;
        times.add(timeNow);
        while (timeNow < endTime) {
            double timeNew = timeNow + dt;
            if (timeNew <= endTime) {
                times.add(timeNew);
                timeNow = timeNew;
            } else {
                times.add(endTime);
                timeNow = timeNew;
            }
        }
        
        //Set up array for state variables
        DoubleArrayList x = new DoubleArrayList();
        x.add(inits.get(0));
        x.add(inits.get(1));
        x.add(inits.get(2));
        x.add(inits.get(3));
        x.add(inits.get(4));
        x.add(inits.get(5));
        
	for (int n = 0; n < times.size(); ++n) {
            
            timeNow = times.getQuick(n);
            
            //Compute betaSeas for focal pop
            double seasNow = ((timeNow + deltaFocal)/365.25) - Math.floor(timeNow/365.25);
            double seasHeight = Math.cos(2*Math.PI * seasNow);
            if (seasHeight < 0) {
                seasHeight = 0;
            }
            double betaSeasFocal = beta * (1 + alphaFocal*seasHeight);
            
            //Compute betaSeas for global pop
            seasNow = ((timeNow + deltaGlobal)/365.25) - Math.floor(timeNow/365.25);
            seasHeight = Math.cos(2*Math.PI * seasNow);
            if (seasHeight < 0) {
                seasHeight = 0;
            }
            double betaSeasGlobal = beta * (1 + alphaGlobal*seasHeight);
            
            
            //stndNormRand = stndNorm.nextDouble();
            //eta = stndNormRand / Math.sqrt(dt); //used for this particle and time only

            //INFECTIONS IN FOCAL POP FROM FOCAL POP
            //fTerm = fNoise * betaSeasFocal * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * eta * dt / matrix.getQuick(4, j, xIndex);
            dSfIf = (betaSeasFocal * x.getQuick(0) * x.getQuick(2) * dt / x.getQuick(4)); // + fTerm;

            //INFECTIONS IN FOCAL POP FROM GLOBAL POP
            dSfIg = betaSeasFocal * M * x.getQuick(0) * x.getQuick(3) * dt / x.getQuick(4);
 

            //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
            //fTerm = fNoise * betaSeasGlobal * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * eta * dt / matrix.getQuick(5, j, xIndex);
            dSgIg = (betaSeasGlobal * x.getQuick(1) * x.getQuick(3) * dt / x.getQuick(5)); //+ fTerm;

            //INFECTIONS IN GLOBAL POP FROM FOCAL POP
            dSgIf = betaSeasGlobal * M * x.getQuick(1) * x.getQuick(2) * dt / x.getQuick(5);

            //RECOVERIES
            dIfRf = nu * x.getQuick(2) * dt;
            dIgRg = nu * x.getQuick(3) * dt;

            //BIRTHS
            dNfB = mu * x.getQuick(4) * dt;
            dNgB = mu * x.getQuick(5) * dt;

            //DEATHS IN FOCAL POP
            Rf = x.getQuick(4) - x.getQuick(0) - x.getQuick(2);
            Rg = x.getQuick(5) - x.getQuick(1) - x.getQuick(3);
            dSfD = mu * x.getQuick(0) * dt;
            dIfD = mu * x.getQuick(2) * dt;
            dRfD = mu * Rf * dt;

            dSgD = mu * x.getQuick(1) * dt;
            dIgD = mu * x.getQuick(3) * dt;
            dRgD = mu * Rg * dt;   

            newSf = x.getQuick(0) + dNfB - dSfIf - dSfIg - dSfD;
            newSg = x.getQuick(1) + dNgB - dSgIg - dSgIf - dSgD;
            newIf = x.getQuick(2) + dSfIf + dSfIg - dIfRf - dIfD;
            newIg = x.getQuick(3) + dSgIg + dSgIf - dIgRg - dIgD;
            newRf = Rf + dIfRf - dRfD;
            newRg = Rg + dIgRg - dRgD;
            newNf = newSf + newIf + newRf;
            newNg = newSg + newIg + newRg;

            x.setQuick(0,newSf);
            x.setQuick(1,newSg);
            x.setQuick(2,newIf);
            x.setQuick(3,newIg);
            x.setQuick(4,newNf);
            x.setQuick(5,newNg);
	
        }
        return x;
    }
    
       
}
