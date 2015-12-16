import cern.colt.list.*;
import java.util.ArrayList;
import cern.jet.random.AbstractDistribution;
import cern.colt.matrix.DoubleMatrix3D;

/**
 * Specifies the epidemiological model
 * Semi spatial has two coupled subpopulations but the global population is assumed to be in endemic equilibrium
 * Migration can either be reversible or non-reversible
 * 
 * @author David
 */
public class EpiModelSemiSpatial {
    

    DoubleArrayList params = new DoubleArrayList();
    ArrayList<Integer> estParams = new ArrayList<Integer>(); //indexes of the estimated params
    ArrayList<Integer> fixedParams = new ArrayList<Integer>(); //indexes of the fixed params
    
    public void setInitParams()
    {
        //Set params to initial values and determine their order in params list
        double beta = 1.0; params.add(beta); estParams.add(0);
        double alpha = 0.256218; params.add(alpha); estParams.add(1);
        double delta = 155.46; params.add(delta); estParams.add(2);
        double M = 0.004292; params.add(M); estParams.add(3);
        double IgEff = 4826.0; params.add(IgEff); estParams.add(4);
        double NFocal = 8.0e06; params.add(NFocal); fixedParams.add(5);
        double mu = 1/(70*365.25); params.add(mu); fixedParams.add(6);
        double nu = 0.2; params.add(nu); fixedParams.add(7);
        double fNoise = 0.0; params.add(fNoise); fixedParams.add(8);
        double MSwitch = 731582; params.add(MSwitch); fixedParams.add(9);
                        
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
        
        //Order init conditions 1) Sf 2) If 3) Nf 4) IgEff
        
        DoubleArrayList currInits = new DoubleArrayList();
        
        double beta = params.get(0); 
        double alpha = params.get(1);
        double M = params.get(3);
        double IgEff = params.get(4);
        double Nf = params.get(5);
        double mu = params.get(6);
        double nu = params.get(7);
        double averageBeta = beta * (1 + alpha * 0.3183); //0.3183 is sin(pi/2) / pi;
        double R0 = averageBeta / (mu + nu);
        
        //Compute equilbrium S and I values for focal population
        //Have to solve the following quadratic to compute equilibrium values
        double aConst = (averageBeta * mu + averageBeta * nu);
        double bConst = (averageBeta * IgEff * M * mu) - (averageBeta * mu * Nf) + (mu * mu *  Nf) - (averageBeta * IgEff * M * nu) + (mu * nu * Nf);
        double cConst = averageBeta * IgEff * M * mu * Nf;
        double IfHat = -bConst + Math.sqrt(bConst*bConst - 4*aConst*cConst) / (2*aConst); //solve quadratic for equlibrium IfHat
        if (IfHat < 0) {
            System.out.println("Init If is a negative number");
        }
        double SfHat = ((mu + nu) * IfHat) / (averageBeta*(IfHat/Nf) + averageBeta*M*(IgEff/Nf)); 
        if (SfHat < 0) {
            System.out.println("Init Sf is a negative number");
        }
        
        currInits.add(SfHat);
        currInits.add(IfHat);
        currInits.add(Nf);
        currInits.add(IgEff);
        
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
        
        //R0 has to be greater than one and less than 30
        double averageBeta = params.get(0) * (1 + params.get(1) * 0.3183);
        double R0 = averageBeta / (0.2 + (1/(70*365.25)));
        if (R0 < 1 | R0 > 30) {
            constraintCheckFail = true;
        }
        
        //M has to be less than one
        if (params.get(3) > 1.0) {
            constraintCheckFail = true;
        }
        
        //deltaFocal has to be less than 365.25
        if (params.get(2) > (365.25)) { //this means delta is less than one year
            constraintCheckFail = true;
        }
        
        return constraintCheckFail;
    }
    
    
    public DoubleMatrix3D updateStates(int jParticles, DoubleMatrix3D matrix, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes, AbstractDistribution stndNorm) 
    {
      
        double eta; double dSfIf; double dSfIg; double dSgIg; double dIfRf; double dIgRg;
        double dSfD; double dIfD; double dRfD; double dSgD; double dIgD; double dRgD; double dNfB; double dNgB; double Rf;
        double newSf; double newIf; double newRf; double newNf; double newIgEff;
        double fTerm;
		
        double beta = params.getQuick(0);
	double alpha = params.getQuick(1);
        double delta = params.getQuick(2);
        double M = params.getQuick(3);
        double IgEff = params.getQuick(4);
        double mu = params.get(6);
        double nu = params.get(7);
	double fNoise = params.get(8);
        double MSwitchPoint = params.get(9);
	double stndNormRand;
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            double dt = xDtTimes.getQuick(n);
            
            double timeNow = xSubTimes.get(n);
            
            //For piece-wise seasonal forcing with a high and low season
            //double currTimeYear = (timeNow / 365.25) - Math.floor(timeNow/365.25);
            //double startHighSeason = delta / 365.25;
            //double endHighSeason = (delta + (365.25/2)) / 365.25;
            //double betaSeas = 0.0;
            //if (endHighSeason < 1.0) {
                //if (currTimeYear > startHighSeason & currTimeYear < endHighSeason) {
                    //betaSeas = betaHigh;
                //} else {
                    //betaSeas = betaLow;
                //}
            //} else {
                //if (currTimeYear < endHighSeason | currTimeYear > startHighSeason) {
                    //betaSeas = betaHigh;
                //} else {
                    //betaSeas = betaLow;
                //}
            //}
            
            double seasNow = ((timeNow + delta)/365.25) - Math.floor(timeNow/365.25);
            double seasHeight = Math.cos(2*Math.PI * seasNow);
            if (seasHeight < 0) {
                seasHeight = 0;
            }
            double betaSeas = beta * (1 + alpha*seasHeight);
            
            for (int j = 0; j < jParticles; ++j) {
                
                stndNormRand = stndNorm.nextDouble();
                eta = stndNormRand / Math.sqrt(dt); //used for this particle and time only

                //INFECTIONS IN FOCAL POP FROM FOCAL POP
                if (timeNow >= MSwitchPoint) {
                    fTerm = fNoise * betaSeas * matrix.getQuick(0,j,xIndex) * matrix.getQuick(1,j,xIndex) * eta * dt / matrix.getQuick(2, j, xIndex);
                    dSfIf = (betaSeas * matrix.getQuick(0,j,xIndex) * matrix.getQuick(1,j,xIndex) * dt / matrix.getQuick(2,j,xIndex)) + fTerm;
                } else {
                    dSfIf = 0.0;
                }
                
                //INFECTIONS IN FOCAL POP FROM GLOBAL POP
                if (timeNow >= MSwitchPoint) {
                    dSfIg = betaSeas * M * matrix.getQuick(0,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(2,j,xIndex);
                } else {
                    dSfIg = 0.0;
                } 

                //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
                dSgIg = 0.0;
                
                //RECOVERIES
                dIfRf = nu * matrix.getQuick(1,j,xIndex) * dt;
                dIgRg = 0.0;
                
                //BIRTHS
                if (timeNow >= MSwitchPoint) {
                    dNfB = mu * matrix.getQuick(2,j,xIndex) * dt;
                } else {
                    dNfB = 0.0;
                }
                dNgB = 0.0;
                
                //DEATHS IN FOCAL POP
                if (timeNow >= MSwitchPoint) {
                    dSfD = mu * matrix.getQuick(0,j,xIndex) * dt;
                } else {
                    dSfD = 0.0;
                }
                dIfD = mu * matrix.getQuick(1,j,xIndex) * dt;
                Rf = matrix.getQuick(2,j,xIndex) - matrix.getQuick(0,j,xIndex) - matrix.getQuick(1,j,xIndex);
                if (timeNow >= MSwitchPoint) {
                    dRfD = mu * Rf * dt;
                } else {
                    dRfD = 0.0;
                }
                
                //DEATHS IN GLOBAL POP
                dSgD = 0.0;
                dIgD = 0.0;
                dRgD = 0.0;

                newSf = matrix.getQuick(0,j,xIndex) + dNfB - dSfIf - dSfIg - dSfD;
                newIf = matrix.getQuick(1,j,xIndex) + dSfIf + dSfIg - dIfRf - dIfD;
                newRf = Rf + dIfRf - dRfD;
                newNf = newSf + newIf + newRf;
                newIgEff = matrix.getQuick(3,j,xIndex);

                int nextIndex = xIndex + 1;
                matrix.setQuick(0,j,nextIndex,newSf);
                matrix.setQuick(1,j,nextIndex,newIf);
                matrix.setQuick(2,j,nextIndex,newNf);
                matrix.setQuick(3,j,nextIndex,newIgEff);
            }	
        }
        return matrix;
    }
}