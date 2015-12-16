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
public class EpiModel2Pop extends EpiModel {

    //DoubleArrayList params = new DoubleArrayList();
    //ArrayList<Integer> estParams = new ArrayList<Integer>(); //indexes of the estimated params
    //ArrayList<Integer> fixedParams = new ArrayList<Integer>(); //indexes of the fixed params
    //DoubleFactory1D factory1D = DoubleFactory1D.dense;
    //DoubleFactory2D factory2D = DoubleFactory2D.dense;
    
    boolean alive = true;
    Poisson poissDist; Normal stndNorm;
    
    double seasNowFocal; double seasNowGlobal; double seasHeightFocal; double seasHeightGlobal; double betaSeasFocal; double betaSeasGlobal;
    double beta; double alphaFocal; double alphaGlobal; double deltaFocal; double deltaGlobal;
    double M; double NFocal; double NGlobal; double mu; double nu; 
    double fNoise; double stndNormRand; double fTerm; double eta;
    
    double dSfIf; double dSfIg; double dSgIg; double dSgIf; double dIfRf; double dIgRg;
    double dSfD; double dIfD; double dRfD; double dSgD; double dIgD; double dRgD; double dNfB; double dNgB; double Rf; double Rg;
    double newSf; double newSg; double newIf; double newIg; double newRf;  double newRg; double newNf; double newNg;
    
    @Override
    public void setInitParams()
    {
        //Set params to initial values and determine their order in params list
        beta = 3.0 * ((1/(70*365.25)) + (1.0/7.0)); params.add(beta); estParams.add(0);
        alphaFocal = 0.08; params.add(alphaFocal); estParams.add(1);
        alphaGlobal = alphaFocal; params.add(alphaGlobal); fixedParams.add(2);
        deltaFocal = 6.0/12.0; params.add(deltaFocal); estParams.add(3);
        deltaGlobal = 0.0/12.0; params.add(deltaGlobal); estParams.add(4);
        M = 0.20; params.add(M); estParams.add(5); // assumes symmetric import and export rates
        NFocal = 2.0e06; params.add(NFocal); fixedParams.add(6);
        NGlobal = 2.0e06; params.add(NGlobal); fixedParams.add(7);
        mu = 1/(70*365.25); params.add(mu); fixedParams.add(8);
        nu = 1.0/7.0; params.add(nu); fixedParams.add(9);
        fNoise = 0.1; params.add(fNoise); fixedParams.add(10);
                        
    }
    
    @Override
    public void setRandomGenerator() 
    {
      
        cern.jet.random.engine.RandomEngine engine = new cern.jet.random.engine.MersenneTwister(new java.util.Date());
        poissDist = new Poisson(0.0, engine);
        stndNorm = new Normal(0.0, 1.0, engine);
        
    }
    
//    public DoubleArrayList getEstParams() 
//    {
//        DoubleArrayList currEstParams = new DoubleArrayList();
//        int p = estParams.size();
//        for (int i = 0; i < p; i++) {
//            int paramListIndex = estParams.get(i);
//            currEstParams.add(params.getQuick(paramListIndex));
//        }
//        return currEstParams;
//        
//    }
//    
//    public DoubleArrayList getFixedParams() 
//    {
//        DoubleArrayList currFixedParams = new DoubleArrayList();
//        int p = fixedParams.size();
//        int paramListIndex;
//        for (int i = 0; i < p; i++) {
//            paramListIndex = fixedParams.get(i);
//            currFixedParams.add(params.getQuick(paramListIndex));
//        }
//        return currFixedParams;
//        
//    }
    
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
        alphaFocal = params.getQuick(1);
        alphaGlobal = alphaFocal; params.setQuick(2,alphaFocal);
        deltaFocal = params.getQuick(3);
        deltaGlobal = params.getQuick(4);
        M = params.getQuick(5);
        
    }
    
    @Override
    public DoubleArrayList getInitialConditions() 
    {
        
        //Order init conditions 1) Sf 2) Sg 3) If 4) Ig 5) Nf 6) Ng
        
        DoubleArrayList currInits = new DoubleArrayList();
        
//        //Compute betaAvg and R0 in order to find equilibrium conditions without migration
//        double beta = params.get(0); 
//        double alpha1 = params.get(1);
//        double alpha2 = params.get(2);
//        double N1 = params.get(6);
//        double N2 = params.get(7);
//        double mu = params.get(8);
//        double nu = params.get(9);
//        double averageBeta1 = beta * (1 + alpha1 * 0.3183); //0.3183 is sin(pi/2) / pi;
//        double averageBeta2 = beta * (1 + alpha2 * 0.3183);
//        double R01 = averageBeta1 / (mu + nu); //assuming no migration
//        double R02 = averageBeta2 / (mu + nu); //assuming no migration
//        
//        //Compute equilibrium conditons without migration
//        double initS1 = N1/R01; currInits.add(initS1);
//        double initS2 = N2/R02; currInits.add(initS2);
//        double initI1 = mu * N1 * (R01 - 1) / averageBeta1; currInits.add(initI1);
//        double initI2 = mu * N2 * (R02 - 1) / averageBeta2; currInits.add(initI2);
//        currInits.add(N1);
//        currInits.add(N2);
//        
//        double dt = 0.50;
//        double length = 365.25*100;
//        currInits = this.fastSim(currInits, dt, length);
        
        //For 2PopStochRhoL tree2 at t = 5114
        //double initS1 = 665199.938; currInits.add(initS1);
        //double initS2 = 674417.425; currInits.add(initS2);
        //double initI1 = 462.519; currInits.add(initI1);
        //double initI2 = 607.76; currInits.add(initI2);
        //double initN1 = 1999805.0; currInits.add(initN1);
        //double initN2 = 2002141.0; currInits.add(initN2);
        
        //For 2PopStochRhoL tree2 at t = 0
        //double initS1 = 661580.384; currInits.add(initS1);
        //double initS2 = 665264.424; currInits.add(initS2);
        //double initI1 = 283.519; currInits.add(initI1);
        //double initI2 = 375.762; currInits.add(initI2);
        //double initN1 = 2000000.0; currInits.add(initN1);
        //double initN2 = 2000000.0; currInits.add(initN2);
        
        //For 2PopStochRhoM tree1 at t = 5114
        //double initS1 = 664876.265; currInits.add(initS1);
        //double initS2 = 661333.919; currInits.add(initS2);
        //double initI1 = 256.725; currInits.add(initI1);
        //double initI2 = 348.65; currInits.add(initI2);
        //double initN1 = 2000634.0; currInits.add(initN1);
        //double initN2 = 1998892.0; currInits.add(initN2);
        
        //For 2PopStochRhoM tree1 at t = 0
        //double initS1 = 658245.265; currInits.add(initS1);
        //double initS2 = 660180.919; currInits.add(initS2);
        //double initI1 = 233.725; currInits.add(initI1);
        //double initI2 = 468.648; currInits.add(initI2);
        //double initN1 = 2000000.0; currInits.add(initN1);
        //double initN2 = 2000000.0; currInits.add(initN2);
        
        //For 2PopStochRhoH tree1 at t = 5114
        //double initS1 = 651020.753; currInits.add(initS1);
        //double initS2 = 650050.444; currInits.add(initS2);
        //double initI1 = 236.89; currInits.add(initI1);
        //double initI2 = 333.115; currInits.add(initI2);
        //double initN1 = 2000456.0; currInits.add(initN1);
        //double initN2 = 1998877.0; currInits.add(initN2);
        
        //For 2PopStochRhoH tree1 at t = 0
        double initS1 = 663705.753; currInits.add(initS1);
        double initS2 = 663898.443; currInits.add(initS2);
        double initI1 = 297.890; currInits.add(initI1);
        double initI2 = 430.115; currInits.add(initI2);
        double initN1 = 2000000.0; currInits.add(initN1);
        double initN2 = 2000000.0; currInits.add(initN2);
        
        currInits.add(0.0); //dSfIf;
        currInits.add(0.0); //dSfIg;
        currInits.add(0.0); //dSgIg;
        currInits.add(0.0); //dSgIf;
        currInits.add(0.0); //IfRf;
        currInits.add(0.0); //IgRg;
        currInits.add(0.0); //NfB;
        currInits.add(0.0); //NgB;        
        currInits.add(0.0); //SfD;
        currInits.add(0.0); //IfD;
        currInits.add(0.0); //RfD;
        currInits.add(0.0); //SgD;
        currInits.add(0.0); //IgD;
        currInits.add(0.0); //RgD;
        
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
        //M has to be less than one
        if (params.get(5) > 1.0) {
            constraintCheckFail = true;
        }
        
        //deltaFocal has to be less than 365.25
        if (params.get(3) > (1.0)) { //this means delta is less than one year
            constraintCheckFail = true;
        }
        
        //deltaGlobal has to be less than 365.25
        if (params.get(4) > (1.0)) { //this means delta is less than one year
            constraintCheckFail = true;
        }
        
        
        return constraintCheckFail;
    }
    
    @Override
    public DoubleMatrix3D updateStatesTauLeap(int jParticles, DoubleMatrix3D matrix, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
    {
        
        //Spatial SIR with focal and global population
        //Order of state variables is 1) Sf 2) Sg 3) If 4) Ig 5) Nf 6) Ng
      
        //double dSfIf; double dSfIg; double dSgIg; double dSgIf; double dIfRf; double dIgRg;
        //double dSfD; double dIfD; double dRfD; double dSgD; double dIgD; double dRgD; double dNfB; double dNgB; double Rf; double Rg;
        //double newSf; double newSg; double newIf; double newIg; double newRf;  double newRg; double newNf; double newNg;
        //double dTerm; //double fTerm;
		
        //beta = params.getQuick(0);
	//alphaFocal = params.getQuick(1);
        //alphaGlobal = params.getQuick(2);
        //deltaFocal = params.getQuick(3) * 365.25;
        //deltaGlobal = params.getQuick(4) * 365.25;
        //M = params.getQuick(5);
        //mu = params.getQuick(8);
        //nu = params.getQuick(9);
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            double dt = xDtTimes.getQuick(n);
            int nextIndex = xIndex + 1;
            
            double timeNow = xSubTimes.get(n);
            
            //Compute betaSeas for focal pop
            seasNowFocal = ((timeNow/365.25) + deltaFocal) - Math.floor(timeNow/365.25);
            seasHeightFocal = Math.cos(2*Math.PI * seasNowFocal);
            betaSeasFocal = beta * (1 + alphaFocal*seasHeightFocal);
            
            //Compute betaSeas for global pop
            seasNowGlobal = ((timeNow/365.25) + deltaGlobal) - Math.floor(timeNow/365.25);
            seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
            betaSeasGlobal = beta * (1 + alphaGlobal*seasHeightGlobal);
            
            for (int j = 0; j < jParticles; ++j) {
                
                //double nScaler = 1.0;

                //INFECTIONS IN FOCAL POP FROM FOCAL POP
                dSfIf = poissDist.nextInt((1-M) * betaSeasFocal * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));
                
                //INFECTIONS IN FOCAL POP FROM GLOBAL POP
                dSfIg = poissDist.nextInt(betaSeasFocal * M * matrix.getQuick(0,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));

                //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
                dSgIg = poissDist.nextInt((1-M) * betaSeasGlobal * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(5,j,xIndex));

                //INFECTIONS IN GLOBAL POP FROM FOCAL POP
                dSgIf = poissDist.nextInt(betaSeasGlobal * M * matrix.getQuick(1,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(5,j,xIndex));
                
                //RECOVERIES
                dIfRf = poissDist.nextInt(nu * matrix.getQuick(2,j,xIndex) * dt);
                dIgRg = poissDist.nextInt(nu * matrix.getQuick(3,j,xIndex) * dt);
                
                //BIRTHS
                dNfB = poissDist.nextInt(mu * matrix.getQuick(4,j,xIndex) * dt);
                dNgB = poissDist.nextInt(mu * matrix.getQuick(5,j,xIndex) * dt);
                
                //DEATHS IN FOCAL POP
                Rf = matrix.getQuick(4,j,xIndex) - matrix.getQuick(0,j,xIndex) - matrix.getQuick(2,j,xIndex);
                Rg = matrix.getQuick(5,j,xIndex) - matrix.getQuick(1,j,xIndex) - matrix.getQuick(3,j,xIndex);
                
                dSfD = poissDist.nextInt(mu * matrix.getQuick(0,j,xIndex) * dt);
                dIfD = poissDist.nextInt(mu * matrix.getQuick(2,j,xIndex) * dt);
                dRfD = poissDist.nextInt(mu * Rf * dt);
                
                dSgD = poissDist.nextInt(mu * matrix.getQuick(1,j,xIndex) * dt);
                dIgD = poissDist.nextInt(mu * matrix.getQuick(3,j,xIndex) * dt);
                dRgD = poissDist.nextInt(mu * Rg * dt);

                newSf = matrix.getQuick(0,j,xIndex) + dNfB - dSfIf - dSfIg - dSfD;
                newSg = matrix.getQuick(1,j,xIndex) + dNgB - dSgIg - dSgIf - dSgD;
                newIf = matrix.getQuick(2,j,xIndex) + dSfIf + dSfIg - dIfRf - dIfD;
                if (newIf < 0.0) {
                    newIf = 0.0;
                }
                newIg = matrix.getQuick(3,j,xIndex) + dSgIg + dSgIf - dIgRg - dIgD;
                if (newIg < 0.0) {
                    newIg = 0.0;
                }
                newRf = Rf + dIfRf - dRfD;
                newRg = Rg + dIgRg - dRgD;
                newNf = newSf + newIf + newRf;
                newNg = newSg + newIg + newRg;
                
                matrix.setQuick(0,j,nextIndex,newSf);
                matrix.setQuick(1,j,nextIndex,newSg);
                matrix.setQuick(2,j,nextIndex,newIf);
                matrix.setQuick(3,j,nextIndex,newIg);
                matrix.setQuick(4,j,nextIndex,newNf);
                matrix.setQuick(5,j,nextIndex,newNg);
                matrix.setQuick(6,j,nextIndex,dSfIf);
                matrix.setQuick(7,j,nextIndex,dSfIg);
                matrix.setQuick(8,j,nextIndex,dSgIg);
                matrix.setQuick(9,j,nextIndex,dSgIf);
                matrix.setQuick(10,j,nextIndex,dIfRf);
                matrix.setQuick(11,j,nextIndex,dIgRg);
                matrix.setQuick(12,j,nextIndex,dNfB);
                matrix.setQuick(13,j,nextIndex,dNgB);
                matrix.setQuick(14,j,nextIndex,dSfD);
                matrix.setQuick(15,j,nextIndex,dIfD);
                matrix.setQuick(16,j,nextIndex,dRfD);
                matrix.setQuick(17,j,nextIndex,dSgD);
                matrix.setQuick(18,j,nextIndex,dIgD);
                matrix.setQuick(19,j,nextIndex,dRgD);
            }	
        }
        return matrix;
    }
    
    @Override
    public DoubleMatrix3D updateStatesDeterministic(int jParticles, DoubleMatrix3D matrix, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
    {
        
        //Spatial SIR with focal and global population
        //Order of state variables is 1) Sf 2) Sg 3) If 4) Ig 5) Nf 6) Ng
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            double dt = xDtTimes.getQuick(n);
            int nextIndex = xIndex + 1;
            
            double timeNow = xSubTimes.get(n);
            
            //Compute betaSeas for focal pop
            seasNowFocal = ((timeNow/365.25) + deltaFocal) - Math.floor(timeNow/365.25);
            seasHeightFocal = Math.cos(2*Math.PI * seasNowFocal);
            betaSeasFocal = beta * (1 + alphaFocal*seasHeightFocal);
            
            //Compute betaSeas for global pop
            seasNowGlobal = ((timeNow/365.25) + deltaGlobal) - Math.floor(timeNow/365.25);
            seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
            betaSeasGlobal = beta * (1 + alphaGlobal*seasHeightGlobal);
            
            for (int j = 0; j < jParticles; ++j) {
                
                //double nScaler = 1.0;

                //INFECTIONS IN FOCAL POP FROM FOCAL POP
                dSfIf = ((1-M) * betaSeasFocal * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));
                
                //INFECTIONS IN FOCAL POP FROM GLOBAL POP
                dSfIg = (betaSeasFocal * M * matrix.getQuick(0,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));

                //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
                dSgIg = ((1-M) * betaSeasGlobal * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(5,j,xIndex));

                //INFECTIONS IN GLOBAL POP FROM FOCAL POP
                dSgIf = (betaSeasGlobal * M * matrix.getQuick(1,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(5,j,xIndex));
                
                //RECOVERIES
                dIfRf = (nu * matrix.getQuick(2,j,xIndex) * dt);
                dIgRg = (nu * matrix.getQuick(3,j,xIndex) * dt);
                
                //BIRTHS
                dNfB = (mu * matrix.getQuick(4,j,xIndex) * dt);
                dNgB = (mu * matrix.getQuick(5,j,xIndex) * dt);
                
                //DEATHS IN FOCAL POP
                Rf = matrix.getQuick(4,j,xIndex) - matrix.getQuick(0,j,xIndex) - matrix.getQuick(2,j,xIndex);
                Rg = matrix.getQuick(5,j,xIndex) - matrix.getQuick(1,j,xIndex) - matrix.getQuick(3,j,xIndex);
                
                dSfD = (mu * matrix.getQuick(0,j,xIndex) * dt);
                dIfD = (mu * matrix.getQuick(2,j,xIndex) * dt);
                dRfD = (mu * Rf * dt);
                
                dSgD = (mu * matrix.getQuick(1,j,xIndex) * dt);
                dIgD = (mu * matrix.getQuick(3,j,xIndex) * dt);
                dRgD = (mu * Rg * dt);

                newSf = matrix.getQuick(0,j,xIndex) + dNfB - dSfIf - dSfIg - dSfD;
                newSg = matrix.getQuick(1,j,xIndex) + dNgB - dSgIg - dSgIf - dSgD;
                newIf = matrix.getQuick(2,j,xIndex) + dSfIf + dSfIg - dIfRf - dIfD;
                if (newIf < 0.0) {
                    newIf = 0.0;
                }
                newIg = matrix.getQuick(3,j,xIndex) + dSgIg + dSgIf - dIgRg - dIgD;
                if (newIg < 0.0) {
                    newIg = 0.0;
                }
                newRf = Rf + dIfRf - dRfD;
                newRg = Rg + dIgRg - dRgD;
                newNf = newSf + newIf + newRf;
                newNg = newSg + newIg + newRg;
                
                matrix.setQuick(0,j,nextIndex,newSf);
                matrix.setQuick(1,j,nextIndex,newSg);
                matrix.setQuick(2,j,nextIndex,newIf);
                matrix.setQuick(3,j,nextIndex,newIg);
                matrix.setQuick(4,j,nextIndex,newNf);
                matrix.setQuick(5,j,nextIndex,newNg);
                matrix.setQuick(6,j,nextIndex,dSfIf);
                matrix.setQuick(7,j,nextIndex,dSfIg);
                matrix.setQuick(8,j,nextIndex,dSgIg);
                matrix.setQuick(9,j,nextIndex,dSgIf);
                matrix.setQuick(10,j,nextIndex,dIfRf);
                matrix.setQuick(11,j,nextIndex,dIgRg);
                matrix.setQuick(12,j,nextIndex,dNfB);
                matrix.setQuick(13,j,nextIndex,dNgB);
                matrix.setQuick(14,j,nextIndex,dSfD);
                matrix.setQuick(15,j,nextIndex,dIfD);
                matrix.setQuick(16,j,nextIndex,dRfD);
                matrix.setQuick(17,j,nextIndex,dSgD);
                matrix.setQuick(18,j,nextIndex,dIgD);
                matrix.setQuick(19,j,nextIndex,dRgD);
            }	
        }
        return matrix;
    }
    
    @Override
    public DoubleMatrix3D updateStates(int jParticles, DoubleMatrix3D matrix, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
    {
        
        //Spatial SIR with focal and global population
        //Order of state variables is 1) Sf 2) Sg 3) If 4) Ig 5) Nf 6) Ng
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            double dt = xDtTimes.getQuick(n);
            int nextIndex = xIndex + 1;
            //if (dt < 0) {
                //nextIndex = xIndex - 1;
            //}
            
            double timeNow = xSubTimes.get(n);
            
            //Compute betaSeas for focal pop
            double seasNow = ((timeNow/365.25) + deltaFocal) - Math.floor(timeNow/365.25);
            double seasHeight = Math.cos(2*Math.PI * seasNow);
            //if (seasHeight < 0) {
                //seasHeight = 0;
            //}
            betaSeasFocal = beta * (1 + alphaFocal*seasHeight);
            
            //Compute betaSeas for global pop
            seasNow = ((timeNow/365.25) + deltaGlobal) - Math.floor(timeNow/365.25);
            seasHeight = Math.cos(2*Math.PI * seasNow);
            //if (seasHeight < 0) {
                //seasHeight = 0;
            //}
            betaSeasGlobal = beta * (1 + alphaGlobal*seasHeight);
            
            for (int j = 0; j < jParticles; ++j) {
                
                //double nScaler = 1.0;

                //INFECTIONS IN FOCAL POP FROM FOCAL POP
                stndNormRand = stndNorm.nextDouble();
                eta = stndNormRand / Math.sqrt(dt); //used for this particle and time only
                fTerm = fNoise * (1-M) * betaSeasFocal * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * eta * dt / matrix.getQuick(4, j, xIndex);
                //dTerm = stndNorm.nextDouble() * Math.sqrt((1-M) * betaSeasFocal * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(4,j,xIndex) * nScaler);
                dSfIf = ((1-M) * betaSeasFocal * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(4,j,xIndex)) + fTerm;
                if (dSfIf < 0.0) {
                    dSfIf = 0.0;
                }
                
                //INFECTIONS IN FOCAL POP FROM GLOBAL POP
                //dTerm = stndNorm.nextDouble() * Math.sqrt(betaSeasFocal * M * matrix.getQuick(0,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(4,j,xIndex) * nScaler);
                dSfIg = (betaSeasFocal * M * matrix.getQuick(0,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));
                if (dSfIg < 0.0) {
                    dSfIg = 0.0;
                }

                //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
                stndNormRand = stndNorm.nextDouble();
                eta = stndNormRand / Math.sqrt(dt); //used for this particle and time only
                fTerm = fNoise * (1-M) * betaSeasGlobal * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * eta * dt / matrix.getQuick(5, j, xIndex);
                //dTerm = stndNorm.nextDouble() * Math.sqrt((1-M) * betaSeasGlobal * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(5,j,xIndex) * nScaler);
                dSgIg = ((1-M) * betaSeasGlobal * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(5,j,xIndex)) + fTerm;
                if (dSgIg < 0.0) {
                    dSgIg = 0.0;
                }

                //INFECTIONS IN GLOBAL POP FROM FOCAL POP
                //dTerm = stndNorm.nextDouble() * Math.sqrt(betaSeasGlobal * M * matrix.getQuick(1,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(5,j,xIndex) * nScaler);
                dSgIf = (betaSeasGlobal * M * matrix.getQuick(1,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(5,j,xIndex));
                if (dSgIf < 0.0) {
                    dSgIf = 0.0;
                }
                
                //RECOVERIES
                //dTerm = stndNorm.nextDouble() * Math.sqrt(nu * matrix.getQuick(2,j,xIndex) * dt * nScaler);
                dIfRf = (nu * matrix.getQuick(2,j,xIndex) * dt);
                if (dIfRf < 0.0) {
                    dIfRf = 0.0;
                }
                //dTerm = stndNorm.nextDouble() * Math.sqrt(nu * matrix.getQuick(3,j,xIndex) * dt * nScaler);
                dIgRg = (nu * matrix.getQuick(3,j,xIndex) * dt);
                if (dIgRg < 0.0) {
                    dIgRg = 0.0;
                }
                
                //BIRTHS
                //dTerm = stndNorm.nextDouble() * Math.sqrt(mu * matrix.getQuick(4,j,xIndex) * dt * nScaler);
                dNfB = (mu * matrix.getQuick(4,j,xIndex) * dt);
                if (dNfB < 0.0) {
                    dNfB = 0.0;
                }
                //dTerm = stndNorm.nextDouble() * Math.sqrt(mu * matrix.getQuick(5,j,xIndex) * dt * nScaler);
                dNgB = (mu * matrix.getQuick(5,j,xIndex) * dt);
                if (dNgB < 0.0) {
                    dNgB = 0.0;
                }
                
                //DEATHS IN FOCAL POP
                Rf = matrix.getQuick(4,j,xIndex) - matrix.getQuick(0,j,xIndex) - matrix.getQuick(2,j,xIndex);
                Rg = matrix.getQuick(5,j,xIndex) - matrix.getQuick(1,j,xIndex) - matrix.getQuick(3,j,xIndex);
                
                //dTerm = stndNorm.nextDouble() * Math.sqrt(mu * matrix.getQuick(0,j,xIndex) * dt * nScaler);
                dSfD = (mu * matrix.getQuick(0,j,xIndex) * dt);
                if (dSfD < 0.0) {
                    dSfD = 0.0;
                }
                //dTerm = stndNorm.nextDouble() * Math.sqrt(mu * matrix.getQuick(2,j,xIndex) * dt * nScaler);
                dIfD = (mu * matrix.getQuick(2,j,xIndex) * dt);
                if (dIfD < 0.0) {
                    dIfD = 0.0;
                }
                //dTerm = stndNorm.nextDouble() * Math.sqrt(mu * Rf * dt * nScaler);
                dRfD = (mu * Rf * dt);
                if (dRfD < 0.0) {
                    dRfD = 0.0;
                }
                
                //dTerm = stndNorm.nextDouble() * Math.sqrt(mu * matrix.getQuick(1,j,xIndex) * dt * nScaler);
                dSgD = (mu * matrix.getQuick(1,j,xIndex) * dt);
                if (dSgD < 0.0) {
                    dSgD = 0.0;
                }
                //dTerm = stndNorm.nextDouble() * Math.sqrt(mu * matrix.getQuick(3,j,xIndex) * dt * nScaler);
                dIgD = (mu * matrix.getQuick(3,j,xIndex) * dt);
                if (dIgD < 0.0) {
                    dIgD = 0.0;
                }
                //dTerm = stndNorm.nextDouble() * Math.sqrt(mu * Rg * dt * nScaler);
                dRgD = (mu * Rg * dt);
                if (dRgD < 0.0) {
                    dRgD = 0.0;
                }

                newSf = matrix.getQuick(0,j,xIndex) + dNfB - dSfIf - dSfIg - dSfD;
                newSg = matrix.getQuick(1,j,xIndex) + dNgB - dSgIg - dSgIf - dSgD;
                newIf = matrix.getQuick(2,j,xIndex) + dSfIf + dSfIg - dIfRf - dIfD;
                if (newIf < 0.0) {
                    newIf = 0.0;
                }
                newIg = matrix.getQuick(3,j,xIndex) + dSgIg + dSgIf - dIgRg - dIgD;
                if (newIg < 0.0) {
                    newIg = 0.0;
                }
                newRf = Rf + dIfRf - dRfD;
                newRg = Rg + dIgRg - dRgD;
                newNf = newSf + newIf + newRf;
                newNg = newSg + newIg + newRg;
                
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
    
    public double computeTransProb(DoubleMatrix1D xTimeNow, DoubleMatrix1D xTimePlusOne, double dt, double timeNow) 
    {
        double transProb = 1.0;

        //beta = params.getQuick(0);
	//alphaFocal = params.getQuick(1);
        //alphaGlobal = params.getQuick(2);
        //deltaFocal = params.getQuick(3) * 365.25;
        //deltaGlobal = params.getQuick(4) * 365.25;
        //M = params.getQuick(5);
        //mu = params.getQuick(8);
        //nu = params.getQuick(9);
        
        //For f pop
        seasNowFocal = ((timeNow/365.25) + deltaFocal) - Math.floor(timeNow/365.25);
        seasHeightFocal = Math.cos(2*Math.PI * seasNowFocal);
        betaSeasFocal = beta * (1 + alphaFocal*seasHeightFocal);
        dNfB = mu * xTimeNow.getQuick(4) * dt;
        dSfIf = ((1-M) * betaSeasFocal * xTimeNow.getQuick(0) * xTimeNow.getQuick(2) * dt / xTimeNow.getQuick(4));
        dSfIg = betaSeasFocal * M * xTimeNow.getQuick(0) * xTimeNow.getQuick(3) * dt / xTimeNow.getQuick(4);
        dSfD = mu * xTimeNow.getQuick(0) * dt;
        dIfRf = nu * xTimeNow.getQuick(2) * dt;
        dIfD = mu * xTimeNow.getQuick(2) * dt;
        Rf = xTimeNow.getQuick(4) - xTimeNow.getQuick(2) - xTimeNow.getQuick(0);
        dRfD = mu * Rf * dt;
        //double nextR = xTimePlusOne.getQuick(4) - xTimePlusOne.getQuick(2) - xTimePlusOne.getQuick(0);
        
        poissDist.setMean(dNfB);
        transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(12));
        poissDist.setMean(dSfIf);
        transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(6));
        poissDist.setMean(dSfIg);
        transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(7));
        poissDist.setMean(dSfD);
        transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(14));
        poissDist.setMean(dIfRf);
        transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(10));
        poissDist.setMean(dIfD);
        transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(15));
        poissDist.setMean(dRfD);
        transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(16));
        
        //For g pop
        seasNowGlobal = ((timeNow/365.25) + deltaGlobal) - Math.floor(timeNow/365.25);
        seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
        betaSeasGlobal = beta * (1 + alphaGlobal*seasHeightGlobal);
        dNgB = mu * xTimeNow.getQuick(5) * dt;
        dSgIg = ((1-M) * betaSeasGlobal * xTimeNow.getQuick(1) * xTimeNow.getQuick(3) * dt / xTimeNow.getQuick(5));
        dSgIf = betaSeasGlobal * M * xTimeNow.getQuick(1) * xTimeNow.getQuick(2) * dt / xTimeNow.getQuick(5);
        dSgD = mu * xTimeNow.getQuick(1) * dt;
        dIgRg = nu * xTimeNow.getQuick(3) * dt;
        dIgD = mu * xTimeNow.getQuick(3) * dt;
        Rg = xTimeNow.getQuick(5) - xTimeNow.getQuick(3) - xTimeNow.getQuick(1);
        dRgD = mu * Rg * dt;
        //nextR = xTimePlusOne.getQuick(5) - xTimePlusOne.getQuick(3) - xTimePlusOne.getQuick(1);
        
        poissDist.setMean(dNgB);
        transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(13));
        poissDist.setMean(dSgIg);
        transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(8));
        poissDist.setMean(dSgIf);
        transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(9));
        poissDist.setMean(dSgD);
        transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(17));
        poissDist.setMean(dIgRg);
        transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(11));
        poissDist.setMean(dIgD);
        transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(18));
        poissDist.setMean(dRgD);
        transProb *= poissDist.pdf((int) xTimePlusOne.getQuick(19));

        if (transProb <= 0.0) {
            transProb = Double.MIN_VALUE;
            //System.out.println("Pop state trans prob is < 0.0");
        }
        if (Double.isNaN(transProb)) {
            transProb = Double.MIN_VALUE;
            System.out.println("Pop state trans prob is NaN");
        }
        return transProb;
    }
    
    
    public double computeTransProbOld(DoubleMatrix1D xTimeNow, DoubleMatrix1D xTimePlusOne, double dt, double timeNow) 
    {
        double transProb = 1.0;

        double beta = params.getQuick(0);
	double alphaFocal = params.getQuick(1);
        double alphaGlobal = params.getQuick(2);
        double deltaFocal = params.getQuick(3) * 365.25;
        double deltaGlobal = params.getQuick(4) * 365.25;
        double M = params.getQuick(5);
        double mu = params.getQuick(8);
        double nu = params.getQuick(9);
        
        //For f pop
        double seasNow = ((timeNow + deltaFocal)/365.25) - Math.floor(timeNow/365.25);
        double seasHeight = Math.cos(2*Math.PI * seasNow);
        double betaSeas = beta * (1 + alphaFocal*seasHeight);
        double dNB = mu * xTimeNow.getQuick(4) * dt;
        double dSfIf = ((1-M) * betaSeas * xTimeNow.getQuick(0) * xTimeNow.getQuick(2) * dt / xTimeNow.getQuick(4));
        double dSfIg = betaSeas * M * xTimeNow.getQuick(0) * xTimeNow.getQuick(3) * dt / xTimeNow.getQuick(4);
        double dSD = mu * xTimeNow.getQuick(0) * dt;
        double dIR = nu * xTimeNow.getQuick(2) * dt;
        double dID = mu * xTimeNow.getQuick(2) * dt;
        double currR = xTimeNow.getQuick(4) - xTimeNow.getQuick(2) - xTimeNow.getQuick(0);
        double dRD = mu * currR * dt;
        double nextR = xTimePlusOne.getQuick(4) - xTimePlusOne.getQuick(2) - xTimePlusOne.getQuick(0);
        
        DoubleMatrix1D means = factory1D.make(6);
        means.setQuick(0, xTimeNow.getQuick(0) + dNB - dSfIf - dSfIg - dSD); //this is Sf
        means.setQuick(1, xTimeNow.getQuick(2) + dSfIf + dSfIg - dIR - dID); //this is If
        means.setQuick(2, currR + dIR - dRD);
        
        DoubleMatrix1D x = factory1D.make(6);
        x.setQuick(0, xTimePlusOne.getQuick(0));
        x.setQuick(1, xTimePlusOne.getQuick(2));
        x.setQuick(2, nextR);
        
        DoubleMatrix2D covMatrix = factory2D.make(6,6);
        double nScaler = 1.0;
        
        //Don't actually need to take abs values here
        
        covMatrix.setQuick(0,0, nScaler * (Math.abs(dNB) + Math.abs(dSfIf) + Math.abs(dSfIg) + Math.abs(dSD)));
        covMatrix.setQuick(1,1, nScaler * (Math.abs(dSfIf) + Math.abs(dSfIg) + Math.abs(dIR) + Math.abs(dID)));
        covMatrix.setQuick(2,2, nScaler * (Math.abs(dIR) + Math.abs(dRD)));
        covMatrix.setQuick(0,1, -1.0 * nScaler * (Math.abs(dSfIf) + Math.abs(dSfIg)));
        covMatrix.setQuick(1,0, -1.0 * nScaler * (Math.abs(dSfIf) + Math.abs(dSfIg)));
        covMatrix.setQuick(1,2, -1.0 * nScaler * (Math.abs(dIR)));
        covMatrix.setQuick(2,1, -1.0 * nScaler * (Math.abs(dIR)));
        
        //Compute prob for Sf (sum of the birth, death and the 2 infection processes
//        double expectedX = xTimeNow.getQuick(0) + dNB - dSfIf - dSfIg - dSD;
//        double devX = expectedX - xTimePlusOne.getQuick(0);
//        double varX = dNB + dSfIf + dSfIg + dSD;
//        double exponent = -Math.pow(devX,2.0) / (2.0 * varX);
//        transProb *= Math.exp(exponent)/Math.sqrt(2 * Math.PI * varX);
//        
//        //Compute prob for If (sum of recovery, death and 2 infection processes)
//        expectedX = xTimeNow.getQuick(2) + dSfIf + dSfIg - dIR - dID;
//        devX = expectedX - xTimePlusOne.getQuick(2);
//        varX = dSfIf + dSfIg + dIR + dID;
//        exponent = -Math.pow(devX,2.0) / (2.0 * varX);
//        transProb *= Math.exp(exponent)/Math.sqrt(2 * Math.PI * varX);
//        
//        //Compute prob for Rf (sum of recovery and death processes)
//        expectedX = currR + dIR - dRD;
//        devX = expectedX - nextR;
//        varX = dIR + dRD;
//        exponent = -Math.pow(devX,2.0) / (2.0 * varX);
//        transProb *= Math.exp(exponent)/Math.sqrt(2 * Math.PI * varX);
        
        //For g pop
        seasNow = ((timeNow + deltaGlobal)/365.25) - Math.floor(timeNow/365.25);
        seasHeight = Math.cos(2*Math.PI * seasNow);
        betaSeas = beta * (1 + alphaGlobal*seasHeight);
        dNB = mu * xTimeNow.getQuick(5) * dt;
        dSfIf = ((1-M) * betaSeas * xTimeNow.getQuick(1) * xTimeNow.getQuick(3) * dt / xTimeNow.getQuick(5));
        dSfIg = betaSeas * M * xTimeNow.getQuick(1) * xTimeNow.getQuick(2) * dt / xTimeNow.getQuick(5);
        dSD = mu * xTimeNow.getQuick(1) * dt;
        dIR = nu * xTimeNow.getQuick(3) * dt;
        dID = mu * xTimeNow.getQuick(3) * dt;
        currR = xTimeNow.getQuick(5) - xTimeNow.getQuick(3) - xTimeNow.getQuick(1);
        dRD = mu * currR * dt;
        nextR = xTimePlusOne.getQuick(5) - xTimePlusOne.getQuick(3) - xTimePlusOne.getQuick(1);
        
        means.setQuick(3, xTimeNow.getQuick(1) + dNB - dSfIf - dSfIg - dSD); //this is Sf
        means.setQuick(4, xTimeNow.getQuick(3) + dSfIf + dSfIg - dIR - dID); //this is If
        means.setQuick(5, currR + dIR - dRD);
        
        x.setQuick(3, xTimePlusOne.getQuick(1));
        x.setQuick(4, xTimePlusOne.getQuick(3));
        x.setQuick(5, nextR);
        
        covMatrix.setQuick(3,3, nScaler * (Math.abs(dNB) + Math.abs(dSfIf) + Math.abs(dSfIg) + Math.abs(dSD)));
        covMatrix.setQuick(4,4, nScaler * (Math.abs(dSfIf) + Math.abs(dSfIg) + Math.abs(dIR) + Math.abs(dID)));
        covMatrix.setQuick(5,5, nScaler * (Math.abs(dIR) + Math.abs(dRD)));
        covMatrix.setQuick(3,4, -1.0 * nScaler * (Math.abs(dSfIf) + Math.abs(dSfIg)));
        covMatrix.setQuick(4,3, -1.0 * nScaler * (Math.abs(dSfIf) + Math.abs(dSfIg)));
        covMatrix.setQuick(4,5, -1.0 * nScaler * (Math.abs(dIR)));
        covMatrix.setQuick(5,4, -1.0 * nScaler * (Math.abs(dIR)));
        
        //Compute prob for Sg (sum of the birth, death and the 2 infection processes
//        expectedX = xTimeNow.getQuick(1) + dNB - dSfIf - dSfIg - dSD;
//        devX = expectedX - xTimePlusOne.getQuick(1);
//        varX = dNB + dSfIf + dSfIg + dSD;
//        exponent = -Math.pow(devX,2.0) / (2.0 * varX);
//        transProb *= Math.exp(exponent)/Math.sqrt(2 * Math.PI * varX);
//        
//        //Compute prob for Ig (sum of recovery, death and 2 infection processes)
//        expectedX = xTimeNow.getQuick(3) + dSfIf + dSfIg - dIR - dID;
//        devX = expectedX - xTimePlusOne.getQuick(3);
//        varX = dSfIf + dSfIg + dIR + dID;
//        exponent = -Math.pow(devX,2.0) / (2.0 * varX);
//        transProb *= Math.exp(exponent)/Math.sqrt(2 * Math.PI * varX);
//        
//        //Compute prob for Rg (sum of recovery and death processes)
//        expectedX = currR + dIR - dRD;
//        devX = expectedX - nextR;
//        varX = dIR + dRD;
//        exponent = -Math.pow(devX,2.0) / (2.0 * varX);
//        transProb *= Math.exp(exponent)/Math.sqrt(2 * Math.PI * varX);
        
        
        transProb = MatrixUtils.getMultiVarNormProb(x, means, covMatrix);
        //System.out.println("transProb = " + transProb);
        
        if (transProb <= 0.0) {
            transProb = Double.MIN_VALUE;
            //System.out.println("Pop state trans prob is < 0.0");
        }
        if (Double.isNaN(transProb)) {
            transProb = Double.MIN_VALUE;
            System.out.println("Pop state trans prob is NaN");
        }
        return transProb;
    }
    
    
    public double computeTransProbOlder(DoubleMatrix1D xTimeNow, DoubleMatrix1D xTimePlusOne, double dt, double timeNow) 
    {
        double transProb = 1.0;

        double beta = params.getQuick(0);
	double alphaFocal = params.getQuick(1);
        double alphaGlobal = params.getQuick(2);
        double deltaFocal = params.getQuick(3) * 365.25;
        double deltaGlobal = params.getQuick(4) * 365.25;
        double M = params.getQuick(5);
        double mu = params.getQuick(8);
        double nu = params.getQuick(9);
        double fNoise = params.getQuick(10);
        //double standDev = fNoise * Math.sqrt(dt);
        
        //For f pop
        double seasNow = ((timeNow + deltaFocal)/365.25) - Math.floor(timeNow/365.25);
        double seasHeight = Math.cos(2*Math.PI * seasNow);
        double betaSeas = beta * (1 + alphaFocal*seasHeight);
        double dNfB = mu * xTimeNow.getQuick(4) * dt;
        double dSfIf = ((1-M) * betaSeas * xTimeNow.getQuick(0) * xTimeNow.getQuick(2) * dt / xTimeNow.getQuick(4));
        double dSfIg = betaSeas * M * xTimeNow.getQuick(0) * xTimeNow.getQuick(3) * dt / xTimeNow.getQuick(4);
        double dSfD = mu * xTimeNow.getQuick(0) * dt;
        //double dIfRf = nu * xTimeNow.getQuick(2) * dt;
        //double dIfD = mu * xTimeNow.getQuick(2) * dt;
        double expectedSfPlusOne = xTimeNow.getQuick(0) + dNfB - dSfIf - dSfIg - dSfD;
        double deviationFromExpectation = expectedSfPlusOne - xTimePlusOne.getQuick(0);
        //double scaledDev = deviationFromExpectation / ((fNoise * (1-M) * betaSeas * xTimeNow.getQuick(0) * xTimeNow.getQuick(2) * Math.sqrt(dt))/xTimeNow.getQuick(4));
        double standDev = ((fNoise * (1-M) * betaSeas * xTimeNow.getQuick(0) * xTimeNow.getQuick(2) * Math.sqrt(dt))/xTimeNow.getQuick(4));
        //double standDev = 1.0;
        double exponent = -Math.pow(deviationFromExpectation,2.0) / (2.0 * Math.pow(standDev,2.0));
        transProb *= Math.exp(exponent)/(standDev * Math.sqrt(2*Math.PI));
        if (transProb <= 0.0) {
            transProb = Double.MIN_VALUE;
        }
//        double expectedIfPlusOne = xTimeNow.getQuick(2) + dSfIf + dSfIg - dIfRf - dIfD;
//        deviationFromExpectation = expectedIfPlusOne - xTimePlusOne.getQuick(2);
//        scaledDev = deviationFromExpectation / ((fNoise * (1-M) * betaSeas * xTimeNow.getQuick(0) * xTimeNow.getQuick(2) * Math.sqrt(dt))/xTimeNow.getQuick(4));
//        exponent = -Math.pow(scaledDev,2.0) / (2.0 * Math.pow(standDev,2.0));
//        transProb *= Math.exp(exponent)/(standDev * Math.sqrt(2*Math.PI));
//        if (transProb <= 0.0) {
//            transProb = Double.MIN_VALUE;
//        }
        
        
        //For g pop
        seasNow = ((timeNow + deltaGlobal)/365.25) - Math.floor(timeNow/365.25);
        seasHeight = Math.cos(2*Math.PI * seasNow);
        betaSeas = beta * (1 + alphaGlobal*seasHeight);
        double dNgB = mu * xTimeNow.getQuick(5) * dt;
        double dSgIg = ((1-M) * betaSeas * xTimeNow.getQuick(1) * xTimeNow.getQuick(3) * dt / xTimeNow.getQuick(5));
        double dSgIf = betaSeas * M * xTimeNow.getQuick(1) * xTimeNow.getQuick(2) * dt / xTimeNow.getQuick(5);
        double dSgD = mu * xTimeNow.getQuick(1) * dt;
        //double dIgRg = nu * xTimeNow.getQuick(3) * dt;
        //double dIgD = mu * xTimeNow.getQuick(3) * dt;
        double expectedSgPlusOne = xTimeNow.getQuick(1) + dNgB - dSgIg - dSgIf - dSgD;
        deviationFromExpectation = expectedSgPlusOne - xTimePlusOne.getQuick(1);
        //scaledDev = deviationFromExpectation / ((fNoise * (1-M) * betaSeas * xTimeNow.getQuick(1) * xTimeNow.getQuick(3) * Math.sqrt(dt))/xTimeNow.getQuick(5));
        standDev = ((fNoise * (1-M) * betaSeas * xTimeNow.getQuick(1) * xTimeNow.getQuick(3) * Math.sqrt(dt))/xTimeNow.getQuick(5));
        exponent = -Math.pow(deviationFromExpectation,2.0) / (2.0 * Math.pow(standDev,2.0));
        transProb *= Math.exp(exponent)/(standDev * Math.sqrt(2*Math.PI));
        if (transProb <= 0.0) {
            transProb = Double.MIN_VALUE;
        }
//        double expectedIgPlusOne = xTimeNow.getQuick(3) + dSgIg + dSgIf - dIgRg - dIgD;
//        deviationFromExpectation = expectedIgPlusOne - xTimePlusOne.getQuick(3);
//        scaledDev = deviationFromExpectation / ((fNoise * (1-M) * betaSeas * xTimeNow.getQuick(1) * xTimeNow.getQuick(3) * Math.sqrt(dt))/xTimeNow.getQuick(5));
//        exponent = -Math.pow(scaledDev,2.0) / (2.0 * Math.pow(standDev,2.0));
//        transProb *= Math.exp(exponent)/(standDev * Math.sqrt(2*Math.PI));

        
        //System.out.println("Trans prob = " + transProb);
        return transProb;
    }
    
       
}
