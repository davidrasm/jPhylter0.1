import cern.colt.list.*;
import java.util.ArrayList;
import cern.jet.random.*;
import cern.colt.matrix.*;
//import cern.jet.random.AbstractDistribution;
//import cern.colt.matrix.DoubleMatrix1D;
//import cern.colt.matrix.DoubleMatrix3D;

/**
 * Specifies the epidemiological model
 * Two-population model with vector in focal pop
 * @author David
 */
public class EpiModel2PopVectored extends EpiModel {

    Poisson poissDist;
    
    double seasNowFocal; double seasNowGlobal; double seasHeightFocal; double seasHeightGlobal; double betaSeasFF; double betaSeasGF; double betaSeasGG; double betaSeasFG;
    double betaFF; double betaGF; double betaGG; double betaFG; double alphaFocal; double alphaGlobal; double deltaFocal; double deltaGlobal;
    double NFocal; double NGlobal; double mu; double nu; //double fNoise;
    
    double dSfIf; double dSfIg; double dSgIg; double dSgIf; double dIfRf; double dIgRg;
    double dSfD; double dIfD; double dRfD; double dSgD; double dIgD; double dRgD; double dNfB; double dNgB; double Rf; double Rg;
    double newSf; double newSg; double newIf; double newIg; double newRf;  double newRg; double newNf; double newNg;
    
    double M; double seasM; double muVec;
    
    double initSFocal;
    
    @Override
    public void setInitParams()
    {
        //Set params to initial values and determine their order in params list
        //betaFF = 3.0 * ((1/(70*365.25)) + (1.0/7.0)); params.add(betaFF); estParams.add(0);
        betaFF = 0.22; params.add(betaFF); estParams.add(0);
        betaGF = 0.0; params.add(betaGF); fixedParams.add(1);
        betaGG = 0.16; params.add(betaGG); fixedParams.add(2);
        betaFG = betaGF; params.add(betaFG); fixedParams.add(3);
        
        alphaFocal = 0.08; params.add(alphaFocal); estParams.add(4);
        alphaGlobal = 0.0; params.add(alphaGlobal); fixedParams.add(5);
        deltaFocal = 6.0/12.0; params.add(deltaFocal); estParams.add(6);
        deltaGlobal = 6.0/12.0; params.add(deltaGlobal); fixedParams.add(7);
        NFocal = 10.0e06; params.add(NFocal); fixedParams.add(8);
        NGlobal = 25.0e06; params.add(NGlobal); fixedParams.add(9);
        mu = 1/(60*365.25); params.add(mu); fixedParams.add(10);
        nu = 1.0/7.0; params.add(nu); fixedParams.add(11);
        
        //Additional vector params;
        M = 1.0; params.add(M); estParams.add(12);
        muVec = 1.0 / 7.0; params.add(muVec); fixedParams.add(13);
        
        //Fix initSFocal
        initSFocal = 0.40; params.add(initSFocal); estParams.add(14); //estimated as the fraction, not the absolute number
                        
    }
    
    @Override
    public void setRandomGenerator() 
    {
      
        cern.jet.random.engine.RandomEngine engine = new cern.jet.random.engine.MersenneTwister(new java.util.Date());
        poissDist = new Poisson(0.0, engine);
        
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
        betaFF = params.getQuick(0);
        //betaGF = params.getQuick(1);
        //betaGG = params.getQuick(2);
        //betaFG = params.getQuick(1); params.setQuick(3, betaFG);
        alphaFocal = params.getQuick(4);
        deltaFocal = params.getQuick(6);
        
        M = params.getQuick(12);
        initSFocal = params.getQuick(14);
        
    }
    
    public DoubleArrayList getInitialConditions() 
    {
        
        //Order init conditions 1) Sf 2) Sg 3) If 4) Ig 5) Nf 6) Ng
        
        DoubleArrayList currInits = new DoubleArrayList();
       
        //Compute equilibrium conditons without migration
        //double R0Focal = betaFF / (mu + nu); //assuming no migration and no vector
        double R0Focal = (betaFF*betaFF * M) / (muVec * (mu + nu));
        double R0Global = betaGG / (mu + nu); //assuming no migration
        
        double initS1 = NFocal/R0Focal; currInits.add(initS1);
        double initS2 = NGlobal/R0Global; currInits.add(initS2);
        double initI1 = mu * NFocal * (R0Focal - 1) / betaFF; currInits.add(initI1);
        double initI2 = mu * NGlobal * (R0Global - 1) / betaGG; currInits.add(initI2);
        currInits.add(NFocal);
        currInits.add(NGlobal);
        
        double dt = 0.25;
        double length = 365.25*200;
        currInits = this.fastSim(currInits, dt, length);
        
        //If estimating init S1
        initS1 = initSFocal * NFocal; currInits.set(0, initS1); //estimating Sf init
        
        //True inits for dengue sim data in DengueSimDataNoNoiseM001.tre (t=14.0*365.25)
        //double initS1 = 664246.057117843; currInits.add(initS1);
        //double initS2 = 1320132.38209387; currInits.add(initS2);
        //double initI1 = 179.547431957558; currInits.add(initI1);
        //double initI2 = 93.4958120618454; currInits.add(initI2);
        //double initN1 = 2.0e06; currInits.add(initN1);
        //double initN2 = 2.0e06; currInits.add(initN2);
        
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
        
        //Need to add trackers for If and Ig (or Iv)
        currInits.add(0.0); //cumulative incidence in focal pop
        currInits.add(0.0); //cumulative incidence in global pop
        
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
        
        //If estimating init S1
        if (params.get(14) > 1.0) {
            constraintCheckFail = true;
        }
        
        //deltaFocal has to be less than 365.25
        if (params.get(6) > (1.0)) { //this means delta is less than one year
            constraintCheckFail = true;
        }
        
        
        return constraintCheckFail;
    }
    
    //Method not fully updated for this model
    @Override
    public DoubleMatrix3D updateStatesTauLeap(int jParticles, DoubleMatrix3D matrix, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
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
            betaSeasFF = betaFF * (1 + alphaFocal*seasHeightFocal);
            betaSeasGF = betaGF * (1 + alphaFocal*seasHeightFocal);
            
            //If estimating seasonal M
            seasM = M * (1 + alphaFocal*seasHeightFocal);
            
            //Compute betaSeas for global pop
            seasNowGlobal = ((timeNow/365.25) + deltaGlobal) - Math.floor(timeNow/365.25);
            seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
            betaSeasGG = betaGG * (1 + alphaGlobal*seasHeightGlobal);
            betaSeasFG = betaFG * (1 + alphaGlobal*seasHeightGlobal);
            
            for (int j = 0; j < jParticles; ++j) {
                
                //Get quasi-equilibrium vector states
                //double vecN = NFocal * seasM;
                //double vecS = vecN * muVec / ((betaFF/NFocal) * matrix.getQuick(2,j,xIndex) + muVec);
                double vecI = seasM * NFocal * (betaFF/NFocal) * muVec * matrix.getQuick(2,j,xIndex) / ((betaFF/NFocal)*matrix.getQuick(2,j,xIndex)*muVec + muVec*muVec);
                

                //INFECTIONS IN FOCAL POP FROM FOCAL POP
                dSfIf = poissDist.nextInt(betaSeasFF * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));
                
                //INFECTIONS IN FOCAL POP FROM GLOBAL POP
                dSfIg = poissDist.nextInt(betaSeasGF * matrix.getQuick(0,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));

                //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
                dSgIg = poissDist.nextInt(betaSeasGG * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(5,j,xIndex));

                //INFECTIONS IN GLOBAL POP FROM FOCAL POP
                dSgIf = poissDist.nextInt(betaSeasFG * matrix.getQuick(1,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(5,j,xIndex));
                
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
            //betaSeasFF = betaFF * (1 + alphaFocal*seasHeightFocal);
            betaSeasGF = betaGF * (1 + alphaFocal*seasHeightFocal); 
            
            //If estimating seasonal M
            seasM = M * (1 + alphaFocal*seasHeightFocal);
            
            //Compute betaSeas for global pop
            seasNowGlobal = ((timeNow/365.25) + deltaGlobal) - Math.floor(timeNow/365.25);
            seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
            betaSeasGG = betaGG * (1 + alphaGlobal*seasHeightGlobal);
            betaSeasFG = betaFG * (1 + alphaGlobal*seasHeightGlobal);
            
            for (int j = 0; j < jParticles; ++j) {
                
                //Get quasi-equilibrium vector states
                //double vecN = NFocal * seasM;
                //double vecS = vecN * muVec / ((betaFF/NFocal) * matrix.getQuick(2,j,xIndex) + muVec);
                double vecI = seasM * NFocal * (betaFF/NFocal) * muVec * matrix.getQuick(2,j,xIndex) / ((betaFF/NFocal)*matrix.getQuick(2,j,xIndex)*muVec + muVec*muVec);

                //INFECTIONS IN FOCAL POP FROM FOCAL POP
                //dSfIf = (betaSeasFF * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));
                dSfIf = (betaFF * matrix.getQuick(0,j,xIndex) * vecI * dt / matrix.getQuick(4,j,xIndex));
                
                //INFECTIONS IN FOCAL POP FROM GLOBAL POP
                dSfIg = (betaSeasGF * matrix.getQuick(0,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));

                //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
                dSgIg = (betaSeasGG * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(5,j,xIndex));

                //INFECTIONS IN GLOBAL POP FROM FOCAL POP
                dSgIf = (betaSeasFG * matrix.getQuick(1,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(5,j,xIndex));
                
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
                
                //Track cumulative incidence
                matrix.setQuick(20,j,nextIndex, (dSfIf + dSfIg));
                matrix.setQuick(21,j,nextIndex, (dSgIg + dSgIf));
                
            }	
        }
        return matrix;
    }
    
    //Defunct method - needs to be update for this model
    /**
    public DoubleMatrix3D updateStates(int jParticles, DoubleMatrix3D matrix, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes, AbstractDistribution stndNorm) 
    {
        
        //Spatial SIR with focal and global population
        //Order of state variables is 1) Sf 2) Sg 3) If 4) Ig 5) Nf 6) Ng
      
        double eta; double dSfIf; double dSfIg; double dSgIg; double dSgIf; double dIfRf; double dIgRg;
        double dSfD; double dIfD; double dRfD; double dSgD; double dIgD; double dRgD; double dNfB; double dNgB; double Rf; double Rg;
        double newSf; double newSg; double newIf; double newIg; double newRf;  double newRg; double newNf; double newNg;
        double dTerm; //double fTerm;
		
        double beta = params.getQuick(0);
	double alphaFocal = params.getQuick(1);
        double alphaGlobal = params.getQuick(2);
        double deltaFocal = params.getQuick(3) * 365.25;
        double deltaGlobal = params.getQuick(4) * 365.25;
        double M = params.getQuick(5);
        double mu = params.getQuick(8);
        double nu = params.getQuick(9);
        //double fNoise = params.getQuick(10);
	//double stndNormRand;
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            double dt = xDtTimes.getQuick(n);
            int nextIndex = xIndex + 1;
            //if (dt < 0) {
                //nextIndex = xIndex - 1;
            //}
            
            double timeNow = xSubTimes.get(n);
            
            //Compute betaSeas for focal pop
            double seasNow = ((timeNow + deltaFocal)/365.25) - Math.floor(timeNow/365.25);
            double seasHeight = Math.cos(2*Math.PI * seasNow);
            //if (seasHeight < 0) {
                //seasHeight = 0;
            //}
            double betaSeasFocal = beta * (1 + alphaFocal*seasHeight);
            
            //Compute betaSeas for global pop
            seasNow = ((timeNow + deltaGlobal)/365.25) - Math.floor(timeNow/365.25);
            seasHeight = Math.cos(2*Math.PI * seasNow);
            //if (seasHeight < 0) {
                //seasHeight = 0;
            //}
            double betaSeasGlobal = beta * (1 + alphaGlobal*seasHeight);
            
            for (int j = 0; j < jParticles; ++j) {
                
                double nScaler = 1.0;

                //INFECTIONS IN FOCAL POP FROM FOCAL POP
                //stndNormRand = stndNorm.nextDouble();
                //eta = stndNormRand / Math.sqrt(Math.abs(dt)); //used for this particle and time only
                //fTerm = fNoise * (1-M) * betaSeasFocal * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * eta * dt / matrix.getQuick(4, j, xIndex);
                dTerm = stndNorm.nextDouble() * Math.sqrt((1-M) * betaSeasFocal * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(4,j,xIndex) * nScaler);
                dSfIf = ((1-M) * betaSeasFocal * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(4,j,xIndex)) + dTerm;
                if (dSfIf < 0.0) {
                    dSfIf = 0.0;
                }
                
                //INFECTIONS IN FOCAL POP FROM GLOBAL POP
                dTerm = stndNorm.nextDouble() * Math.sqrt(betaSeasFocal * M * matrix.getQuick(0,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(4,j,xIndex) * nScaler);
                dSfIg = (betaSeasFocal * M * matrix.getQuick(0,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(4,j,xIndex)) + dTerm;
                if (dSfIg < 0.0) {
                    dSfIg = 0.0;
                }

                //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
                //stndNormRand = stndNorm.nextDouble();
                //eta = stndNormRand / Math.sqrt(Math.abs(dt)); //used for this particle and time only
                //fTerm = fNoise * (1-M) * betaSeasGlobal * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * eta * dt / matrix.getQuick(5, j, xIndex);
                dTerm = stndNorm.nextDouble() * Math.sqrt((1-M) * betaSeasGlobal * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(5,j,xIndex) * nScaler);
                dSgIg = ((1-M) * betaSeasGlobal * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(5,j,xIndex)) + dTerm;
                if (dSgIg < 0.0) {
                    dSgIg = 0.0;
                }

                //INFECTIONS IN GLOBAL POP FROM FOCAL POP
                dTerm = stndNorm.nextDouble() * Math.sqrt(betaSeasGlobal * M * matrix.getQuick(1,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(5,j,xIndex) * nScaler);
                dSgIf = (betaSeasGlobal * M * matrix.getQuick(1,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(5,j,xIndex)) + dTerm;
                if (dSgIf < 0.0) {
                    dSgIf = 0.0;
                }
                
                //RECOVERIES
                dTerm = stndNorm.nextDouble() * Math.sqrt(nu * matrix.getQuick(2,j,xIndex) * dt * nScaler);
                dIfRf = (nu * matrix.getQuick(2,j,xIndex) * dt) + dTerm;
                if (dIfRf < 0.0) {
                    dIfRf = 0.0;
                }
                dTerm = stndNorm.nextDouble() * Math.sqrt(nu * matrix.getQuick(3,j,xIndex) * dt * nScaler);
                dIgRg = (nu * matrix.getQuick(3,j,xIndex) * dt) + dTerm;
                if (dIgRg < 0.0) {
                    dIgRg = 0.0;
                }
                
                //BIRTHS
                dTerm = stndNorm.nextDouble() * Math.sqrt(mu * matrix.getQuick(4,j,xIndex) * dt * nScaler);
                dNfB = (mu * matrix.getQuick(4,j,xIndex) * dt) + dTerm;
                if (dNfB < 0.0) {
                    dNfB = 0.0;
                }
                dTerm = stndNorm.nextDouble() * Math.sqrt(mu * matrix.getQuick(5,j,xIndex) * dt * nScaler);
                dNgB = (mu * matrix.getQuick(5,j,xIndex) * dt) + dTerm;
                if (dNgB < 0.0) {
                    dNgB = 0.0;
                }
                
                //DEATHS IN FOCAL POP
                Rf = matrix.getQuick(4,j,xIndex) - matrix.getQuick(0,j,xIndex) - matrix.getQuick(2,j,xIndex);
                Rg = matrix.getQuick(5,j,xIndex) - matrix.getQuick(1,j,xIndex) - matrix.getQuick(3,j,xIndex);
                
                dTerm = stndNorm.nextDouble() * Math.sqrt(mu * matrix.getQuick(0,j,xIndex) * dt * nScaler);
                dSfD = (mu * matrix.getQuick(0,j,xIndex) * dt) + dTerm;
                if (dSfD < 0.0) {
                    dSfD = 0.0;
                }
                dTerm = stndNorm.nextDouble() * Math.sqrt(mu * matrix.getQuick(2,j,xIndex) * dt * nScaler);
                dIfD = (mu * matrix.getQuick(2,j,xIndex) * dt) + dTerm;
                if (dIfD < 0.0) {
                    dIfD = 0.0;
                }
                dTerm = stndNorm.nextDouble() * Math.sqrt(mu * Rf * dt * nScaler);
                dRfD = (mu * Rf * dt) + dTerm;
                if (dRfD < 0.0) {
                    dRfD = 0.0;
                }
                
                dTerm = stndNorm.nextDouble() * Math.sqrt(mu * matrix.getQuick(1,j,xIndex) * dt * nScaler);
                dSgD = (mu * matrix.getQuick(1,j,xIndex) * dt) + dTerm;
                if (dSgD < 0.0) {
                    dSgD = 0.0;
                }
                dTerm = stndNorm.nextDouble() * Math.sqrt(mu * matrix.getQuick(3,j,xIndex) * dt * nScaler);
                dIgD = (mu * matrix.getQuick(3,j,xIndex) * dt) + dTerm;
                if (dIgD < 0.0) {
                    dIgD = 0.0;
                }
                dTerm = stndNorm.nextDouble() * Math.sqrt(mu * Rg * dt * nScaler);
                dRgD = (mu * Rg * dt) + dTerm;
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
    } */
    
    public DoubleArrayList fastSim(DoubleArrayList inits, double dt, double endTime)
    {
        
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
            seasNowFocal = ((timeNow/365.25) + deltaFocal) - Math.floor(timeNow/365.25);
            seasHeightFocal = Math.cos(2*Math.PI * seasNowFocal);
            betaSeasFF = betaFF * (1 + alphaFocal*seasHeightFocal);
            betaSeasGF = betaGF * (1 + alphaFocal*seasHeightFocal);
            
            //If estimating seasonal M
            seasM = M * (1 + alphaFocal*seasHeightFocal);
            
            //Compute betaSeas for global pop
            seasNowGlobal = ((timeNow/365.25) + deltaGlobal) - Math.floor(timeNow/365.25);
            seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
            betaSeasGG = betaGG * (1 + alphaGlobal*seasHeightGlobal);
            betaSeasFG = betaFG * (1 + alphaGlobal*seasHeightGlobal);            
            
            //stndNormRand = stndNorm.nextDouble();
            //eta = stndNormRand / Math.sqrt(dt); //used for this particle and time only
            
            //Get quasi-equilibrium vector states
            //double vecN = NFocal * seasM;
            //double vecS = vecN * muVec / ((betaFF/NFocal) * matrix.getQuick(2,j,xIndex) + muVec);
            double vecI = seasM * NFocal * (betaFF/NFocal) * muVec * x.getQuick(2) / ((betaFF/NFocal)*x.getQuick(2)*muVec + muVec*muVec);

            //INFECTIONS IN FOCAL POP FROM FOCAL POP
            //fTerm = fNoise * betaSeasFocal * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * eta * dt / matrix.getQuick(4, j, xIndex);
            //dSfIf = (betaSeasFF * x.getQuick(0) * x.getQuick(2) * dt / x.getQuick(4)); // + fTerm;
            dSfIf = (betaFF * x.getQuick(0) * vecI * dt / x.getQuick(4));

            //INFECTIONS IN FOCAL POP FROM GLOBAL POP
            dSfIg = betaSeasGF * x.getQuick(0) * x.getQuick(3) * dt / x.getQuick(4);
 

            //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
            //fTerm = fNoise * betaSeasGlobal * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * eta * dt / matrix.getQuick(5, j, xIndex);
            dSgIg = (betaSeasGG * x.getQuick(1) * x.getQuick(3) * dt / x.getQuick(5)); //+ fTerm;

            //INFECTIONS IN GLOBAL POP FROM FOCAL POP
            dSgIf = betaSeasFG * x.getQuick(1) * x.getQuick(2) * dt / x.getQuick(5);

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
    
    //Need to update this method for this model
    public double computeTransProb(DoubleMatrix1D xTimeNow, DoubleMatrix1D xTimePlusOne, double dt, double timeNow) 
    {
        double transProb = 1.0;
        
        //For f pop
        seasNowFocal = ((timeNow/365.25) + deltaFocal) - Math.floor(timeNow/365.25);
        seasHeightFocal = Math.cos(2*Math.PI * seasNowFocal);
        betaSeasFF = betaFF * (1 + alphaFocal*seasHeightFocal);
        betaSeasGF = betaGF * (1 + alphaFocal*seasHeightFocal);
        dNfB = mu * xTimeNow.getQuick(4) * dt;
        dSfIf = (betaSeasFF * xTimeNow.getQuick(0) * xTimeNow.getQuick(2) * dt / xTimeNow.getQuick(4));
        dSfIg = betaSeasGF * xTimeNow.getQuick(0) * xTimeNow.getQuick(3) * dt / xTimeNow.getQuick(4);
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
        betaSeasGG = betaGG * (1 + alphaGlobal*seasHeightGlobal);
        betaSeasFG = betaFG * (1 + alphaGlobal*seasHeightGlobal);
        dNgB = mu * xTimeNow.getQuick(5) * dt;
        dSgIg = (betaSeasGG * xTimeNow.getQuick(1) * xTimeNow.getQuick(3) * dt / xTimeNow.getQuick(5));
        dSgIf = betaSeasFG * xTimeNow.getQuick(1) * xTimeNow.getQuick(2) * dt / xTimeNow.getQuick(5);
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

       
}