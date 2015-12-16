import cern.colt.list.*;
import java.util.ArrayList;
import cern.jet.random.*;
import cern.colt.matrix.*;
//import cern.jet.random.AbstractDistribution;
//import cern.colt.matrix.DoubleMatrix1D;
//import cern.colt.matrix.DoubleMatrix3D;

/**
 * Specifies the epidemiological model
 * Katia's model with seasonal birth rate B
 * Corrected coalescent rates
 * This is model used in MBE paper for vectored and combined vector/spatial model
 * @author David
 */

public class EpiModelDengueVectored extends EpiModel {

    Poisson poissDist;
    
    double seasNowFocal; double seasNowGlobal; double seasHeightFocal; double seasHeightGlobal; double betaSeasFF; double betaSeasGF; double betaSeasGG; double betaSeasFG;
    double betaHV; double betaVH; double betaGF; double betaGG; double betaFG; double alphaFocal; double alphaGlobal; double deltaFocal; double deltaGlobal;
    double NFocal; double NGlobal; double mu; double nu; //double fNoise;
    
    double dSfIf; double dSfIg; double dSgIg; double dSgIf; double dIfRf; double dIgRg;
    double dSfD; double dIfD; double dRfD; double dSgD; double dIgD; double dRgD; double dNfB; double dNgB; double Rf; double Rg;
    double newSf; double newSg; double newIf; double newIg; double newRf;  double newRg; double newNf; double newNg;
    
    double M; double muVec;
    double B; double seasB; double seasM;
    
    double initSFocal;
    
    @Override
    public void setInitParams()
    {
        //Set params to initial values and determine their order in params list
        //betaFF = 0.22; params.add(betaFF); estParams.add(0);
        betaHV = 0.2770; params.add(betaHV); estParams.add(0);
        betaVH = betaHV; params.add(betaVH); fixedParams.add(1);
        betaGF = 0.0019; params.add(betaGF); estParams.add(2);
        betaGG = 0.1595; params.add(betaGG); estParams.add(3);
        betaFG = betaGF; params.add(betaFG); fixedParams.add(4);
        
        alphaFocal = 0.11; params.add(alphaFocal); estParams.add(5);
        alphaGlobal = 0.0; params.add(alphaGlobal); fixedParams.add(6);
        deltaFocal = 6.0/12.0; params.add(deltaFocal); estParams.add(7);
        deltaGlobal = 6.0/12.0; params.add(deltaGlobal); fixedParams.add(8);
        NFocal = 10.0e06; params.add(NFocal); fixedParams.add(9);
        NGlobal = 25.0e06; params.add(NGlobal); fixedParams.add(10);
        mu = 1/(60*365.25); params.add(mu); fixedParams.add(11);
        nu = 1.0/7.0; params.add(nu); fixedParams.add(12);
        
        //Additional vector params;
        M = 0.65; params.add(M); estParams.add(13);
        muVec = 1.0 / 7.0; params.add(muVec); fixedParams.add(14);
        
        //Fix initSFocal
        initSFocal = 0.4041; params.add(initSFocal); estParams.add(15); //estimated as the fraction, not the absolute number
                        
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
        betaHV = params.getQuick(0);
        betaVH = params.getQuick(0); params.setQuick(1, betaVH);
        betaGF = params.getQuick(2);
        betaGG = params.getQuick(3);
        betaFG = params.getQuick(2); params.setQuick(4, betaFG);
        alphaFocal = params.getQuick(5);
        deltaFocal = params.getQuick(7);
        
        M = params.getQuick(13);
        initSFocal = params.getQuick(15);
        
    }
    
    public DoubleArrayList getInitialConditions() 
    {
        
        //Order init conditions 1) Sf 2) Sg 3) If 4) Ig 5) Nf 6) Ng
        
        DoubleArrayList currInits = new DoubleArrayList();
       
        //Compute equilibrium conditons without migration
        //double R0Focal = betaFF / (mu + nu); //assuming no migration and no vector
        double R0Focal = (betaHV*betaVH * M) / (muVec * (mu + nu));
        double betaFocal = R0Focal * (mu + nu); 
        double R0Global = betaGG / (mu + nu); //assuming no migration
        
        double initS1 = NFocal/R0Focal;
        double initS2 = NGlobal/R0Global;
        double initI1 = mu * NFocal * (R0Focal - 1) / (betaFocal); 
        double initI2 = mu * NGlobal * (R0Global - 1) / betaGG;
        
        double initSv = (muVec * M * NFocal) / ((betaHV/NFocal) * initI1 + muVec);
        double initIv = (betaHV * initSv * (initI1/NFocal)) / muVec;
        
        currInits.add(initS1);
        currInits.add(initS2);
        currInits.add(initSv);
        currInits.add(initI1);
        currInits.add(initI2);
        currInits.add(initIv);
        currInits.add(NFocal);
        currInits.add(NGlobal);
        
        double dt = 0.25;
        double length = 365.25*200;
        currInits = this.fastSim(currInits, dt, length);
        
        //If estimating init S1
        initS1 = initSFocal * NFocal; currInits.set(0, initS1); //estimating Sf init
        
        //True inits for VectorDengueTree_M1R3_062513.tre at 14 years
        //double initS1 = 666807.299; currInits.add(initS1);
        //double initS2 = 1320132.38209387; currInits.add(initS2);
        //double initSv = 1839394.279; currInits.add(initSv);
        //double initI1 = 368.142; currInits.add(initI1);
        //double initI2 = 93.4958120618454; currInits.add(initI2);
        //double initIv = 611.639; currInits.add(initIv);
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
        //if (params.get(15) > 1.0) {
            //constraintCheckFail = true;
        //}
        
        //deltaFocal has to be less than 365.25
        if (params.get(7) > (1.0)) { //this means delta is less than one year
            constraintCheckFail = true;
        }
        
        
        return constraintCheckFail;
    }
    
    //Method not fully updated for this model
    /**
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
    } */
    
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
            //B = muVec * M * NFocal;
            seasM = M * (1 + alphaFocal*seasHeightFocal);
            
            //Compute betaSeas for global pop
            seasNowGlobal = ((timeNow/365.25) + deltaGlobal) - Math.floor(timeNow/365.25);
            seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
            betaSeasGG = betaGG * (1 + alphaGlobal*seasHeightGlobal);
            betaSeasFG = betaFG * (1 + alphaGlobal*seasHeightGlobal);
            
            for (int j = 0; j < jParticles; ++j) {
                
                //Get quasi-equilibrium vector states
                //double vecI = seasM * NFocal * (betaFF/NFocal) * muVec * matrix.getQuick(2,j,xIndex) / ((betaFF/NFocal)*matrix.getQuick(2,j,xIndex)*muVec + muVec*muVec);
                //double vecI = (seasB * betaHV * (matrix.getQuick(2,j,xIndex) / NFocal)) / ((betaHV/NFocal)*matrix.getQuick(2,j,xIndex)*muVec + muVec*muVec);

                //INFECTIONS IN FOCAL HUMANS POP FROM FOCAL VECTORS POP
                //dSfIf = (betaSeasFF * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));
                dSfIf = (betaVH * matrix.getQuick(0,j,xIndex) * matrix.getQuick(5,j,xIndex) * dt / matrix.getQuick(6,j,xIndex));
                
                //INFECTIONS IN FOCAL VECTORS FROM FOCAL HUMANS
                double dSvIf = (betaHV * matrix.getQuick(2,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(6,j,xIndex));
                
                //INFECTIONS IN FOCAL POP FROM GLOBAL POP
                dSfIg = (betaSeasGF * matrix.getQuick(0,j,xIndex) * matrix.getQuick(4,j,xIndex) * dt / matrix.getQuick(6,j,xIndex));

                //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
                dSgIg = (betaSeasGG * matrix.getQuick(1,j,xIndex) * matrix.getQuick(4,j,xIndex) * dt / matrix.getQuick(7,j,xIndex));

                //INFECTIONS IN GLOBAL POP FROM FOCAL POP
                dSgIf = (betaSeasFG * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(7,j,xIndex));
                
                //RECOVERIES
                dIfRf = (nu * matrix.getQuick(3,j,xIndex) * dt);
                dIgRg = (nu * matrix.getQuick(4,j,xIndex) * dt);
                
                //BIRTHS
                dNfB = (mu * matrix.getQuick(6,j,xIndex) * dt);
                dNgB = (mu * matrix.getQuick(7,j,xIndex) * dt);
                
                //DEATHS IN HUMANS
                Rf = matrix.getQuick(6,j,xIndex) - matrix.getQuick(0,j,xIndex) - matrix.getQuick(3,j,xIndex);
                Rg = matrix.getQuick(7,j,xIndex) - matrix.getQuick(1,j,xIndex) - matrix.getQuick(4,j,xIndex);
                
                dSfD = (mu * matrix.getQuick(0,j,xIndex) * dt);
                dIfD = (mu * matrix.getQuick(3,j,xIndex) * dt);
                dRfD = (mu * Rf * dt);
                
                dSgD = (mu * matrix.getQuick(1,j,xIndex) * dt);
                dIgD = (mu * matrix.getQuick(4,j,xIndex) * dt);
                dRgD = (mu * Rg * dt);
                
                //DEATHS IN MOSQUITOES
                double dIvD = muVec * matrix.getQuick(5,j,xIndex) * dt;
                
                //UPDATES
                newSf = matrix.getQuick(0,j,xIndex) + dNfB - dSfIf - dSfIg - dSfD;
                newSg = matrix.getQuick(1,j,xIndex) + dNgB - dSgIg - dSgIf - dSgD;
                newIf = matrix.getQuick(3,j,xIndex) + dSfIf + dSfIg - dIfRf - dIfD;
                if (newIf < 0.0) {
                    newIf = 0.0;
                }
                newIg = matrix.getQuick(4,j,xIndex) + dSgIg + dSgIf - dIgRg - dIgD;
                if (newIg < 0.0) {
                    newIg = 0.0;
                }
                newRf = Rf + dIfRf - dRfD;
                newRg = Rg + dIgRg - dRgD;
                newNf = newSf + newIf + newRf;
                newNg = newSg + newIg + newRg;
                
                double newIv = matrix.getQuick(5,j,xIndex) + dSvIf - dIvD;
                double newNv = seasM * matrix.getQuick(6,j,xIndex);
                double newSv = newNv - newIv;
                
                matrix.setQuick(0,j,nextIndex,newSf);
                matrix.setQuick(1,j,nextIndex,newSg);
                matrix.setQuick(2,j,nextIndex,newSv);
                matrix.setQuick(3,j,nextIndex,newIf);
                matrix.setQuick(4,j,nextIndex,newIg);
                matrix.setQuick(5,j,nextIndex,newIv);
                matrix.setQuick(6,j,nextIndex,newNf);
                matrix.setQuick(7,j,nextIndex,newNg);
                matrix.setQuick(8,j,nextIndex,dSfIf);
                matrix.setQuick(9,j,nextIndex,dSfIg);
                matrix.setQuick(10,j,nextIndex,dSgIg);
                matrix.setQuick(11,j,nextIndex,dSgIf);
                matrix.setQuick(12,j,nextIndex,dIfRf);
                matrix.setQuick(13,j,nextIndex,dIgRg);
                matrix.setQuick(14,j,nextIndex,dNfB);
                matrix.setQuick(15,j,nextIndex,dNgB);
                matrix.setQuick(16,j,nextIndex,dSfD);
                matrix.setQuick(17,j,nextIndex,dIfD);
                matrix.setQuick(18,j,nextIndex,dRfD);
                matrix.setQuick(19,j,nextIndex,dSgD);
                matrix.setQuick(20,j,nextIndex,dIgD);
                matrix.setQuick(21,j,nextIndex,dRgD);
                
                //Track cumulative incidence
                matrix.setQuick(22,j,nextIndex, (dSfIf + dSfIg));
                matrix.setQuick(23,j,nextIndex, (dSgIg + dSgIf));
                
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
        x.add(inits.get(6));
        x.add(inits.get(7));
        
	for (int n = 0; n < times.size(); ++n) {
            
            timeNow = times.getQuick(n);
            
            //Compute betaSeas for focal pop
            seasNowFocal = ((timeNow/365.25) + deltaFocal) - Math.floor(timeNow/365.25);
            seasHeightFocal = Math.cos(2*Math.PI * seasNowFocal);
            //betaSeasFF = betaFF * (1 + alphaFocal*seasHeightFocal);
            betaSeasGF = betaGF * (1 + alphaFocal*seasHeightFocal);
            
            //If estimating seasonal M
            //B = muVec * M * NFocal;
            //seasB = B * (1 + alphaFocal*seasHeightFocal);
            seasM = M * (1 + alphaFocal*seasHeightFocal);
            
            //Compute betaSeas for global pop
            seasNowGlobal = ((timeNow/365.25) + deltaGlobal) - Math.floor(timeNow/365.25);
            seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
            betaSeasGG = betaGG * (1 + alphaGlobal*seasHeightGlobal);
            betaSeasFG = betaFG * (1 + alphaGlobal*seasHeightGlobal);            
            
            //INFECTIONS IN FOCAL HUMANS POP FROM FOCAL VECTORS POP
            //dSfIf = (betaSeasFF * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));
            dSfIf = (betaVH * x.getQuick(0) * x.getQuick(5) * dt / x.getQuick(6));

            //INFECTIONS IN FOCAL VECTORS FROM FOCAL HUMANS
            double dSvIf = (betaHV * x.getQuick(2) * x.getQuick(3) * dt / x.getQuick(6));

            //INFECTIONS IN FOCAL POP FROM GLOBAL POP
            dSfIg = (betaSeasGF * x.getQuick(0) * x.getQuick(4) * dt / x.getQuick(6));

            //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
            dSgIg = (betaSeasGG * x.getQuick(1) * x.getQuick(4) * dt / x.getQuick(7));

            //INFECTIONS IN GLOBAL POP FROM FOCAL POP
            dSgIf = (betaSeasFG * x.getQuick(1) * x.getQuick(3) * dt / x.getQuick(7));

            //RECOVERIES
            dIfRf = (nu * x.getQuick(3) * dt);
            dIgRg = (nu * x.getQuick(4) * dt);

            //BIRTHS
            dNfB = (mu * x.getQuick(6) * dt);
            dNgB = (mu * x.getQuick(7) * dt);

            //DEATHS IN HUMANS
            Rf = x.getQuick(6) - x.getQuick(0) - x.getQuick(3);
            Rg = x.getQuick(7) - x.getQuick(1) - x.getQuick(4);

            dSfD = (mu * x.getQuick(0) * dt);
            dIfD = (mu * x.getQuick(3) * dt);
            dRfD = (mu * Rf * dt);

            dSgD = (mu * x.getQuick(1) * dt);
            dIgD = (mu * x.getQuick(4) * dt);
            dRgD = (mu * Rg * dt);

            //DEATHS IN MOSQUITOES
            double dIvD = muVec * x.getQuick(5) * dt;

            //UPDATES
            newSf = x.getQuick(0) + dNfB - dSfIf - dSfIg - dSfD;
            newSg = x.getQuick(1) + dNgB - dSgIg - dSgIf - dSgD;
            newIf = x.getQuick(3) + dSfIf + dSfIg - dIfRf - dIfD;
            if (newIf < 0.0) {
                newIf = 0.0;
            }
            newIg = x.getQuick(4) + dSgIg + dSgIf - dIgRg - dIgD;
            if (newIg < 0.0) {
                newIg = 0.0;
            }
            newRf = Rf + dIfRf - dRfD;
            newRg = Rg + dIgRg - dRgD;
            newNf = newSf + newIf + newRf;
            newNg = newSg + newIg + newRg;

            double newIv = x.getQuick(5) + dSvIf - dIvD;
            double newNv = seasM * x.getQuick(6);
            double newSv = newNv - newIv;

            x.setQuick(0,newSf);
            x.setQuick(1,newSg);
            x.setQuick(2,newSv);
            x.setQuick(3,newIf);
            x.setQuick(4,newIg);
            x.setQuick(5,newIv);
            x.setQuick(6,newNf);
            x.setQuick(7,newNg);
	
        }
        return x;
    }
    
    //Need to update this method for this model
    /**
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
    } */

       
}
