import cern.colt.matrix.*;
import cern.colt.list.*;
import java.util.ArrayList;

/**
 *
 * @author David
 */
public class StructCoalModelAll {
    
    DoubleMatrix2D F;
    DoubleMatrix2D G;
    DoubleMatrix1D Y;
    int states;
    int transCount;
    int nonTransCount;
    ArrayList<ArrayList<Integer>> transCountArray = new ArrayList<ArrayList<Integer>>();
    double currTime;

    
    public void make() 
    {   
        
        /**
         * 
         * 
         * NEED TO UPDATE THIS FOR EACH MODEL!!
         * 
         */
        
        states = 2; 
        
        DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
        F = factory2D.make(states, states);
        F.assign(0);
        G = factory2D.make(states, states);
        G.assign(0);
        
        DoubleFactory1D factory1D;
        factory1D = DoubleFactory1D.dense;
        Y = factory1D.make(states);
        Y.assign(0);
    }
    
    /** For SEIR model
    public void updateF(DoubleMatrix1D xCurr, DoubleArrayList theta) 
    {
        
        double beta = theta.get(0) * theta.get(3);
        double totalN = xCurr.get(3);
        double fIE = beta * xCurr.get(0) * xCurr.get(2) / totalN;
        
        F.setQuick(0,0,0.0);
        F.setQuick(0,1,0.0);
        F.setQuick(1,0,fIE);
        F.setQuick(1,1,0.0);
    }
    
    public void updateG(DoubleMatrix1D xCurr, DoubleArrayList theta)
    {
        //theta.get(1) ==  sigma
        double gEI = theta.get(1) * xCurr.get(1); //E * sigma
        
        G.setQuick(0,0,0.0);
        G.setQuick(0,1,gEI);
        G.setQuick(1,0,0.0);
        G.setQuick(1,1,0.0);
    }
    
    public void updateY(DoubleMatrix1D xCurr) 
    {
        Y.setQuick(0,xCurr.get(1)); //infected juveniles
        Y.setQuick(1,xCurr.get(2)); //infected adults    
    }*/
     
    
     
    
    /** For age-structured model
    public void updateF(DoubleMatrix1D xCurr, DoubleArrayList theta) 
    {
        
        double totalN = xCurr.get(4) + xCurr.get(5);
        double fJJ = theta.get(0) * xCurr.get(0) * xCurr.get(2) / totalN; //rate of transmission from J -> J
        double fJA = theta.get(1) * xCurr.get(1) * xCurr.get(2) / totalN; //rate of transmission from J -> A
        double fAJ = theta.get(1) * xCurr.get(0) * xCurr.get(3) / totalN;
        double fAA = theta.get(1) * xCurr.get(1) * xCurr.get(3) / totalN;
        
        F.setQuick(0,0,fJJ);
        F.setQuick(0,1,fJA);
        F.setQuick(1,0,fAJ);
        F.setQuick(1,1,fAA);
    }
    
    public void updateG(DoubleMatrix1D xCurr, DoubleArrayList theta)
    {
        //G rates are always zero
    }
    
    public void updateY(DoubleMatrix1D xCurr) 
    {
        Y.setQuick(0,xCurr.get(2)); //infected juveniles
        Y.setQuick(1,xCurr.get(3)); //infected adults    
    } */
    
    
    
    
    
    /** Semi-spatial model with global population at endemic equilbrium
    public void updateF(DoubleMatrix1D xCurr, DoubleArrayList theta, double currTime) 
    {
        
        double beta = theta.get(0);
        double alpha = theta.get(1);
        double delta = theta.get(2);
        double M = theta.get(3);
        double IgEff = xCurr.get(3);
        double mu = theta.get(6);
        double nu = theta.get(7);
        
        //For peice-wise seasonal forcing with a high and low season
        //double currTimeYear = (currTime / 365.25) - Math.floor(currTime/365.25);
        //double startHighSeason = delta / 365.25;
        //double endHighSeason = (delta + (365.25/2)) / 365.25;
        //double betaNow = 0.0;
        //if (endHighSeason < 1.0) {
            //if (currTimeYear > startHighSeason & currTimeYear < endHighSeason) {
                //betaNow = betaHigh;
            //} else {
                //betaNow = betaLow;
            //}
        //} else {
            //if (currTimeYear < endHighSeason | currTimeYear > startHighSeason) {
                //betaNow = betaHigh;
            //} else {
                //betaNow = betaLow;
            //}
        //}
        
        //For sinuisoidal forcing
        double seasNow = ((currTime + delta)/365.25) - Math.floor(currTime/365.25);
        double seasHeight = Math.cos(2*Math.PI * seasNow);
        if (seasHeight < 0) {
            seasHeight = 0;
        }
        double betaNow = beta * (1 + alpha*seasHeight);
        
        double fFF = betaNow * xCurr.get(0) * xCurr.get(1) / xCurr.get(2); //rate of transmission within focal pop
        double fGF = betaNow * M * xCurr.get(0) * IgEff / xCurr.get(2); //rate of transmission from J -> A
        double fGG = IgEff * (nu + mu);
        
        //fFG adds export from the focal pop
        double R0 = beta/(mu+nu);
        double NgEff = IgEff / (mu * ((R0 - 1)/beta));
        double SgEff = NgEff / R0; 
        double fFG = beta * M * SgEff * xCurr.get(1) / NgEff;
        //double fFG = 0;
        
        F.setQuick(0,0,fFF);
        F.setQuick(0,1,fFG);
        F.setQuick(1,0,fGF);
        F.setQuick(1,1,fGG);
    }
    
    public void updateG(DoubleMatrix1D xCurr, DoubleArrayList theta)
    {
        //G Matrix is always zero
    }
    
    public void updateY(DoubleMatrix1D xCurr) 
    {
        Y.setQuick(0,xCurr.get(1)); //infected juveniles
        Y.setQuick(1,xCurr.get(3)); //infected adults    
    }*/
    
    

    
    //Full spatial SIR model with seasonality in focal and global pops
    /**
    public void updateF(DoubleMatrix1D xCurr, DoubleArrayList theta, double currTime) 
    {
        
        //double totalN = xCurr.get(3);
        double betaFF = theta.getQuick(0);
        double betaGF = theta.getQuick(1);
        double betaGG = theta.getQuick(2);
        double betaFG = theta.getQuick(3);
        double alphaFocal = theta.get(4);
        double alphaGlobal = theta.get(5);
        double deltaFocal = theta.get(6);
        double deltaGlobal = theta.get(7);
        
        //Compute betaSeas for focal pop
        double seasNowFocal = ((currTime/365.25) + deltaFocal) - Math.floor(currTime/365.25);
        double seasHeightFocal = Math.cos(2*Math.PI * seasNowFocal);
        double betaSeasFF = betaFF * (1 + alphaFocal*seasHeightFocal);
        double betaSeasGF = betaGF * (1 + alphaFocal*seasHeightFocal); 

        //Compute betaSeas for global pop
        double seasNowGlobal = ((currTime/365.25) + deltaGlobal) - Math.floor(currTime/365.25);
        double seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
        double betaSeasGG = betaGG * (1 + alphaGlobal*seasHeightGlobal);
        double betaSeasFG = betaFG * (1 + alphaGlobal*seasHeightGlobal);
        
        
        double fFF = betaSeasFF * xCurr.get(0) * xCurr.get(2) / xCurr.get(4); //rate of transmission within focal pop
        double fFG = betaSeasFG * xCurr.get(1) * xCurr.get(2) / xCurr.get(5);
        double fGF = betaSeasGF * xCurr.get(0) * xCurr.get(3) / xCurr.get(4); //rate of transmission from J -> A
        double fGG = betaSeasGG * xCurr.get(1) * xCurr.get(3) / xCurr.get(5);
        
        F.setQuick(0,0,fFF);
        F.setQuick(0,1,fFG);
        F.setQuick(1,0,fGF);
        F.setQuick(1,1,fGG);
    }*/
    
    
    
    
    //Vectored spatial SIR model with seasonality in focal and global pops
    /**
    public void updateF(DoubleMatrix1D xCurr, DoubleArrayList theta, double currTime) 
    {
        
        //double totalN = xCurr.get(3);
        double betaFF = theta.getQuick(0);
        double betaGF = theta.getQuick(1);
        double betaGG = theta.getQuick(2);
        double betaFG = theta.getQuick(3);
        double alphaFocal = theta.get(4);
        double alphaGlobal = theta.get(5);
        double deltaFocal = theta.get(6);
        double deltaGlobal = theta.get(7);
        double M = theta.get(12);
        double muVec = theta.get(13);
        
        //Compute betaSeas for focal pop
        double seasNowFocal = ((currTime/365.25) + deltaFocal) - Math.floor(currTime/365.25);
        double seasHeightFocal = Math.cos(2*Math.PI * seasNowFocal);
        //double betaSeasFF = betaFF * (1 + alphaFocal*seasHeightFocal);
        double betaSeasGF = betaGF * (1 + alphaFocal*seasHeightFocal);
        
        double seasM = M * (1 + alphaFocal*seasHeightFocal);

        //Compute betaSeas for global pop
        double seasNowGlobal = ((currTime/365.25) + deltaGlobal) - Math.floor(currTime/365.25);
        double seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
        double betaSeasGG = betaGG * (1 + alphaGlobal*seasHeightGlobal);
        double betaSeasFG = betaFG * (1 + alphaGlobal*seasHeightGlobal);
        
        //double fFF = betaSeasFF * xCurr.get(0) * xCurr.get(2) / xCurr.get(4); //rate of transmission within focal pop
        double vecI = seasM * xCurr.get(4) * (betaFF/xCurr.get(4)) * muVec * xCurr.get(2) / ((betaFF/xCurr.get(4))*xCurr.get(2)*muVec + muVec*muVec);
        double vecS = seasM * xCurr.get(4) * muVec / ((betaFF/xCurr.get(4)) * xCurr.get(2) + muVec);
        double fVH = betaFF * (xCurr.get(0)/xCurr.get(4)) * vecI;
        double fHV = betaFF * (vecS/xCurr.get(4)) * xCurr.get(2);
        double totalI = vecI + xCurr.get(2);
        //double stateProbTerm = 2 * (vecI/totalI) * (xCurr.get(2)/totalI); 
        double fFF = ((fVH + fHV) * xCurr.get(2) * xCurr.get(2)) / (totalI*totalI);
        double fFG = betaSeasFG * xCurr.get(1) * xCurr.get(2) / xCurr.get(5);
        double fGF = betaSeasGF * xCurr.get(0) * xCurr.get(3) / xCurr.get(4); //rate of transmission from J -> A
        double fGG = betaSeasGG * xCurr.get(1) * xCurr.get(3) / xCurr.get(5);
        
        F.setQuick(0,0,fFF);
        F.setQuick(0,1,fFG);
        F.setQuick(1,0,fGF);
        F.setQuick(1,1,fGG);
    } 
     
    public void updateG(DoubleMatrix1D xCurr, DoubleArrayList theta)
    {
        //G Matrix is always zero
    }
    
    public void updateY(DoubleMatrix1D xCurr) 
    {
        Y.setQuick(0,xCurr.get(2)); //infections in focal pop
        Y.setQuick(1,xCurr.get(3)); //infected in global pop    
    } 
    */
    

}
