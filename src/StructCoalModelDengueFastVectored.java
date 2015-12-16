import cern.colt.matrix.*;
import cern.colt.list.*;

/**
 *
 * @author David
 */
public class StructCoalModelDengueFastVectored extends StructCoalModel {
    
    @Override
    public void make() 
    {   
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
    
    @Override
    public void updateF(DoubleMatrix1D xCurr, DoubleArrayList theta, double currTime) 
    {
        
        //double totalN = xCurr.get(3);
        double betaHV = theta.getQuick(0);
        double betaVH = theta.getQuick(1);
        double betaGF = theta.getQuick(2);
        double betaGG = theta.getQuick(3);
        double betaFG = theta.getQuick(4);
        double alphaFocal = theta.get(5);
        double alphaGlobal = theta.get(6);
        double deltaFocal = theta.get(7);
        double deltaGlobal = theta.get(8);
        double M = theta.get(13);
        double muVec = theta.get(14);
        
        //Compute betaSeas for focal pop
        double seasNowFocal = ((currTime/365.25) + deltaFocal) - Math.floor(currTime/365.25);
        double seasHeightFocal = Math.cos(2*Math.PI * seasNowFocal);
        //double betaSeasFF = betaFF * (1 + alphaFocal*seasHeightFocal);
        double betaSeasGF = betaGF * (1 + alphaFocal*seasHeightFocal);
        
        double B = muVec * M * xCurr.get(4);
        double seasB = B * (1 + alphaFocal*seasHeightFocal);

        //Compute betaSeas for global pop
        double seasNowGlobal = ((currTime/365.25) + deltaGlobal) - Math.floor(currTime/365.25);
        double seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
        double betaSeasGG = betaGG * (1 + alphaGlobal*seasHeightGlobal);
        double betaSeasFG = betaFG * (1 + alphaGlobal*seasHeightGlobal);
        
        //double fFF = betaSeasFF * xCurr.get(0) * xCurr.get(2) / xCurr.get(4); //rate of transmission within focal pop
        double vecI = (seasB * betaHV * (xCurr.getQuick(2) / xCurr.get(4))) / ((betaHV/xCurr.get(4))*xCurr.get(2)*muVec + muVec*muVec);
        double vecS = seasB / ((betaHV/xCurr.get(4)) * xCurr.get(2) + muVec);
        //double vecI = seasM * xCurr.get(4) * (betaFF/xCurr.get(4)) * muVec * xCurr.get(2) / ((betaFF/xCurr.get(4))*xCurr.get(2)*muVec + muVec*muVec);
        //double vecS = seasM * xCurr.get(4) * muVec / ((betaFF/xCurr.get(4)) * xCurr.get(2) + muVec);
        double fVH = betaVH * (xCurr.get(0)/xCurr.get(4)) * vecI;
        double fHV = betaHV * (vecS/xCurr.get(4)) * xCurr.get(2);
        double totalI = vecI + xCurr.get(2);
        //double stateProbTerm = 2 * (vecI/totalI) * (xCurr.get(2)/totalI); 
        
        
        //Old way of computing coalescent rate
        //double fFF = ((fVH + fHV) * xCurr.get(2) * xCurr.get(2)) / (totalI*totalI);
        
        //New way with correct lineage state probs
        double humI = xCurr.get(2);
        double stateProbTerm = ((fVH/humI)*(fHV/vecI))/(Math.pow(((fVH/humI)+(fHV/vecI)),2)); //otherwise multiplied by 2
        double fFF = ((fVH + fHV) / (vecI*humI)) * stateProbTerm * (humI * humI); 
        
        double fFG = betaSeasFG * xCurr.get(1) * xCurr.get(2) / xCurr.get(5);
        double fGF = betaSeasGF * xCurr.get(0) * xCurr.get(3) / xCurr.get(4); //rate of transmission from J -> A
        double fGG = betaSeasGG * xCurr.get(1) * xCurr.get(3) / xCurr.get(5);
        
        F.setQuick(0,0,fFF);
        F.setQuick(0,1,fFG);
        F.setQuick(1,0,fGF);
        F.setQuick(1,1,fGG);
    } 
    
    @Override
    public void updateG(DoubleMatrix1D xCurr, DoubleArrayList theta)
    {
        //G Matrix is always zero
    }
    
    @Override
    public void updateY(DoubleMatrix1D xCurr) 
    {
        Y.setQuick(0,xCurr.get(2)); //infections in focal pop
        Y.setQuick(1,xCurr.get(3)); //infected in global pop    
    }
    
}
