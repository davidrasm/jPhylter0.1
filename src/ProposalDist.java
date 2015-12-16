import cern.colt.list.*;
import cern.colt.matrix.*;
import cern.jet.random.*;
import cern.colt.matrix.linalg.*;
import cern.colt.matrix.linalg.Algebra;

/**
 * modify nextProposal method to change proposal density (currently MVN)
 * @author David
 */

public class ProposalDist {
    
    DoubleMatrix2D covMatrix;
    AbstractDistribution stndNorm;
    
    public void setRandomGenerator() 
    {
      
        cern.jet.random.engine.RandomEngine engine = new cern.jet.random.engine.MersenneTwister(new java.util.Date());
        stndNorm = new Normal(0.0, 1.0, engine);
        
    }
    
    public void setCov(double[][] covArray, int numParams)
    {
        
        DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
	covMatrix = factory2D.make(numParams, numParams);
        
        for (int j = 0; j < numParams; ++j) {
            for (int k = 0; k < numParams; ++k) {
                covMatrix.set(j,k,covArray[j][k]);
            }
        }
        
    }
    
    public DoubleArrayList nextProposal(DoubleArrayList theta, EpiModel epi)
    {
        
        //This is just a multivariate normal proposal density
        
        //Make mu vector
        int thetaSize = theta.size();
	DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
	DoubleMatrix2D thetaVec = factory2D.make(thetaSize,1);
	for(int param = 0; param < thetaSize; ++param) {
            thetaVec.set(param, 0, theta.get(param));
	}
		
	//Get Cholesky Decomposition
	CholeskyDecomposition cholesky = new CholeskyDecomposition(covMatrix);
	DoubleMatrix2D LMatrix = cholesky.getL();
	
        boolean constraintCheckFail = true;
        DoubleArrayList newThetaList = new DoubleArrayList(thetaSize);
        DoubleMatrix2D newTheta = factory2D.make(thetaSize,1);
        
        while (constraintCheckFail) {
            
            constraintCheckFail = false;
            
            //Get uniform random variates
            DoubleMatrix2D uniRands = factory2D.make(thetaSize,1);
            for (int j = 0; j < thetaSize; ++j) {
                uniRands.set(j,0, stndNorm.nextDouble());
            }
	
            //Multiply LMatrix by uniform random variates
            DoubleMatrix2D genVec = factory2D.make(thetaSize,1);
            genVec = Algebra.DEFAULT.mult(LMatrix, uniRands);
	
            //Add means
            double newVal;
            for (int k = 0; k < thetaSize; ++k) {
                newVal = genVec.get(k,0) + theta.get(k);
                newTheta.set(k,0,newVal);
            }
            
            //Rescale delta params
//            int paramIndex = 2; //index in estParams (may be different than index in params)
//            if (newTheta.get(paramIndex,0) < 0.0) {
//                while (newTheta.get(paramIndex,0) < -1.0) {
//                    newTheta.set(paramIndex,0,(newTheta.get(paramIndex,0) + 1));
//                }
//                newTheta.set(paramIndex,0,(1 + newTheta.get(paramIndex,0)));
//            }
//            if (newTheta.get(paramIndex,0) > 1.0) {
//                while (newTheta.get(paramIndex,0) > 2.0) {
//                    newTheta.set(paramIndex,0,(newTheta.get(paramIndex,0) - 1));
//                }
//                newTheta.set(paramIndex,0,(0.0 + (newTheta.get(paramIndex,0) - 1)));
//            }
//            
//            //Rescale delta params
//            paramIndex = 3; //index in estParams (may be different than index in params)
//            if (newTheta.get(paramIndex,0) < 0.0) {
//                while (newTheta.get(paramIndex,0) < -1.0) {
//                    newTheta.set(paramIndex,0,(newTheta.get(paramIndex,0) + 1));
//                }
//                newTheta.set(paramIndex,0,(1 + newTheta.get(paramIndex,0)));
//            }
//            if (newTheta.get(paramIndex,0) > 1.0) {
//                while (newTheta.get(paramIndex,0) > 2.0) {
//                    newTheta.set(paramIndex,0,(newTheta.get(paramIndex,0) - 1));
//                }
//                newTheta.set(paramIndex,0,(0.0 + (newTheta.get(paramIndex,0) - 1)));
//            }
            
            
            newThetaList.removeAll(newThetaList);
            for (int m = 0; m < thetaSize; ++m) {
                newThetaList.add(newTheta.get(m,0));
            }
            
            epi.updateEstParams(newThetaList);
            constraintCheckFail = epi.checkParamConstraints();
               
        }
        
	return newThetaList;
        
    }
    
}
