import cern.colt.matrix.*;
import cern.colt.list.*;
import cern.colt.matrix.linalg.*;
import cern.jet.stat.Descriptive;

public class MatrixUtils
{

    //Matrix exponentiation
    public static DoubleMatrix2D exp(DoubleMatrix2D A)
    {
        int size = A.rows();
        EigenvalueDecomposition eigen = new EigenvalueDecomposition(A);
	DoubleMatrix2D V = eigen.getV();
	DoubleMatrix2D invV = Algebra.DEFAULT.inverse(V);
	DoubleMatrix2D D = eigen.getD();
		
	DoubleMatrix2D expD = D.copy();
	for (int j = 0; j < size; ++j) {
            expD.setQuick(j,j,(Math.exp(D.getQuick(j,j))));
	}
		
	DoubleMatrix2D E = Algebra.DEFAULT.mult((Algebra.DEFAULT.mult(V,expD)),invV);
	return E;	
    }
    
    public static DoubleMatrix2D getCovMatrix(DoubleMatrix2D data)
    {
        int dataCols = data.columns();
        
        DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
        DoubleMatrix2D covMatrix = factory2D.make(dataCols, dataCols); 
       
        for (int i = 0; i < dataCols; ++i) {
            for (int j = i; j < dataCols; ++j) {
                
                DoubleMatrix1D data1Matrix = data.viewColumn(i);
                DoubleMatrix1D data2Matrix = data.viewColumn(j);
                double [] data1Array = data1Matrix.toArray();
                double [] data2Array = data2Matrix.toArray();
                DoubleArrayList data1List = new DoubleArrayList();
                DoubleArrayList data2List = new DoubleArrayList();
                data1List.elements(data1Array);
                data2List.elements(data2Array);
                
                double pairCov = Descriptive.covariance(data1List, data2List);
                
                covMatrix.set(i,j,pairCov);
                if (i != j) {
                    covMatrix.set(j,i,pairCov);
                }                
            }
        }
        
        return covMatrix;
        
    }
    
    public static double getMultiVarNormProb(DoubleMatrix1D x, DoubleMatrix1D mu, DoubleMatrix2D sigma) 
    {
        double prob = 1.0;
        double detSigma = Algebra.DEFAULT.det(sigma);
        if (detSigma == 0.0) { //cov matrix is singular
            prob = 0.0;
        } else {
        
            int dim = x.size();
            DoubleMatrix1D diff = x.like();
            for (int n = 0; n < dim; n++) {
                diff.setQuick(n, (x.getQuick(n) - mu.getQuick(n)));
            }
            DoubleMatrix2D invSigma = Algebra.DEFAULT.inverse(sigma);
            DoubleMatrix1D sigmaDiff = Algebra.DEFAULT.mult(invSigma, diff);
            double diffSigmaDiff = Algebra.DEFAULT.mult(diff, sigmaDiff);
            double exponent = Math.exp(-0.5 * diffSigmaDiff);
            double denom = Math.sqrt(Math.pow((2*Math.PI),dim) * detSigma);
            prob = exponent / denom;
        }
        
        return prob;
    }
    
}
