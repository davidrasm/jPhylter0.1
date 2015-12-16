import java.util.ArrayList;
import cern.colt.matrix.*;
import cern.colt.list.DoubleArrayList;
import java.io.*; //IOException;

/**
 *
 * @author David
 */
public class AutoOptimizer {
    
    boolean autoOptimize;
    ArrayList<Integer> runPoints = new ArrayList<Integer>();
    DoubleMatrix2D samples;
    int sampleParams;
    int currSample;
    int currStage;
    
    
    public void setRunPoints(ArrayList<Integer> points, int estParams)
    {
        sampleParams = estParams;
        runPoints.add(0);
        runPoints.addAll(points);
        currStage = 0;
    }
    
    public void getSampleMatrix() {
        
        int nextRunPoint = runPoints.get(currStage + 1);
        int lastRunPoint = runPoints.get(currStage);
        int nSamples = nextRunPoint - lastRunPoint;
        DoubleFactory2D factory2D;
        factory2D = DoubleFactory2D.dense;
        samples = factory2D.make(nSamples, sampleParams);
        currSample = 0;
        
    }
    
    public void addSample(DoubleArrayList sample) {
    
        int sampleSize = sample.size();
        for(int n = 0; n < sampleSize; n++) {
            samples.set(currSample, n, sample.getQuick(n));
        }
        currSample++;
        
    }
    
    public DoubleMatrix2D autoOptimize() throws IOException
    { 
        DoubleMatrix2D covMatrix = MatrixUtils.getCovMatrix(samples);
        currStage++;
        System.out.println();
        if (currStage < (runPoints.size() - 1)) {
            this.getSampleMatrix();
        } else {
            //Terminate autoOptimization and print out final covMatrix
            //WriteMatrixToCSV csvWriter = new WriteMatrixToCSV();
            //String covMatrixFile = "PHYLTerFinalProposalCov";
            //csvWriter.writeCSV(covMatrix, covMatrixFile);
            autoOptimize = false;
            
        }
        return covMatrix;
        
    }
    
    
}
