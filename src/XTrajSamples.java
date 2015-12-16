import cern.colt.matrix.*;
import cern.colt.list.*;
import java.util.ArrayList;

/**
 *
 * @author David
 */
public class XTrajSamples {
    
    DoubleArrayList sample;
    DoubleArrayList sampleTimes;
    
    public void getSampleTimes(double startTime, double endTime, double sampleFreq) 
    {    
        sampleTimes = new DoubleArrayList();
        sampleTimes.add(startTime);
        double timeNow = startTime + sampleFreq;
        while (timeNow < endTime) {
            sampleTimes.add(timeNow);
            timeNow = timeNow + sampleFreq;    
        }
        if (sampleTimes.contains(endTime)) {
            //Don't add it
        } else {
            sampleTimes.add(endTime);
        } 
    }
    
    public DoubleArrayList getSample(int stateVar, ZVectors dataZ, XTrajectory xTraj) 
    {    
        sample = new DoubleArrayList();
        int thisZLoc; int xSampleLoc; double xTimeSample;
        for (int n = 0; n < sampleTimes.size(); ++n) {
            thisZLoc = dataZ.absoluteTimes.indexOf(sampleTimes.get(n));
            xSampleLoc = thisZLoc - xTraj.zLocStart;
            xTimeSample = xTraj.full.getQuick(stateVar, xSampleLoc);
            sample.add(xTimeSample);   
        }
        return sample;
        //System.out.println();
    }
    
    public DoubleArrayList getSampleIncidence(int stateVar, ZVectors dataZ, XTrajectory xTraj) 
    {    
        sample = new DoubleArrayList();
        double sumIncidence;
        int thisZLoc = 0; int prevZLoc = 0; int xSampleLoc = 0;
        
        for (int n = 0; n < sampleTimes.size(); ++n) {
            
            sumIncidence = 0.0;
            thisZLoc = dataZ.absoluteTimes.indexOf(sampleTimes.get(n));
            if (n > 0) {
                prevZLoc = dataZ.absoluteTimes.indexOf(sampleTimes.get(n-1)) + 1;
            } else {
                prevZLoc = dataZ.absoluteTimes.indexOf(dataZ.startTime);
            }
            for (int t = prevZLoc; t <= thisZLoc; t++) {
                xSampleLoc = t - xTraj.zLocStart;
                sumIncidence += xTraj.full.getQuick(stateVar, xSampleLoc);
            }
            //System.out.println();
            
            //intervalIncidences = xTraj.full.g
            //int xSampleLoc = thisZLoc - xTraj.zLocStart;
            //double xTimeSample = xTraj.full.getQuick(stateVar, xSampleLoc);
            sample.add(sumIncidence);
            
        }
        return sample;
        //System.out.println();
    }
    
}
