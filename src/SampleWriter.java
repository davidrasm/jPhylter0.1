import cern.colt.list.*;
//import java.util.ArrayList;
//import java.util.Collections;
/**
 *
 * @author David
 */
public class SampleWriter {
    
    String sampleString;

    public void samplesToString(DoubleArrayList theta) 
    {
        
        int cols = theta.size();
        sampleString = "";
        for (int param = 0; param < (cols - 1); ++param) {
            sampleString = sampleString + Double.toString(theta.getQuick(param)) + ", ";
            //thetaSamples.set(sampleLoc, param, thetaEstNow.get(param));
        }
        sampleString = sampleString + Double.toString(theta.getQuick(cols - 1));
        //System.out.println();
    }
}
