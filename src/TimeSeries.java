import cern.colt.list.*;
//import cern.colt.matrix.*;
import java.util.ArrayList;

/**
 *
 * @author David
 */
public class TimeSeries {
    
    boolean active;
    DoubleArrayList times;
    DoubleArrayList observations;
    
    public void setTimeSeries(ArrayList<ArrayList<String>> epiData)
    {
        int totalTimes = epiData.get(0).size();
        times = new DoubleArrayList();
        observations = new DoubleArrayList();
        for (int n = 0; n < totalTimes; n++) {
            times.add(Double.valueOf(epiData.get(0).get(n)));
            observations.add(Double.valueOf(epiData.get(1).get(n)));
        }
    }
    
}
