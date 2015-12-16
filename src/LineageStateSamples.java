import java.util.ArrayList;
import cern.colt.list.*;

/**
 *
 * @author David
 */
public class LineageStateSamples {
    
    ArrayList<ArrayList<ArrayList<Integer>>> array = new ArrayList<ArrayList<ArrayList<Integer>>>();
    ArrayList<ArrayList<ArrayList<Double>>> posteriorProbs = new ArrayList<ArrayList<ArrayList<Double>>>();
    ArrayList<ArrayList<Double>> posteriorProbsSingleState = new ArrayList<ArrayList<Double>>();
    ArrayList<ArrayList<Integer>> discreteProbs = new ArrayList<ArrayList<Integer>>();
    ArrayList<DoubleArrayList> absoluteTimes = new ArrayList<DoubleArrayList>();
    ArrayList<Double> lineageStateArrayTimes = new ArrayList<Double>();
    ArrayList<ArrayList<Integer>> lineageStateArrayMap;
    ArrayList<ArrayList<Double>> entropies = new ArrayList<ArrayList<Double>>();
    ArrayList<ArrayList<Integer>> discreteEntropies = new ArrayList<ArrayList<Integer>>();
    
    int totalTimes;
    int states;
    int lineages;
    int samples = 0;
    
    public void getArray(ParticleLineageStateArray s)
    {
        //Set up sample array
        states = s.states;
        lineages = s.lineages;
        lineageStateArrayTimes = s.absoluteTimes;
        totalTimes = s.absoluteTimes.size();
        lineageStateArrayMap = s.arrayMap;
        
        for (int t = 0; t < totalTimes; t++) {
            array.add(new ArrayList<ArrayList<Integer>>());
            int linesAtTime = s.arrayMap.get(t).size();
            for (int l = 0; l < linesAtTime; l++) {
                array.get(t).add(new ArrayList<Integer>());
                for (int st = 0; st < states; st++) {
                    array.get(t).get(l).add(0);
                }
            }
        }
        
    }
    
    public void addSample(ArrayList<ArrayList<Integer>> sampleArray)
    {
        //Add current sample to sample array
        int linesAtTime = 0;
        int lineageState = 0;
        for (int t = 0; t < totalTimes; t++) {
            linesAtTime = array.get(t).size();
            for (int l = 0; l < linesAtTime; l++) {
                lineageState = sampleArray.get(t).get(l);
                array.get(t).get(l).set(lineageState, (array.get(t).get(l).get(lineageState) + 1));
            }
        }
        samples++;        
    }
    
    public void computeSingleStatePosteriorProbs(TreeNode[] tree, ZVectors dataZ)
    {   
        double samplesAsDouble = (double) samples;
        for (int lin = 0; lin < lineages; lin++) {
            
            int parentNodeIndex = tree[lin].parentNode - 1;
            double lineageStartTime = tree[parentNodeIndex].nodeHeight;
            double lineageEndTime = tree[lin].nodeHeight;
            int lineageStartZIndex = dataZ.absoluteTimes.indexOf(lineageStartTime);
            int lineageEndZIndex = dataZ.absoluteTimes.indexOf(lineageEndTime);
            if (lineageStartZIndex > lineageEndZIndex) {
                lineageStartZIndex = lineageEndZIndex; //this happens at the root
            }
            DoubleArrayList timesForLineage = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(lineageStartZIndex,lineageEndZIndex);
            posteriorProbsSingleState.add(new ArrayList<Double>());
            absoluteTimes.add(timesForLineage);
            
            for (int t = 0; t < timesForLineage.size(); t++) {
                double time = timesForLineage.get(t);
                int lineageStateArrayTimeIndex = lineageStateArrayTimes.indexOf(time);
                if (lineageStateArrayTimeIndex == -1) {
                    //Don't have states for this time
                    double probLineageState = 1.0 / states;
                    //for (int st = 0; st < states; st++) {
                        posteriorProbsSingleState.get(lin).add(probLineageState);
                    //}
                } else {
                    int lineageIndexAtTime = lineageStateArrayMap.get(lineageStateArrayTimeIndex).indexOf(lin);  
                    //for (int st = 0; st < states; st++) {
                        double sampsInState = (double) array.get(lineageStateArrayTimeIndex).get(lineageIndexAtTime).get(0);
                        double probLineageState = sampsInState / samplesAsDouble;
                        posteriorProbsSingleState.get(lin).add(probLineageState);
                    //}
                }
            }
        }
        
    }
    
    public void discretizeSingleStatePosteriorProbs(int classes)
    {
        //For all lineages at all times
        for (int lin = 0; lin < lineages; lin++) {
            discreteProbs.add(new ArrayList<Integer>());
            int lineageTimes = absoluteTimes.get(lin).size();
            for (int time = 0; time < lineageTimes; time++) {
                int probClass = 0;
                double prob = posteriorProbsSingleState.get(lin).get(time);
                if (prob > 0.0) {
                    probClass = (int) Math.ceil(prob * classes); //index from 1!    
                    discreteProbs.get(lin).add(probClass);
                } else {
                    discreteProbs.get(lin).add(1);
                }
            }
        }    
    }
    
    public void computePosteriorProbs(TreeNode[] tree, ZVectors dataZ)
    {   
        double samplesAsDouble = (double) samples;
        for (int lin = 0; lin < lineages; lin++) {
            
            int parentNodeIndex = tree[lin].parentNode - 1;
            double lineageStartTime = tree[parentNodeIndex].nodeHeight;
            double lineageEndTime = tree[lin].nodeHeight;
            int lineageStartZIndex = dataZ.absoluteTimes.indexOf(lineageStartTime);
            int lineageEndZIndex = dataZ.absoluteTimes.indexOf(lineageEndTime);
            if (lineageStartZIndex > lineageEndZIndex) {
                lineageStartZIndex = lineageEndZIndex; //this happens at the root
            }
            DoubleArrayList timesForLineage = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(lineageStartZIndex,lineageEndZIndex);
            posteriorProbs.add(new ArrayList<ArrayList<Double>>());
            absoluteTimes.add(timesForLineage);
            
            for (int t = 0; t < timesForLineage.size(); t++) {
                
                posteriorProbs.get(lin).add(new ArrayList<Double>());
                double time = timesForLineage.get(t);
                int lineageStateArrayTimeIndex = lineageStateArrayTimes.indexOf(time);
                if (lineageStateArrayTimeIndex == -1) {
                    //Don't have states for this time
                    double probLineageState = 1.0 / states;
                    for (int st = 0; st < states; st++) {
                        posteriorProbs.get(lin).get(t).add(probLineageState);
                    }
                } else {
                    int lineageIndexAtTime = lineageStateArrayMap.get(lineageStateArrayTimeIndex).indexOf(lin);  
                    for (int st = 0; st < states; st++) {
                        double sampsInState = (double) array.get(lineageStateArrayTimeIndex).get(lineageIndexAtTime).get(st);
                        double probLineageState = sampsInState / samplesAsDouble;
                        posteriorProbs.get(lin).get(t).add(probLineageState);
                    }
                }
            }
        }
        
    }
    
    public void computeLineageStateEntropies(TreeNode[] tree, ZVectors dataZ)
    {   
        double sumEntropy; double stateProb; double entropyThisState;
        for (int lin = 0; lin < lineages; lin++) {
            
            int totalTimesForLineage = absoluteTimes.get(lin).size();
            entropies.add(new ArrayList<Double>());
            for (int t = 0; t < totalTimesForLineage; t++) {
                sumEntropy = 0.0;
                for (int state = 0; state < states; state++) {
                    stateProb = posteriorProbs.get(lin).get(t).get(state); 
                    entropyThisState =  -stateProb * (Math.log(stateProb)/Math.log(2.0)); //only for first particle
                    sumEntropy += entropyThisState;
                }
                entropies.get(lin).add(sumEntropy); 
            }
        }
        
    }
    
    public void discretizeEntropies(int classes)
    {
        //For all lineages at all times
        for (int lin = 0; lin < lineages; lin++) {
            discreteEntropies.add(new ArrayList<Integer>());
            int lineageTimes = absoluteTimes.get(lin).size();
            for (int time = 0; time < lineageTimes; time++) {
                int entropyClass = 0;
                double entropy = entropies.get(lin).get(time);
                if (entropy > 0.0) {
                    entropyClass = (int) Math.ceil(entropy * classes); //index from 1!    
                    discreteEntropies.get(lin).add(entropyClass);
                } else {
                    discreteEntropies.get(lin).add(1);
                }
            }
        }    
    }
    
    
}
