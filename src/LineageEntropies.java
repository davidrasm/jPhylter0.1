import cern.colt.list.*;
import java.util.ArrayList;

/**
 *
 * @author David
 */
public class LineageEntropies {
    
    ArrayList<DoubleArrayList> absoluteTimes = new ArrayList<DoubleArrayList>();
    ArrayList<ArrayList<Double>> entropies = new ArrayList<ArrayList<Double>>();
    ArrayList<ArrayList<Integer>> discreteEntropies = new ArrayList<ArrayList<Integer>>();
    ArrayList<Integer> activeLineages = new ArrayList<Integer>();
    int states;
    int lineages;
    
    
    public void getArrays(int lns, int sts) 
    {
        
        states = sts;
        lineages = lns;    
        for (int lin = 0; lin < lineages; lin++) {
            absoluteTimes.add(new DoubleArrayList());
            entropies.add(new ArrayList<Double>());    
        }
        
    }
    
    public void addSamples(ZVectors dataZ, int zLoc, TreeNode[] tree) 
    {
       
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        double eventTime = dataZ.absoluteTimes.get(zLoc);
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                int newLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                activeLineages.add(newLineage);
                absoluteTimes.get(newLineage).add(eventTime);
                
                //For entropies
                entropies.get(newLineage).add(0.0);
                
                //For state probabilities
//                int newState = tree[newLineage].nodeLabel - 1;
//                if (newState == 0) {
//                    entropies.get(newLineage).add(1.0);
//                } else {
//                    entropies.get(newLineage).add(0.0);
//                }
            }
        }
        
    }
    
    public void updateEntropies(double time, LineageStateProbs stateProbs) 
    {
        for (int lin = 0; lin < activeLineages.size(); lin++) {
            int i = activeLineages.get(lin);
            absoluteTimes.get(i).add(time);
            
            //For entropies
            double sumEntropy = 0.0;
            for (int state = 0; state < states; state++) {
                double stateProb = stateProbs.matrix.get(0, i, state); 
                double entropyThisState =  -stateProb * (Math.log(stateProb)/Math.log(2.0)); //only for first particle
                sumEntropy += entropyThisState;
            }
            entropies.get(i).add(sumEntropy); 
            
            //For state probabilities
            //entropies.get(i).add(stateProbs.matrix.get(0,i,0));
            
        }
        
    }
    
    public void updateEntropiesCoal(ZVectors dataZ, int zLoc, TreeNode[] tree, LineageStateProbs stateProbs) 
    {
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        double eventTime = dataZ.absoluteTimes.get(zLoc);
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 1) { //if coalescent event
                
                int coalNode = dataZ.nodePointers.get(zLoc).get(event);
                int daughterLineage1Index = tree[coalNode].childNodes[0] - 1;
                int daughterLineage2Index = tree[coalNode].childNodes[1] - 1;
                
                int listIndex1 = activeLineages.indexOf(daughterLineage1Index);
                if (listIndex1 != -1) {
                    activeLineages.remove(listIndex1);
                }
                int listIndex2 = activeLineages.indexOf(daughterLineage2Index);
                if (listIndex2 != -1) {
                    activeLineages.remove(listIndex2);
                }
                
                activeLineages.add(coalNode);
                absoluteTimes.get(coalNode).add(eventTime);
                
                //For entropies
                double sumEntropy = 0.0;
                for (int state = 0; state < states; state++) {
                    double stateProb = stateProbs.matrix.get(0, coalNode, state); 
                    double entropyThisState =  stateProb * (Math.log(stateProb)/Math.log(2.0)); //only for first particle
                    sumEntropy += entropyThisState;
                }
                entropies.get(coalNode).add(sumEntropy);
                
                //For state probabilities
                //entropies.get(coalNode).add(stateProbs.matrix.get(0, coalNode, 0));
                
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