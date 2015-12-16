import cern.colt.list.*;
import java.util.ArrayList;
import cern.colt.matrix.*;
import java.util.Collections;

/**
 *
 * @author David
 */
public class ParticleLineageStateArray {
    
    //Should maybe have another array that holds the current lineage states for each particle
    //ArrayList<DoubleArrayList> absoluteTimes = new ArrayList<DoubleArrayList>();
    ArrayList<Double> absoluteTimes = new ArrayList<Double>();
    ArrayList<ArrayList<ArrayList<Integer>>> array;
    ArrayList<ArrayList<Integer>> arrayMap = new ArrayList<ArrayList<Integer>> ();
    ArrayList<ArrayList<Integer>> lineageStateCounts = new ArrayList<ArrayList<Integer>>();  
    ArrayList<ArrayList<Integer>> sampleArray;
    ArrayList<ArrayList<Integer>> trueArray;
    //ArrayList<ArrayList<Integer>> sampleArrayCopy;
    ArrayList<ArrayList<Integer>> sampleArrayLineageStates;
    ArrayList<ArrayList<Integer>> trueArrayLineageStates;
    ArrayList<Integer> activeLineages;
    int states;
    int lineages;
    int particles;
    int endParticle;
    int rootTimeIndex;
    ArrayList<Integer> lineagesAdded = new ArrayList<Integer>(); //just for debugging
    ArrayList<Integer> lineagesRemoved = new ArrayList<Integer>(); //just for debugging
    
    public void getTrueStatesFromTree(TreeNode[] tree, ZVectors dataZ) 
    {
        
        lineages = tree.length;
        int rootLineage = (lineages+1)/2;
        ArrayList<DoubleArrayList> allLineageAbsoluteTimes = new ArrayList<DoubleArrayList>();
        ArrayList<ArrayList<Integer>> allLineageStatesAtAbsoluteTimes = new ArrayList<ArrayList<Integer>>();
        String delimiter = ":";
        char COMMA = ',';
        int totalLineageStateChangeEvents = 0;

        for (int lineage = 0; lineage < lineages; ++lineage) {
            
            double nodeTime = tree[lineage].nodeHeight;
            String annotation = tree[lineage].nodeAnnotation;
            ArrayList<Integer> lineageIntervalStates = new ArrayList<Integer>(); //don't need to store these for all lineages
            ArrayList<Double> lineageIntervalTimes = new ArrayList<Double>(); //don't need to store these for all lineage
            lineageIntervalTimes.add(nodeTime);
            
            if (annotation != null) { //some nodes (like the root) may not have annotation
                
                String[] lineSegments = annotation.split(delimiter);
                //totalLineageStateChangeEvents += (lineSegments.length - 1);
                double cumulativeSegmentLengths = 0.0;
                for (int seg = 0; seg < lineSegments.length; seg++) {
                    int indexOfComma = lineSegments[seg].indexOf(COMMA, 0);
                    String segmentStateString = lineSegments[seg].substring(0,indexOfComma);
                    int segmentState = Integer.parseInt(segmentStateString) - 1; //minus one because these are indexed from one in the input file
                    lineageIntervalStates.add(segmentState);
                    String segmentLengthString = lineSegments[seg].substring(indexOfComma+1);
                    double segmentLength = Double.parseDouble(segmentLengthString);
                    if (segmentLength == 0.0) {
                        System.out.println("Segments found with length zero");
                    }
                    cumulativeSegmentLengths += segmentLength;
                    double segmentEndTime = nodeTime - cumulativeSegmentLengths;
                    lineageIntervalTimes.add(segmentEndTime);
                }
                
                double lineageStartTime = lineageIntervalTimes.get(lineSegments.length);
                double lineageEndTime = lineageIntervalTimes.get(0);
                int zLocStart = dataZ.absoluteTimes.indexOf(lineageStartTime);
                int zLocEnd = dataZ.absoluteTimes.indexOf(lineageEndTime);
                //if (zLocStart < 0) {
                    //System.out.println("Replaced lineage " + lineage + " start time");
                    //zLocStart = 1;
                //}
                DoubleArrayList lineAbsoluteTimes = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zLocStart, zLocEnd);
                ArrayList<Integer> lineStatesAtAbsTimes = new ArrayList<Integer>();
                int currSegment = 0;
                for(int time = lineAbsoluteTimes.size() - 1; time > 0; time--) {
                    double currTime = lineAbsoluteTimes.get(time);
                    if (currTime > lineageIntervalTimes.get(currSegment+1)) { //then stay in same segment
                        lineStatesAtAbsTimes.add(lineageIntervalStates.get(currSegment));
                    } else {
                        currSegment++;
                        lineStatesAtAbsTimes.add(lineageIntervalStates.get(currSegment));
                    }
                }
                lineStatesAtAbsTimes.add(lineageIntervalStates.get(currSegment)); //add state at final time
                Collections.reverse(lineStatesAtAbsTimes);
                allLineageAbsoluteTimes.add(lineAbsoluteTimes);
                allLineageStatesAtAbsoluteTimes.add(lineStatesAtAbsTimes);
 
            } else {
                System.out.println("Missing annotation for lineage: " + lineage);
                DoubleArrayList lineAbsoluteTimes = new DoubleArrayList();
                ArrayList<Integer> lineStatesAtAbsTimes = new ArrayList<Integer>();
                allLineageAbsoluteTimes.add(lineAbsoluteTimes);
                allLineageStatesAtAbsoluteTimes.add(lineStatesAtAbsTimes);
            }
        }
        //System.out.println("Total lineage state events observed in SimMap tree: " + totalLineageStateChangeEvents);
        
        sampleArray = new ArrayList<ArrayList<Integer>>();
        int totalTimes = arrayMap.size();
        int lineIndex; int lineAbsTimeIndex; int lineState;
        for (int k = 0; k < totalTimes; k++) {
            sampleArray.add(new ArrayList<Integer>());
            double timeNow = absoluteTimes.get(k);
            for (int lin = 0; lin < arrayMap.get(k).size(); lin++) { //for each lineage in tree at this time
                
                lineIndex = arrayMap.get(k).get(lin);
                lineAbsTimeIndex = allLineageAbsoluteTimes.get(lineIndex).indexOf(timeNow);
                lineState = -1;
                if (lineIndex == rootLineage) {
                    System.out.println("Missing a lineage state for lineage " + lineIndex);
                    if (lineIndex == rootLineage) {
                        
                        int daughter1IndexInMap = arrayMap.get(k).indexOf(tree[rootLineage].childNodes[0] - 1);
                        int daughter2IndexInMap = arrayMap.get(k).indexOf(tree[rootLineage].childNodes[1] - 1);
                        int daughter1State = sampleArray.get(k).get(daughter1IndexInMap);
                        int daughter2State = sampleArray.get(k).get(daughter2IndexInMap);
                        if (daughter1State == daughter2State) {
                            lineState = daughter1State; //sampleArray.get(k).add(daughter1State);
                        } else {
                            if (Math.random() < 0.5) {
                                lineState = daughter1State; //sampleArray.get(k).add(daughter1State);
                            } else {
                                lineState = daughter2State; //sampleArray.get(k).add(daughter2State);
                            }
                        }
                    }
                } else {
                    lineState = allLineageStatesAtAbsoluteTimes.get(lineIndex).get(lineAbsTimeIndex);
                }
                
                sampleArray.get(k).add(lineState);
                if (k > 0) {
                    int prevIndexInArray = arrayMap.get(k-1).indexOf(lineIndex);
                    if (prevIndexInArray != -1) {
                        int prevState = sampleArray.get(k-1).get(prevIndexInArray);
                        if (lineState != prevState) {
                            totalLineageStateChangeEvents++;
                        }
                    }
                }
            }
        }
        System.out.println("Total lineage state events observed in SimMap tree: " + totalLineageStateChangeEvents);
        this.getSampleArrayLineageStates();
        
//        sampleArrayCopy = new ArrayList<ArrayList<Integer>>();
//        for (int k = 0; k < totalTimes; k++) {
//            sampleArrayCopy.add(new ArrayList<Integer>());
//            for (int lin = 0; lin < arrayMap.get(k).size(); lin++) {
//                sampleArrayCopy.get(k).add(sampleArray.get(k).get(lin));
//            }
//        }      

    }
    
//    public void copyLineageStateSampleArray()
//    {
//        int totalTimes = arrayMap.size();
//        for (int k = 0; k < totalTimes; k++) {
//            for (int lin = 0; lin < arrayMap.get(k).size(); lin++) {
//                sampleArrayCopy.get(k).set(lin, sampleArray.get(k).get(lin));
//            }
//        }  
//    }
    
    //Defunct method
    public void getArrayOld(int jParticles, int lines, int sts) 
    {
        
        states = sts;
        particles = jParticles;
        lineages = lines;
        
        //First dimension is particles, second dimesnion is lineages, third dimension is times
        for (int j = 0; j < particles; j++) {
            array.add(new ArrayList<ArrayList<Integer>>());
            for (int k = 0; k < lineages; k++) {
                array.get(j).add(new ArrayList<Integer>());
            }
        }
        
        //Set up lineageStateCounts
        for (int j = 0; j < particles; j++) {
            lineageStateCounts.add(new ArrayList<Integer>());
            for (int k = 0; k < states; k++) {
                lineageStateCounts.get(j).add(0); //set initial number to zero
            }
        }
                
    }
    
    public void getArray(int jParticles, int totalTimes, int sts) 
    {
        particles = jParticles;
        states = sts;
        activeLineages = new ArrayList<Integer>();
        
        array = new ArrayList<ArrayList<ArrayList<Integer>>>();
        //First dimension is particles, second dimesnion is lineages, third dimension is times
        for (int j = 0; j < particles; j++) {
            array.add(new ArrayList<ArrayList<Integer>>());
            for (int k = 0; k < totalTimes; k++) {
                array.get(j).add(new ArrayList<Integer>());
            }
        }
        
        //Set up lineageStateCounts
        lineageStateCounts = new ArrayList<ArrayList<Integer>>();
        for (int j = 0; j < particles; j++) {
            lineageStateCounts.add(new ArrayList<Integer>());
            for (int k = 0; k < states; k++) {
                lineageStateCounts.get(j).add(0); //set initial number to zero
            }
        }
                
    }
    
    public void getConditionalArray() 
    {   
        array.get(particles-1).removeAll(array.get(particles-1));
        array.get(particles-1).addAll(sampleArray);
    }
    
    //Defunct method
    public void getAbsoluteTimes(DoubleArrayList xTimes, ZVectors dataZ, int zLocEnd, TreeNode[] tree) 
    {
//
//        for (int lin = 0; lin < lineages; lin++) {
//            absoluteTimes.add(new DoubleArrayList());
//        }
//        
//        ArrayList<Integer> currLineages = new ArrayList<Integer>();
//        int totalTimes = xTimes.size(); int zLocCurr; int backCounter = 0;
//        for (int time = (totalTimes-1); time > 0; time--) {
//            
//            zLocCurr = zLocEnd - backCounter;
//            double eventTime = dataZ.absoluteTimes.get(zLocCurr);
//            double nextTime = dataZ.absoluteTimes.get(zLocCurr - 1);
//            
//            //Add time for new samples
//            if (dataZ.omegaEvents.get(zLocCurr).contains(2)) { //if this interval begins with a sampling event
//                int totalEvents = dataZ.omegaEvents.get(zLocCurr).size();
//                for (int event = 0; event < totalEvents; event++) {
//                    if (dataZ.omegaEvents.get(zLocCurr).get(event) == 2) { //if sampling event
//                        int newLineage = dataZ.nodePointers.get(zLocCurr).get(event); //indexed from zero
//                        currLineages.add(newLineage);
//                        absoluteTimes.get(newLineage).add(eventTime);
//                    }
//                }
//            }
//            
//            //Add times for all currLineages
//            for (int line = 0; line < currLineages.size(); line++) {
//                int i = currLineages.get(line);
//                absoluteTimes.get(i).add(nextTime);
//            }
//            
//            //Add times for lineages entering through coalescent events
//            if (dataZ.omegaEvents.get(zLocCurr - 1).contains(1)) {
//                int totalEvents = dataZ.omegaEvents.get(zLocCurr - 1).size();
//                ArrayList<Integer> coalNodesAtEvent = new ArrayList<Integer>();
//                for (int event = 0; event < totalEvents; event++) {
//                    if (dataZ.omegaEvents.get(zLocCurr - 1).get(event)==1) {
//                        int coalNode = dataZ.nodePointers.get(zLocCurr - 1).get(event);
//                        coalNodesAtEvent.add(coalNode);
//                    }
//                }
//                int coalEvents = coalNodesAtEvent.size();
//                for (int event = 0; event < coalEvents; event++) {
//                    int coalNode = coalNodesAtEvent.get(event);
//                    int daughterLineage1Index = tree[coalNode].childNodes[0] - 1;
//                    int daughterLineage2Index = tree[coalNode].childNodes[1] - 1;
//                    absoluteTimes.get(coalNode).add(nextTime);
//                    //Remove daughter lineages from and add parent lineage to activeLineages
//                    int listIndex1 = currLineages.indexOf(daughterLineage1Index);
//                    if (listIndex1 != -1) {
//                        currLineages.remove(listIndex1);
//                    } else {
//                        System.out.println("Cannot remove daughter lineage");
//                    }
//
//                    int listIndex2 = currLineages.indexOf(daughterLineage2Index);
//                    if (listIndex2 != -1) {
//                        currLineages.remove(listIndex2);
//                    } else {
//                        System.out.println("Cannot remove daughter lineage");
//                    }
//                    currLineages.add(coalNode);
//                }
//            }
//            
//            backCounter++;
//        }
//        //System.out.println();
//        
   }
    
    public void getArrayMap(DoubleArrayList xTimes, ZVectors dataZ, int zLocEnd, TreeNode[] tree, StructCoalModel coal) 
    {
        //Set up lineage map for state array. Lineages in this array are ordered as they are in the array.
        lineages = tree.length;
        int rootNode = (tree.length + 1)/2;
        states = coal.Y.size();
        for (int k = 0; k < xTimes.size(); k++) {
            arrayMap.add(new ArrayList<Integer>());
        }
        ArrayList<Integer> currLineages = new ArrayList<Integer>();
        int totalTimes = xTimes.size(); int zLocCurr; int backCounter = 0;
        for (int time = (totalTimes-1); time > 0; time--) {
            
            zLocCurr = zLocEnd - backCounter;
            absoluteTimes.add(dataZ.absoluteTimes.get(zLocCurr));
            
            //Add new lineages to arrayMap
            if (dataZ.omegaEvents.get(zLocCurr).contains(2)) { //if this interval begins with a sampling event
                int totalEvents = dataZ.omegaEvents.get(zLocCurr).size();
                for (int event = 0; event < totalEvents; event++) {
                    if (dataZ.omegaEvents.get(zLocCurr).get(event) == 2) { //if sampling event
                        int newLineage = dataZ.nodePointers.get(zLocCurr).get(event); //indexed from zero
                        currLineages.add(newLineage);
                        arrayMap.get(time).add(newLineage);
                    }
                }
            }
            
            //Update lineages in arrayMap for the next time step
            for (int line = 0; line < currLineages.size(); line++) {
                int i = currLineages.get(line);
                arrayMap.get(time-1).add(i); //lineages are added to the arrayMap in the same order as they will be in activeLineages
            }
            
            //Add lineages to the arrayMap entering at the coalescent event
            if (dataZ.omegaEvents.get(zLocCurr - 1).contains(1)) {
                int totalEvents = dataZ.omegaEvents.get(zLocCurr - 1).size();
                ArrayList<Integer> coalNodesAtEvent = new ArrayList<Integer>();
                for (int event = 0; event < totalEvents; event++) {
                    if (dataZ.omegaEvents.get(zLocCurr - 1).get(event)==1) {
                        int coalNode = dataZ.nodePointers.get(zLocCurr - 1).get(event);
                        coalNodesAtEvent.add(coalNode);
                    }
                }
                int coalEvents = coalNodesAtEvent.size();
                for (int event = 0; event < coalEvents; event++) {
                    int coalNode = coalNodesAtEvent.get(event);
                    int daughterLineage1Index = tree[coalNode].childNodes[0] - 1;
                    int daughterLineage2Index = tree[coalNode].childNodes[1] - 1;
                    arrayMap.get(time-1).add(coalNode);
                    
                    //Remove daughter lineages from and add parent lineage to activeLineages
                    int listIndex1 = currLineages.indexOf(daughterLineage1Index);
                    if (listIndex1 != -1) {
                        currLineages.remove(listIndex1);
                    } else {
                        System.out.println("Cannot remove daughter lineage");
                    }

                    int listIndex2 = currLineages.indexOf(daughterLineage2Index);
                    if (listIndex2 != -1) {
                        currLineages.remove(listIndex2);
                    } else {
                        System.out.println("Cannot remove daughter lineage");
                    }
                    
                    if (coalNode != rootNode) { //Don't add root lineage if hit root
                        currLineages.add(coalNode);
                    } else {
                        rootTimeIndex = time-1;
                    }
                    
                }
            }
            
            backCounter++;
        }
        absoluteTimes.add(dataZ.absoluteTimes.get(zLocEnd - backCounter));
        Collections.reverse(absoluteTimes);
        //System.out.println();
        
    }
    
    //Defunct method
    public void addSamplesOld(ZVectors dataZ, int zLoc, TreeNode[] tree) 
    {
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        //double eventTime = dataZ.absoluteTimes.get(zLoc);
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                int newLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                activeLineages.add(newLineage);
                lineagesAdded.add(newLineage);
                //absoluteTimes.get(newLineage).add(eventTime);
                int newState = tree[newLineage].nodeLabel - 1; //nodeLabels are indexed from one
                for (int j = 0; j < particles; j++) {
                    if (newState == -1) {
                        System.out.println("Hit unknown lineage states");
                    }
                    array.get(j).get(newLineage).add(newState);
                    lineageStateCounts.get(j).set(newState, (lineageStateCounts.get(j).get(newState) + 1));
                }
            }
        }    
    }
    
    public void addSamples(ZVectors dataZ, int zLoc, TreeNode[] tree, int time) 
    {
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        //double eventTime = dataZ.absoluteTimes.get(zLoc);
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                int newLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                activeLineages.add(newLineage);
                //lineagesAdded.add(newLineage);
                int newState = tree[newLineage].nodeLabel - 1; //nodeLabels are indexed from one
                for (int j = 0; j < (particles-1); j++) { //CONDITIONING ON PREVIOUS Lg SO NOT ADDING NEW SAMPLES TO LAST PARTICLE
                    if (newState == -1) {
                        System.out.println("Hit unknown lineage states");
                    }
                    //array.get(j).get(newLineage).add(newState);
                    array.get(j).get(time).add(newState);
                    //arrayMap.get(time).add(newLineage);
                    lineageStateCounts.get(j).set(newState, (lineageStateCounts.get(j).get(newState) + 1));
                }
                //UPDATE LINEAGE STATE COUNTS FOR THE LAST PARTICLE
                lineageStateCounts.get(particles-1).set(newState, (lineageStateCounts.get(particles-1).get(newState) + 1));
            }
        }    
    }
    
    //Defunct method
    public void updateStatesOld(DoubleMatrix1D xCurr, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime) 
    {
//        
//        int numLineages = activeLineages.size(); int timeIndex;
//        int i; int currState; int nextState; double birthRateOut; double migrationRateOut; int linesInLNotInA;
//        DoubleFactory2D factory2D;
//	factory2D = DoubleFactory2D.dense;
//	DoubleMatrix2D transMatrix = factory2D.make(states,states);
//        DoubleMatrix2D normCumTransMatrix = factory2D.make(states,states);
//        ArrayList<Integer> lineagesInAByState;
//        
//        //Can update this outside the loop over particles since only have one xTraj
//        coal.updateF(xCurr, theta, currTime);
//        coal.updateG(xCurr, theta);
//        coal.updateY(xCurr);
//        
//        //For each particle
//        for (int x = 0; x < particles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j
//            
//            //Precompute A(l) terms by summing probs in matrix
//            //ArrayList<Integer> lineagesInAByState = this.getLineagesByState(x, currTime);
//            lineagesInAByState = lineageStateCounts.get(x);
//            
//            //Can precompute lineage transition probabilties here
//            double rowSum; double pMove;
//            for (int k = 0; k < states; k++) {
//                rowSum = 0.0;
//                for (int l = 0; l < states; l++) {
//                    if (l != k) {
//                        if (coal.Y.getQuick(k) > 0.0 & coal.Y.getQuick(l) > 0.0) {
//                            linesInLNotInA = (int) coal.Y.getQuick(l) - lineagesInAByState.get(l);
//                            if (linesInLNotInA > 0) {
//                                birthRateOut = (coal.F.getQuick(l,k)/coal.Y.getQuick(k))*((linesInLNotInA) / coal.Y.getQuick(l));
//                            } else {
//                                birthRateOut = 0.0;
//                            }
//                        } else {
//                            birthRateOut = 0.0;
//                        }
//                        if (coal.Y.getQuick(k) > 0.0) {
//                            migrationRateOut = coal.G.getQuick(l,k)/coal.Y.getQuick(k);
//                        } else {
//                            migrationRateOut = 0.0;
//                        }
//                        pMove = 1 - Math.exp(-((birthRateOut + migrationRateOut)*dtTime));
//                        rowSum += pMove;
//                        transMatrix.setQuick(k, l, pMove);
//                    }   
//                }
//                transMatrix.setQuick(k, k, (1-rowSum));
//            }
//            
//            //Compute normalized cumulative transition probabilities
//            for (int k = 0; k < states; k++) {
//                double cumSum = 0;
//                for (int l = 0; l < states; l++) {
//                    cumSum += transMatrix.getQuick(k,l);
//                    normCumTransMatrix.setQuick(k,l,cumSum);
//                }
//            }
//            //System.out.println();
//            
//            for (int lin = 0; lin < numLineages; lin++) { //for all lineages i
//                
//                i = activeLineages.get(lin);
//                timeIndex = absoluteTimes.get(i).indexOf(currTime);
//                currState = array.get(x).get(i).get(timeIndex);
//                nextState = this.getNextState(currState, normCumTransMatrix);
//                if (currState != nextState) { //update lineageStateCounts
//                    lineageStateCounts.get(x).set(currState, (lineageStateCounts.get(x).get(currState) - 1));
//                    lineageStateCounts.get(x).set(nextState, (lineageStateCounts.get(x).get(nextState) + 1));
//                }
//                array.get(x).get(i).add(nextState);
//                
//            }
//        }
    }
    
    public void updateStates(DoubleMatrix1D xCurr, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime, ArrayList<Integer> events) 
    {
        
        int numLineages = activeLineages.size(); //int timeIndex;
        int currState; int nextState; double birthRateOut; double migrationRateOut; int linesInLNotInA;
        double uniRandNextState; double cumValNextState; int counterNextState; double cumSum;
        DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
	DoubleMatrix2D transMatrix = factory2D.make(states,states);
        DoubleMatrix2D normCumTransMatrix = factory2D.make(states,states);
        DoubleMatrix2D birthRateMatrix = factory2D.make(states,states);
        DoubleMatrix2D migrationRateMatrix = factory2D.make(states,states);
        ArrayList<Integer> lineagesInAByState;
        
        //Can update this outside the for loop over particles since only have one xTraj
        coal.updateF(xCurr, theta, currTime);
        coal.updateG(xCurr, theta);
        coal.updateY(xCurr);
        
        //Precompute some constant values used to compute the transition probabilities
        for (int k = 0; k < states; k++) {
            for (int l = 0; l < states; l++) {
                if (l != k) {
                birthRateOut = (coal.F.getQuick(l,k)/coal.Y.getQuick(k))*(1.0 / coal.Y.getQuick(l));
                birthRateMatrix.setQuick(k,l,birthRateOut);
                migrationRateOut = coal.G.getQuick(l,k)/coal.Y.getQuick(k);
                migrationRateMatrix.setQuick(k,l,migrationRateOut);
                }
            }
        }
        
        //Get previous lineage indexes for current activeLineages
        ArrayList<Integer> prevLineageIndexes = new ArrayList<Integer>();
        for (int lin = 0; lin < numLineages; lin++) {
            int lineageIndex = arrayMap.get(time-1).get(lin);
            prevLineageIndexes.add(arrayMap.get(time).indexOf(lineageIndex));    
        }
        
        //For each particle
        for (int x = 0; x < (particles-1); x++) { //conditioning on previous Lg so do not update the last particle
            
            //Precompute A(l) terms by summing probs in matrix
            lineagesInAByState = lineageStateCounts.get(x);
            
            //Can precompute lineage transition probabilties here
            double rowSum; double pMove;
            for (int k = 0; k < states; k++) {
                rowSum = 0.0;
                for (int l = 0; l < states; l++) {
                    if (l != k) {
                        if (coal.Y.getQuick(k) > 0.0 & coal.Y.getQuick(l) > 0.0) {
                            linesInLNotInA = (int) coal.Y.getQuick(l) - lineagesInAByState.get(l);
                            if (linesInLNotInA > 0) {
                                birthRateOut = birthRateMatrix.getQuick(k,l) * linesInLNotInA;
                            } else {
                                birthRateOut = 0.0;
                            }
                        } else {
                            birthRateOut = 0.0;
                        }
                        if (coal.Y.getQuick(k) > 0.0) {
                            migrationRateOut = migrationRateMatrix.getQuick(k,l);
                        } else {
                            migrationRateOut = 0.0;
                        }
                        pMove = 1 - Math.exp(-((birthRateOut + migrationRateOut)*dtTime));
                        rowSum += pMove;
                        transMatrix.setQuick(k, l, pMove);
                    }   
                }
                transMatrix.setQuick(k, k, (1-rowSum));
            }
            
            //Compute normalized cumulative transition probabilities
            for (int k = 0; k < states; k++) {
                cumSum = 0;
                for (int l = 0; l < states; l++) {
                    cumSum += transMatrix.getQuick(k,l);
                    normCumTransMatrix.setQuick(k,l,cumSum);
                }
            }
            
            for (int ln = 0; ln < numLineages; ln++) { //for all lineages i
                
                currState = array.get(x).get(time).get(prevLineageIndexes.get(ln));
                uniRandNextState = Math.random();
                cumValNextState = 0.0;
                counterNextState = 0;
                while (cumValNextState <= uniRandNextState)
                {
                    cumValNextState = normCumTransMatrix.getQuick(currState, counterNextState);
                    counterNextState++;
                }   
                nextState = counterNextState - 1; 
                array.get(x).get(time-1).add(nextState);
                if (currState != nextState) { //update lineageStateCounts
                    lineageStateCounts.get(x).set(currState, (lineageStateCounts.get(x).get(currState) - 1));
                    lineageStateCounts.get(x).set(nextState, (lineageStateCounts.get(x).get(nextState) + 1));
                }
                
            }
        }
        //UPDATE LINEAGE STATE COUNTS FOR LAST PARTICLE
        for (int ln = 0; ln < numLineages; ln++) { //for all lineages i
            
            int lineageIndex = arrayMap.get(time).get(prevLineageIndexes.get(ln));
            currState = array.get(particles-1).get(time).get(prevLineageIndexes.get(ln));
            int lineageNextIndex = arrayMap.get(time-1).indexOf(lineageIndex);
            nextState = array.get(particles-1).get(time-1).get(lineageNextIndex); 
            if (currState != nextState) { //update lineageStateCounts
                lineageStateCounts.get(particles-1).set(currState, (lineageStateCounts.get(particles-1).get(currState) - 1));
                lineageStateCounts.get(particles-1).set(nextState, (lineageStateCounts.get(particles-1).get(nextState) + 1));
            }
        }
    }
    
    public void updateStatesAfterCoal(int x, StructCoalModel coal, int coalNode, int leftDaughterState, int rightDaughterState, int time) 
    {
        
        //Probabilistically choose a state for the parent of the two daughter lineages
        int coalNodeState = -1;
        
        if (leftDaughterState != rightDaughterState) {
            double fKtoL = coal.F.get(leftDaughterState, rightDaughterState);
            double fLtoK = coal.F.get(rightDaughterState, leftDaughterState);
            double pKtoL = fKtoL / (fKtoL + fLtoK);  //probability left daughter transmitted
            double uniRand = Math.random();
            //This only works if there's just two states
            if (uniRand <= pKtoL) {   
                coalNodeState = leftDaughterState;   
            } else {
            coalNodeState = rightDaughterState;    
            }
        } else {
            coalNodeState = leftDaughterState;
        }
        //array.get(x).get(coalNode).add(coalNodeState);
        array.get(x).get(time-1).add(coalNodeState);
        //Update lineageStateCounts
        lineageStateCounts.get(x).set(leftDaughterState, (lineageStateCounts.get(x).get(leftDaughterState) - 1));
        lineageStateCounts.get(x).set(rightDaughterState, (lineageStateCounts.get(x).get(rightDaughterState) - 1));
        lineageStateCounts.get(x).set(coalNodeState, (lineageStateCounts.get(x).get(coalNodeState) + 1));
 
    }
    
    public void updateStatesAfterCoalForFixedStates(int x, int coalNode, int leftDaughterState, int rightDaughterState, int time) 
    {

        int coalNodeIndex = arrayMap.get(time-1).indexOf(coalNode);
        int coalNodeState = array.get(x).get(time-1).get(coalNodeIndex);
        //Update lineageStateCounts
        lineageStateCounts.get(x).set(leftDaughterState, (lineageStateCounts.get(x).get(leftDaughterState) - 1));
        lineageStateCounts.get(x).set(rightDaughterState, (lineageStateCounts.get(x).get(rightDaughterState) - 1));
        lineageStateCounts.get(x).set(coalNodeState, (lineageStateCounts.get(x).get(coalNodeState) + 1));
 
    }
    
    public void getSampleArray(int particle) 
    {
        //Select one particle to be the sample
        sampleArray = new ArrayList<ArrayList<Integer>>();
        sampleArray = array.get(particle); //this might not be a deep copy
        
    }
    
    public void getNewSampleArray() 
    {
        //This is for block sampling
        sampleArray = new ArrayList<ArrayList<Integer>>();
        for (int t = 0; t < absoluteTimes.size(); t++) {
            sampleArray.add(new ArrayList<Integer>());
        }    
    }
    
    public void updateStateArrayForNextBlock(int time) 
    {
        
        ArrayList<Integer> prevBlockSampleStates = sampleArray.get(time);
        ArrayList<Integer> prevBlockLineageStateCounts = lineageStateCounts.get(endParticle);
        for (int x = 0; x < (particles-1); x++) {
            array.get(x).get(time).removeAll(array.get(x).get(time));
            array.get(x).get(time).addAll(prevBlockSampleStates);
            lineageStateCounts.get(x).removeAll(lineageStateCounts.get(x));
            lineageStateCounts.get(x).addAll(prevBlockLineageStateCounts);
        }
        System.out.println();
        
    }
    
    //Defunct method
    public void getSampleArrayLineageStatesOld(DoubleArrayList xTimes) 
    {
//        //Get number of lineages in each state at each time
//        int currState; double currTime; int linesInState; int lineTimeIndex;
//        sampleArrayLineageStates = new ArrayList<ArrayList<Integer>>();
//        for (int t = 0; t < xTimes.size(); t++) {
//            currTime = xTimes.get(t);
//            sampleArrayLineageStates.add(new ArrayList<Integer>());
//            for (int s = 0; s < states; s++) {
//                linesInState = 0;
//                for (int lin = 0; lin < lineages; lin++) {
//                    lineTimeIndex = absoluteTimes.get(lin).indexOf(currTime);
//                    if (lineTimeIndex > -1) {
//                        currState = sampleArray.get(lin).get(lineTimeIndex);
//                        if (currState == s) {
//                            linesInState++;
//                        }
//                    } 
//                }
//                sampleArrayLineageStates.get(t).add(linesInState);
//            }
//        }   
    }
    
    public void getSampleArrayLineageStates() 
    {
        //Get number of lineages in each state at each time
        int currState;
        sampleArrayLineageStates = new ArrayList<ArrayList<Integer>>();
        for (int t = 0; t < sampleArray.size(); t++) {
            sampleArrayLineageStates.add(new ArrayList<Integer>());
            for (int s = 0; s < states; s++) {
                sampleArrayLineageStates.get(t).add(0);
            }
            for (int lin = 0; lin < sampleArray.get(t).size(); lin++) {
                currState = sampleArray.get(t).get(lin);
                sampleArrayLineageStates.get(t).set(currState, (sampleArrayLineageStates.get(t).get(currState) + 1));
            }
        }   
    }
    
    
    
    public DoubleArrayList getLineageStatesAsDoubleArray(int state)
    {
        DoubleArrayList stateList = new DoubleArrayList();
        int currState;
        for (int t = 0; t < sampleArrayLineageStates.size(); t++) {
            currState = sampleArrayLineageStates.get(t).get(state);
            stateList.add(currState);
        }
        return stateList;
    }
    
    public void setTrueStates() 
    {
        trueArray = new ArrayList<ArrayList<Integer>>();
        //Copy states to trueArray
        for (int t = 0; t < sampleArray.size(); t++) {
            trueArray.add(new ArrayList<Integer>());
            for (int lin = 0; lin < sampleArray.get(t).size(); lin++) {
                trueArray.get(t).add(sampleArray.get(t).get(lin));
            }  
        }
        
        trueArrayLineageStates = new ArrayList<ArrayList<Integer>>();
        for (int t = 0; t < sampleArray.size(); t++) {
            trueArrayLineageStates.add(new ArrayList<Integer>());
            for (int s = 0; s < states; s++) {
                trueArrayLineageStates.get(t).add(0);
            }
            for (int lin = 0; lin < sampleArray.get(t).size(); lin++) {
                int currState = sampleArray.get(t).get(lin);
                trueArrayLineageStates.get(t).set(currState, (sampleArrayLineageStates.get(t).get(currState) + 1));
            }
        } 
        
    }
    
    
    private int getNextState(int currState, DoubleMatrix2D normCumTransMatrix) 
    {
        //Probabilistically choose new state based on the current transition probabilities
        
        double uniRand = Math.random();
        double cumVal = 0.0;
        int counter = 0;
        while (cumVal <= uniRand)
        {
            cumVal = normCumTransMatrix.getQuick(currState, counter);
            counter++;
        }   
        int nextState = counter - 1;
        return nextState; 
    }
    
    //Defunct method
//    public ArrayList<Integer> getLineagesByState(int x, double currTime) 
//    { 
//        //Could avoid doing this if had extra array that held the number of lineages in each state for each particle
//        
//        ArrayList<Integer> lineagesByState = new ArrayList<Integer>();
//        int numLineages = activeLineages.size(); int timeIndex;
//        int i; int stateSum;
//        for (int k = 0; k < states; k++) {
//            stateSum = 0;
//            for (int lin = 0; lin < numLineages; lin++) {
//                i = activeLineages.get(lin);
//                timeIndex = absoluteTimes.get(i).indexOf(currTime);
//                if (array.get(x).get(i).get(timeIndex) == k) {
//                    stateSum++;
//                }
//            }
//            lineagesByState.add(stateSum);
//        }
//        return lineagesByState;
//    }
    
}
