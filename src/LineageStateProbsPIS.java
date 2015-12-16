import cern.colt.list.*;
import cern.colt.matrix.*;
import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author David
 */
public class LineageStateProbsPIS {
    
    int particles;
    int states;
    int lineages;
    int rootNode;
    ArrayList<ArrayList<ArrayList<Double>>> probs; //times, lineages, states
    ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> particleProbs;
    ArrayList<ArrayList<ArrayList<Double>>> particleMiniProbs;
    ArrayList<Integer> activeLineages = new ArrayList<Integer>();
    //ArrayList<ArrayList<Integer>> activeLineageMap = new ArrayList<ArrayList<Integer>>();
    //ArrayList<ArrayList<Integer>> activeLineageIndexes = new ArrayList<ArrayList<Integer>>();
    ArrayList<ArrayList<Integer>> lineageMap = new ArrayList<ArrayList<Integer>> ();
    ArrayList<Double> absoluteTimes = new ArrayList<Double>();
    DoubleArrayList stateSums = new DoubleArrayList();
    
    //ArrayList<Integer> lineagesRemoved = new ArrayList<Integer>();
    //ArrayList<Integer> lineagesAdded = new ArrayList<Integer>();
    
    DoubleMatrix2D transProbs;
    
    /**
    public void getMatrix(int jParticles, int lines, int sts, DoubleFactory3D factory3D) 
    {
        
        states = sts;
        particles = jParticles;
        lineages = lines;
        matrix = factory3D.make(particles, lineages, sts);
        
    } */
    
    public void setUpProbArray() 
    {
        
        //Set up array for lineage state probs
        probs = new ArrayList<ArrayList<ArrayList<Double>>>();
        int times = lineageMap.size();
        int lineagesAtTime = 0;
        for (int t = 0; t < times; t++) {
            probs.add(new ArrayList<ArrayList<Double>>());
            lineagesAtTime = lineageMap.get(t).size();
            for (int lin = 0; lin < lineagesAtTime; lin++) {
                probs.get(t).add(new ArrayList<Double>());
                for (int s = 0; s < states; s++) {
                    probs.get(t).get(lin).add(0.0);
                }
            }
        }

        for (int st = 0; st < states; st++) {
            stateSums.add(0.0); 
        }
        
        //Set up tranProb matrix
        DoubleFactory2D factory2D;
        factory2D = DoubleFactory2D.dense;
        transProbs = factory2D.make(states, states);
        
    }
    
    public void getLineageMap(DoubleArrayList xTimes, ZVectors dataZ, int zLocEnd, TreeNode[] tree, StructCoalModel coal) 
    {
        //Set up lineage map for state array. Lineages in this array are ordered as they are in the array.
        lineages = tree.length;
        rootNode = (tree.length + 1)/2;
        states = coal.Y.size();
        for (int k = 0; k < xTimes.size(); k++) {
            lineageMap.add(new ArrayList<Integer>());
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
                        lineageMap.get(time).add(newLineage);
                    }
                }
            }
            
            //Update lineages in arrayMap for the next time step
            for (int line = 0; line < currLineages.size(); line++) {
                int i = currLineages.get(line);
                lineageMap.get(time-1).add(i); //lineages are added to the arrayMap in the same order as they will be in activeLineages
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
                    lineageMap.get(time-1).add(coalNode);
                    
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
                        //rootTimeIndex = time-1;
                    }
                    
                }
            }
            
            backCounter++;
        }
        absoluteTimes.add(dataZ.absoluteTimes.get(zLocEnd - backCounter));
        Collections.reverse(absoluteTimes);
        //System.out.println();
        
    }
    
    /***Defunct method***
    public void updateActiveLineageArray() 
    {
        activeLineageMap.add(new ArrayList<Integer>());
        activeLineageMap.get(activeLineageMap.size() - 1).addAll(activeLineages);
        activeLineageIndexes.add(new ArrayList<Integer>());
        for (int i = 0; i < activeLineages.size(); i++) {
            activeLineageIndexes.get(activeLineageMap.size() - 1).add(i);
        }
    } */
    
    public void addSamples(ZVectors dataZ, int zLoc, TreeNode[] tree, int time) 
    {
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                int newLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                int newLineageIndex = lineageMap.get(time).indexOf(newLineage);
                activeLineages.add(newLineage);
                int newState = tree[newLineage].nodeLabel - 1; //nodeLabels are indexed from one
                if (newState == -1) {
                    System.out.println("Hit unknown lineage states");
                }
                probs.get(time).get(newLineageIndex).set(newState, 1.0);
            }
        }    
    }
    
    public void addSamplesWithPriors(ZVectors dataZ, int zLoc, TreeNode[] tree, int time) 
    {
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                int newLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                int newLineageIndex = lineageMap.get(time).indexOf(newLineage);
                activeLineages.add(newLineage);
                for (int st = 0; st < states; st++) {
                    probs.get(time).get(newLineageIndex).set(st, tree[newLineage].statePriors.get(st));
                }
            }
        }
        
    }
    
    public void updateProbs(DoubleMatrix2D xCurr, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime) 
    {
        
        int numLineages = activeLineages.size(); //int timeIndex;
        
        //Can probably update this outside of this method since this will be the same for all
        int pIndex = 0;
        coal.updateF(xCurr.viewColumn(pIndex), theta, currTime);
        coal.updateG(xCurr.viewColumn(pIndex), theta);
        coal.updateY(xCurr.viewColumn(pIndex));
        
        //Get previous lineage indexes for current activeLineages
        ArrayList<Integer> prevLineageIndexes = new ArrayList<Integer>();
        ArrayList<Integer> currLineageIndexes = new ArrayList<Integer>();
        for (int lin = 0; lin < numLineages; lin++) {
            int currLine = activeLineages.get(lin);
            prevLineageIndexes.add(lineageMap.get(time).indexOf(currLine));
            currLineageIndexes.add(lineageMap.get(time-1).indexOf(currLine));   
        }
        
        //Sum over S's to get A's here - this could be done automatically recomputed at the end of each time step
        this.updateStateSums(time);
        
        //Precompute transition probabilties between all states and store these for forwards pass
        double probOut; //prob of transitioning from k -> l backwards in time
        double rateOutByBirth;
        double rateOutByMigration;
        double rowSum; //double rowSumMethod0;
        for (int k = 0; k < states; ++k) { //For each state or subpopulation
            rowSum = 0.0;
            for (int l = 0; l < states; ++l) { //<loop through l here so you can reuse the following four methods again on the up pass
                if (k != l) {
                    
                    //New way for dealing with Y's = 0
                    if (coal.Y.getQuick(k) == 0.0 || coal.Y.getQuick(l) == 0.0) {
                        rateOutByBirth = 0.0;
                    } else {
                        if (coal.Y.getQuick(l) > stateSums.get(l)) {
                            rateOutByBirth = (coal.F.getQuick(l,k)/coal.Y.get(k)) * ((coal.Y.getQuick(l)-stateSums.get(l))/coal.Y.getQuick(l));
                        } else {
                            rateOutByBirth = 0.0;
                        }
                    }
                    if (coal.Y.getQuick(k) == 0.0) {
                        rateOutByMigration = 0.0;
                    } else {
                        rateOutByMigration = (coal.G.getQuick(l,k)/coal.Y.getQuick(k));
                    }

                    probOut = 1 - Math.exp(-(rateOutByBirth + rateOutByMigration)*dtTime);
                    transProbs.setQuick(k, l, probOut);
                    rowSum += probOut;
                    
                }
            }
            //transProbs.setQuick(k, k, (1-rowSum)); //Necessary??
            
        }
        
        //Compute new cumulative probs for each lineage
        //Update this on 03/11/13 so that indexing should always reflect indexes in arrayMap
        double probMassIntoL; double newSk; double newSl;
        int lineIndex; int prevLineageIndex;
        
        for (int lin = 0; lin < numLineages; lin++) {
            
            lineIndex = currLineageIndexes.get(lin);
            prevLineageIndex = prevLineageIndexes.get(lin);
            for (int st = 0; st < states; st++) {
                probs.get(time-1).get(lineIndex).set(st, probs.get(time).get(prevLineageIndex).get(st));
            }
            
            for (int k = 0; k < states; ++k) {
                for (int l = 0; l < states; ++l) {
                    if (k != l) {
                        probMassIntoL = transProbs.getQuick(k, l) * probs.get(time).get(prevLineageIndex).get(k);
                        newSk = probs.get(time-1).get(lineIndex).get(k) - probMassIntoL;
                        probs.get(time-1).get(lineIndex).set(k, newSk);
                        newSl = probs.get(time-1).get(lineIndex).get(l) + probMassIntoL;
                        probs.get(time-1).get(lineIndex).set(l, newSl);
                    }
                }
            }
            
        }
        
    }
    
    public void updateProbsCoalEvent(ZVectors dataZ, int zLoc, TreeNode[] tree, DoubleMatrix2D xCurr, DoubleArrayList theta, int time, StructCoalModel coal, double currTime) 
    {
        
        int pIndex = 0;
        coal.updateF(xCurr.viewColumn(pIndex), theta, currTime); //xCurr from time-1
        coal.updateG(xCurr.viewColumn(pIndex), theta);
        coal.updateY(xCurr.viewColumn(pIndex));
        
        DoubleFactory2D factory2D = DoubleFactory2D.dense;
        DoubleMatrix2D coalRates = factory2D.make(states, states);
        //DoubleMatrix2D otherCoalRates = factory2D.make(states, states);
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        int mapIndexDaughter1; int mapIndexDaughter2; int mapIndexParent;
        
        //Find coalescent events at this time
        ArrayList<Integer> coalNodesAtEvent = new ArrayList<Integer>();
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event)==1) {
                int coalNode = dataZ.nodePointers.get(zLoc).get(event);
                coalNodesAtEvent.add(coalNode);
            }
        }
        
        int coalEvents = coalNodesAtEvent.size();
        for (int event = 0; event < coalEvents; event++) {
                
            int coalNode = coalNodesAtEvent.get(event);
            int daughterLineage1Index = tree[coalNode].childNodes[0] - 1;
            int daughterLineage2Index = tree[coalNode].childNodes[1] - 1;
    
            mapIndexDaughter1 = lineageMap.get(time-1).indexOf(daughterLineage1Index);
            mapIndexDaughter2 = lineageMap.get(time-1).indexOf(daughterLineage2Index);
            mapIndexParent = lineageMap.get(time-1).indexOf(coalNode);
            
            ArrayList<Double> daughter1Probs = probs.get(time-1).get(mapIndexDaughter1);
            ArrayList<Double> daughter2Probs = probs.get(time-1).get(mapIndexDaughter2);
            
            double lambda; double lambdaSum = 0.0;
            for (int k = 0; k < states; k++) {
                for (int l = 0; l < states; l++) {
                    
                    //New way for dealing with Y's = 0
                    if (coal.Y.getQuick(k) == 0.0 || coal.Y.getQuick(l) == 0.0) {
                        lambda = 0.0;
                    } else {
                        //THERE WAS A BUG HERE WITH HOW LAMBDA WAS BEING COMPUTED WITH THE INDEXING OF K AND L IN DAUGTHER PROBS
                        //lambda = (coal.F.getQuick(k,l)/(coal.Y.getQuick(k)*coal.Y.getQuick(l))) * (daughter1Probs.get(k)*daughter2Probs.get(l) + daughter2Probs.get(l)*daughter1Probs.get(k));
                        lambda = (coal.F.getQuick(k,l)/(coal.Y.getQuick(k)*coal.Y.getQuick(l))) * ((daughter1Probs.get(k)*daughter2Probs.get(l)) + (daughter1Probs.get(l)*daughter2Probs.get(k)));
                    }
                    
                    //Old way
                    //lambda = (coal.F.getQuick(k,l)/(coal.Y.getQuick(k)*coal.Y.getQuick(l))) * (daughter1Probs.get(k)*daughter2Probs.get(l) + daughter2Probs.get(l)*daughter1Probs.get(k));
                    
                    coalRates.setQuick(k,l,lambda);
                    lambdaSum += lambda;
                }
            }
                
            //Compute prob of parent of being in each state k after coal event
            //ArrayList<Double> postCoalStateProbs = new ArrayList<Double>(); //for debugging
            double totalProbStateK;
            //double totalProbAllStates = 0.0;
            for (int k = 0; k < states; k++) {
                totalProbStateK = 0.0;
                for (int l = 0; l < states; l++) {
                    totalProbStateK += coalRates.getQuick(k,l);
                }
                totalProbStateK = totalProbStateK / lambdaSum;
                //postCoalStateProbs.add(totalProbStateK);
                probs.get(time-1).get(mapIndexParent).set(k,totalProbStateK);
            }
            
            //Remove daughter lineages from and add parent lineage to activeLineages
            int listIndex1 = activeLineages.indexOf(daughterLineage1Index);
            if (listIndex1 != -1) {
                activeLineages.remove(listIndex1);
            } else {
                System.out.println("Cannot remove daughter lineage");
            }

            int listIndex2 = activeLineages.indexOf(daughterLineage2Index);
            if (listIndex2 != -1) {
                activeLineages.remove(listIndex2);
            } else {
                System.out.println("Cannot remove daughter lineage");
            }
            
            if (coalNode != rootNode) { //Do not add root lineage
                activeLineages.add(coalNode);
            }

        }
    }
    
    public void updateStateSums(int time)
    {
        double currStateSum;
        for (int st = 0; st < states; st++) {
            currStateSum = 0.0;
            for (int lin = 0; lin < activeLineages.size(); lin++) {
                //if (lin >= probs.get(time).size()) {
                    //System.out.println("WTF");
                //}
                //if (st >= probs.get(time).get(lin).size()) {
                    //System.out.println("WTF");
                //}
                currStateSum += probs.get(time).get(lin).get(st);
            }
            stateSums.setQuick(st, currStateSum);
        }
        
        //Old way
        //for (int lin = 0; lin < activeLineages.size(); lin++) {
            //currLineage = activeLineages.get(lin);
            //currLineageIndex = lineageMap.get(time).indexOf(currLineage);
            //for (int st = 0; st < states; st++) {
                //newSum = stateSums.getQuick(st) + probs.get(time).get(currLineageIndex).get(st);
                //stateSums.setQuick(st, newSum);
            //}
        //}
        
    }
    
   public DoubleArrayList getLineagesByState(int time)
    {
        DoubleArrayList linesAk = new DoubleArrayList();
        int numLines = lineageMap.get(time).size();
        for (int st = 0; st < states; st++) {
            double stateSum = 0.0;
            for (int lin = 0; lin < numLines; lin++) {
                stateSum += probs.get(time).get(lin).get(st);
            }
            linesAk.add(stateSum);
        }
        return linesAk;
        
    }
    
    /**
     * THESE METHODS OR FOR COMPUTING PARTICLE-SPECIFIC PROBS 
     */
    
    public void setUpParticleProbArray(ArrayList<ArrayList<Integer>> uniqueParticles) 
    {
        
        //Set up array for particle lineage state probs
        particleProbs = new ArrayList<ArrayList<ArrayList<ArrayList<Double>>>>();
        int times = lineageMap.size();
        int lineagesAtTime = 0; int particlesAtTime;
        int numUniqueParticles = uniqueParticles.get(times-1).size();
        for (int t = 0; t < times; t++) {
            particleProbs.add(new ArrayList<ArrayList<ArrayList<Double>>>());
            //particlesAtTime = uniqueParticles.get(t).size(); 
            for (int p = 0; p < numUniqueParticles; p++) {
                particleProbs.get(t).add(new ArrayList<ArrayList<Double>>());
                lineagesAtTime = lineageMap.get(t).size();
                for (int lin = 0; lin < lineagesAtTime; lin++) {
                    particleProbs.get(t).get(p).add(new ArrayList<Double>());
                    for (int s = 0; s < states; s++) {
                        particleProbs.get(t).get(p).get(lin).add(0.0);
                    }
                }
            }
        }

        for (int st = 0; st < states; st++) {
            stateSums.add(0.0); 
        }
        
        //Set up tranProb matrix
        DoubleFactory2D factory2D;
        factory2D = DoubleFactory2D.dense;
        transProbs = factory2D.make(states, states);
        
    }
    
    public void setUpMiniParticleProbArray(ArrayList<Integer> finalParticles) 
    {
        
        //Set up array for particle lineage state probs (only for unique particles)
        particleMiniProbs = new ArrayList<ArrayList<ArrayList<Double>>>();
        //int times = lineageMap.size();
        int numUniqueParticles = finalParticles.size();
        for (int p = 0; p < numUniqueParticles; p++) {
            particleMiniProbs.add(new ArrayList<ArrayList<Double>>());
        }
        
        for (int st = 0; st < states; st++) {
            stateSums.add(0.0); 
        }
        
        //Set up tranProb matrix
        DoubleFactory2D factory2D;
        factory2D = DoubleFactory2D.dense;
        transProbs = factory2D.make(states, states);
        
    }
    
    public void addParticleSamples(ZVectors dataZ, int zLoc, TreeNode[] tree, int time) 
    {
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        int particlesAtTime = particleProbs.get(time).size();
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                int newLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                int newLineageIndex = lineageMap.get(time).indexOf(newLineage);
                activeLineages.add(newLineage);
                int newState = tree[newLineage].nodeLabel - 1; //nodeLabels are indexed from one
                if (newState == -1) {
                    System.out.println("Hit unknown lineage states");
                }
                for (int p = 0; p < particlesAtTime; p++) {
                    particleProbs.get(time).get(p).get(newLineageIndex).set(newState, 1.0);
                }
            }
        }    
    }
    
    public void addParticleMiniSamples(ZVectors dataZ, int zLoc, TreeNode[] tree, int time) 
    {
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        int uniqueParticles = particleMiniProbs.size();
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                int newLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                //int newLineageMapIndex = lineageMap.get(time).indexOf(newLineage);
                activeLineages.add(newLineage);
                int newActiveLineagesIndex = activeLineages.size() - 1;
                //if (newLineageMapIndex != newActiveLineagesIndex) {
                    //System.out.println("Lineage map and active lineage indexes are off!!");
                //}
                int newState = tree[newLineage].nodeLabel - 1; //nodeLabels are indexed from one
                if (newState == -1) {
                    System.out.println("Hit unknown lineage states");
                }
                for (int p = 0; p < uniqueParticles; p++) {
                    particleMiniProbs.get(p).add(new ArrayList<Double>()); //add new array for lineage
                    for (int st = 0; st < states; st++) {
                        particleMiniProbs.get(p).get(newActiveLineagesIndex).add(0.0);
                    }
                    particleMiniProbs.get(p).get(newActiveLineagesIndex).set(newState, 1.0);
                }
            }
        }    
    }
    
    public void addParticleMiniSamplesWithPriors(ZVectors dataZ, int zLoc, TreeNode[] tree, int time) 
    {
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        int uniqueParticles = particleMiniProbs.size();
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                int newLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                activeLineages.add(newLineage);
                int newActiveLineagesIndex = activeLineages.size() - 1;
                for (int p = 0; p < uniqueParticles; p++) {
                    particleMiniProbs.get(p).add(new ArrayList<Double>()); //add new array for lineage
                    for (int st = 0; st < states; st++) {
                        particleMiniProbs.get(p).get(newActiveLineagesIndex).add(tree[newLineage].statePriors.get(st));
                    }
                }
            }
        }    
    }
    
    public void updateParticleProbs(DoubleMatrix2D xCurr, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime, ArrayList<Integer> ancestors, ArrayList<Integer> uniqueParticles) 
    {
        
        int numLineages = activeLineages.size(); //int timeIndex;
        double probOut; double rateOutByBirth; double rateOutByMigration; double rowSum; 
        double probMassIntoL; double newSk; double newSl;
        int lineIndex; int prevLineageIndex;
        
        //Get previous lineage indexes for current activeLineages
        ArrayList<Integer> prevLineageIndexes = new ArrayList<Integer>();
        ArrayList<Integer> currLineageIndexes = new ArrayList<Integer>();
        for (int lin = 0; lin < numLineages; lin++) {
            int currLine = activeLineages.get(lin);
            prevLineageIndexes.add(lineageMap.get(time).indexOf(currLine));
            currLineageIndexes.add(lineageMap.get(time-1).indexOf(currLine));   
        }
        
        //if (uniqueParticles.get(time-1).size() != uniqueParticles.get(time).size()) {
            //System.out.println("Problem" + currTime);
        //} else {
            //System.out.println("No problem" + currTime);
        //}
        
        /**
         * CHANGE THIS SO ONLY COMPUTING TRANSITION PROBS FOR UNIQUE PARTICLES AT THIS TIME!!
         */
        
        int pIndex;
        int particlesAtTime = uniqueParticles.size();
        for (int x = 0; x < particlesAtTime; x++) {
            
            //For stateProbs at t-1: use x
            //For stateProbs at t: use nextUniqueIndex
            //For xMatrix: use particle index
            pIndex = ancestors.get(uniqueParticles.get(x));
            //childIndexInA = ancestryMatrix.get(time).get(pIndex);
            //nextUniqueX = uniqueParticles.get(time).indexOf(childIndexInA);
            //nextX = uniqueParticles.get(time).indexOf(pIndex);
        
            coal.updateF(xCurr.viewColumn(pIndex), theta, currTime);
            coal.updateG(xCurr.viewColumn(pIndex), theta);
            coal.updateY(xCurr.viewColumn(pIndex));

            //Sum over S's to get A's here - this could be done automatically recomputed at the end of each time step
            //if (nextX == -1) {
                //System.out.println("At-1: " + ancestryMatrix.get(time-1));
                //System.out.println("At: " + ancestryMatrix.get(time));
                //System.out.println("Ut-1: " + uniqueParticles.get(time-1));
                //System.out.println("Ut: " + uniqueParticles.get(time));
                //System.out.println();
            //}
            this.updateParticleStateSums(time,x);

            //Precompute transition probabilties between all states and store these for forwards pass
            for (int k = 0; k < states; ++k) { //For each state or subpopulation
                rowSum = 0.0;
                for (int l = 0; l < states; ++l) { //<loop through l here so you can reuse the following four methods again on the up pass
                    if (k != l) {

                        //New way for dealing with Y's = 0
                        if (coal.Y.getQuick(k) == 0.0 || coal.Y.getQuick(l) == 0.0) {
                            rateOutByBirth = 0.0;
                        } else {
                            if (coal.Y.getQuick(l) > stateSums.get(l)) {
                                rateOutByBirth = (coal.F.getQuick(l,k)/coal.Y.get(k)) * ((coal.Y.getQuick(l)-stateSums.get(l))/coal.Y.getQuick(l));
                            } else {
                                rateOutByBirth = 0.0;
                            }
                        }
                        if (coal.Y.getQuick(k) == 0.0) {
                            rateOutByMigration = 0.0;
                        } else {
                            rateOutByMigration = (coal.G.getQuick(l,k)/coal.Y.getQuick(k));
                        }

                        probOut = 1 - Math.exp(-(rateOutByBirth + rateOutByMigration)*dtTime);
                        transProbs.setQuick(k, l, probOut);
                        rowSum += probOut;

                    }
                }
                transProbs.setQuick(k, k, (1-rowSum)); //Necessary??

            }

            //Compute new cumulative probs for each lineage
            //Update this on 03/11/13 so that indexing should always reflect indexes in arrayMap
            for (int lin = 0; lin < numLineages; lin++) {

                lineIndex = currLineageIndexes.get(lin);
                prevLineageIndex = prevLineageIndexes.get(lin);
                for (int st = 0; st < states; st++) {
                    
                    //double startingProb = particleProbs.get(time).get(nextX).get(prevLineageIndex).get(st);
                    particleProbs.get(time-1).get(x).get(lineIndex).set(st, particleProbs.get(time).get(x).get(prevLineageIndex).get(st));
                    
                }

                for (int k = 0; k < states; ++k) {
                    for (int l = 0; l < states; ++l) {
                        if (k != l) {
                            probMassIntoL = transProbs.getQuick(k, l) * particleProbs.get(time).get(x).get(prevLineageIndex).get(k);
                            newSk = particleProbs.get(time-1).get(x).get(lineIndex).get(k) - probMassIntoL;
                            particleProbs.get(time-1).get(x).get(lineIndex).set(k, newSk);
                            newSl = particleProbs.get(time-1).get(x).get(lineIndex).get(l) + probMassIntoL;
                            particleProbs.get(time-1).get(x).get(lineIndex).set(l, newSl);
                        }
                    }
                }

            }
        
        }
        
    }
    
    public void updateMiniParticleProbs(DoubleMatrix2D xCurr, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime, ArrayList<Integer> ancestors, ArrayList<Integer> uniqueParticles) 
    {
        
        int numLineages = activeLineages.size(); //int timeIndex;
        double probOut; double rateOutByBirth; double rateOutByMigration; double rowSum; 
        double probMassIntoL; double newSk; double newSl;
        //int lineIndex; int prevLineageIndex;
        
        //Get previous lineage indexes for current activeLineages
        /**
         * This is not necessary if order in activeLineages is the same as in the lineageMap
         */
        //ArrayList<Integer> prevLineageIndexes = new ArrayList<Integer>();
        //ArrayList<Integer> currLineageIndexes = new ArrayList<Integer>();
        //for (int lin = 0; lin < numLineages; lin++) {
            //int currLine = activeLineages.get(lin);
            //prevLineageIndexes.add(lineageMap.get(time).indexOf(currLine));
            //currLineageIndexes.add(lineageMap.get(time-1).indexOf(currLine));   
        //}
        
        /**
         * CHANGE THIS SO ONLY COMPUTING TRANSITION PROBS FOR UNIQUE PARTICLES AT THIS TIME!!
         */
        
        DoubleArrayList prevStateProbs = new DoubleArrayList();
        for (int st = 0; st < states; st++) {
            prevStateProbs.add(0.0);
        }
        
        int pIndex;
        int particlesAtTime = uniqueParticles.size();
        for (int x = 0; x < particlesAtTime; x++) {
            
            pIndex = ancestors.get(uniqueParticles.get(x));
        
            coal.updateF(xCurr.viewColumn(pIndex), theta, currTime);
            coal.updateG(xCurr.viewColumn(pIndex), theta);
            coal.updateY(xCurr.viewColumn(pIndex));

            this.updateMiniParticleStateSums(x);

            //Precompute transition probabilties between all states and store these for forwards pass
            for (int k = 0; k < states; ++k) { //For each state or subpopulation
                rowSum = 0.0;
                for (int l = 0; l < states; ++l) { //<loop through l here so you can reuse the following four methods again on the up pass
                    if (k != l) {

                        //New way for dealing with Y's = 0
                        if (coal.Y.getQuick(k) == 0.0 || coal.Y.getQuick(l) == 0.0) {
                            rateOutByBirth = 0.0;
                        } else {
                            if (coal.Y.getQuick(l) > stateSums.get(l)) {
                                rateOutByBirth = (coal.F.getQuick(l,k)/coal.Y.get(k)) * ((coal.Y.getQuick(l)-stateSums.get(l))/coal.Y.getQuick(l));
                            } else {
                                rateOutByBirth = 0.0;
                            }
                        }
                        if (coal.Y.getQuick(k) == 0.0) {
                            rateOutByMigration = 0.0;
                        } else {
                            rateOutByMigration = (coal.G.getQuick(l,k)/coal.Y.getQuick(k));
                        }

                        probOut = 1 - Math.exp(-(rateOutByBirth + rateOutByMigration)*dtTime);
                        transProbs.setQuick(k, l, probOut);
                        rowSum += probOut;

                    }
                }
                //transProbs.setQuick(k, k, (1-rowSum)); //Necessary??

            }

            //Compute new cumulative probs for each lineage
            //Update this on 03/11/13 so that indexing should always reflect indexes in arrayMap
            for (int lin = 0; lin < numLineages; lin++) {

                //lineIndex = currLineageIndexes.get(lin);
                //prevLineageIndex = prevLineageIndexes.get(lin);
                for (int st = 0; st < states; st++) {
                    
                    //Copy over current lineage state probs
                    prevStateProbs.set(st, particleMiniProbs.get(x).get(lin).get(st));
                    
                    //double startingProb = particleProbs.get(time).get(nextX).get(prevLineageIndex).get(st);
                    //particleProbs.get(time-1).get(x).get(lineIndex).set(st, particleProbs.get(time).get(x).get(prevLineageIndex).get(st));
                    
                }

                for (int k = 0; k < states; ++k) {
                    for (int l = 0; l < states; ++l) {
                        if (k != l) {
                            probMassIntoL = transProbs.getQuick(k, l) * prevStateProbs.getQuick(k); //particleProbs.get(time).get(x).get(prevLineageIndex).get(k);
                            newSk = particleMiniProbs.get(x).get(lin).get(k) - probMassIntoL;
                            particleMiniProbs.get(x).get(lin).set(k, newSk);
                            newSl = particleMiniProbs.get(x).get(lin).get(l) + probMassIntoL;
                            particleMiniProbs.get(x).get(lin).set(l, newSl);
                        }
                    }
                }

            }
        
        }
        
    }
    
    public void updateParticleProbsCoalEvent(ZVectors dataZ, int zLoc, TreeNode[] tree, DoubleMatrix2D xCurr, DoubleArrayList theta, int time, StructCoalModel coal, double currTime, ArrayList<Integer> ancestors, ArrayList<Integer> uniqueParticles) 
    {
        
        DoubleFactory2D factory2D = DoubleFactory2D.dense;
        DoubleMatrix2D coalRates = factory2D.make(states, states);
        //DoubleMatrix2D otherCoalRates = factory2D.make(states, states);
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        int coalNode; int daughterLineage1Index; int daughterLineage2Index;
        int mapIndexDaughter1; int mapIndexDaughter2; int mapIndexParent;
        double lambda; double lambdaSum;
        double totalProbStateK;
        
        //Find coalescent events at this time
        ArrayList<Integer> coalNodesAtEvent = new ArrayList<Integer>();
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event)==1) {
                coalNode = dataZ.nodePointers.get(zLoc).get(event);
                coalNodesAtEvent.add(coalNode);
            }
        }
        
        int coalEvents = coalNodesAtEvent.size();
        for (int event = 0; event < coalEvents; event++) {
                
            coalNode = coalNodesAtEvent.get(event);
            daughterLineage1Index = tree[coalNode].childNodes[0] - 1;
            daughterLineage2Index = tree[coalNode].childNodes[1] - 1;
    
            mapIndexDaughter1 = lineageMap.get(time).indexOf(daughterLineage1Index);
            mapIndexDaughter2 = lineageMap.get(time).indexOf(daughterLineage2Index);
            mapIndexParent = lineageMap.get(time).indexOf(coalNode);
            
            //For each particle
            int pIndex;
            //int nextX;
            int particlesAtTime = uniqueParticles.size();
            for (int x = 0; x < particlesAtTime; x++) {
                
                //For stateProbs at t: use x
                //For xMatrix: use particle index
                pIndex = ancestors.get(uniqueParticles.get(x));
                
                coal.updateF(xCurr.viewColumn(pIndex), theta, currTime); //xCurr from time-1
                coal.updateG(xCurr.viewColumn(pIndex), theta);
                coal.updateY(xCurr.viewColumn(pIndex));
                
                ArrayList<Double> daughter1Probs = particleProbs.get(time).get(x).get(mapIndexDaughter1);
                ArrayList<Double> daughter2Probs = particleProbs.get(time).get(x).get(mapIndexDaughter2);
            
                lambdaSum = 0.0;
                for (int k = 0; k < states; k++) {
                    for (int l = 0; l < states; l++) {

                        //New way for dealing with Y's = 0
                        if (coal.Y.getQuick(k) == 0.0 || coal.Y.getQuick(l) == 0.0) {
                            lambda = 0.0;
                        } else {
                            //THERE WAS A BUG HERE WITH HOW LAMBDA WAS BEING COMPUTED WITH THE INDEXING OF K AND L IN DAUGTHER PROBS
                            //lambda = (coal.F.getQuick(k,l)/(coal.Y.getQuick(k)*coal.Y.getQuick(l))) * (daughter1Probs.get(k)*daughter2Probs.get(l) + daughter2Probs.get(l)*daughter1Probs.get(k));
                            lambda = (coal.F.getQuick(k,l)/(coal.Y.getQuick(k)*coal.Y.getQuick(l))) * ((daughter1Probs.get(k)*daughter2Probs.get(l)) + (daughter1Probs.get(l)*daughter2Probs.get(k)));
                        }

                        //Old way
                        //lambda = (coal.F.getQuick(k,l)/(coal.Y.getQuick(k)*coal.Y.getQuick(l))) * (daughter1Probs.get(k)*daughter2Probs.get(l) + daughter2Probs.get(l)*daughter1Probs.get(k));

                        coalRates.setQuick(k,l,lambda);
                        lambdaSum += lambda;
                    }
                }

                //Compute prob of parent of being in each state k after coal event
                //ArrayList<Double> postCoalStateProbs = new ArrayList<Double>(); //for debugging
                //double totalProbAllStates = 0.0;
                for (int k = 0; k < states; k++) {
                    totalProbStateK = 0.0;
                    for (int l = 0; l < states; l++) {
                        totalProbStateK += coalRates.getQuick(k,l);
                    }
                    totalProbStateK = totalProbStateK / lambdaSum;
                    //postCoalStateProbs.add(totalProbStateK);
                    particleProbs.get(time).get(x).get(mapIndexParent).set(k,totalProbStateK);
                }
            
            } //end FOR particles
            
            //Remove daughter lineages from and add parent lineage to activeLineages
            int listIndex1 = activeLineages.indexOf(daughterLineage1Index);
            if (listIndex1 != -1) {
                activeLineages.remove(listIndex1);
            } else {
                System.out.println("Cannot remove daughter lineage");
            }

            int listIndex2 = activeLineages.indexOf(daughterLineage2Index);
            if (listIndex2 != -1) {
                activeLineages.remove(listIndex2);
            } else {
                System.out.println("Cannot remove daughter lineage");
            }
            
            if (coalNode != rootNode) { //Do not add root lineage
                activeLineages.add(coalNode);
            }

        }
    }
    
    public void updateMiniParticleProbsCoalEvent(ZVectors dataZ, int zLoc, TreeNode[] tree, DoubleMatrix2D xCurr, DoubleArrayList theta, int time, StructCoalModel coal, double currTime, ArrayList<Integer> ancestors, ArrayList<Integer> uniqueParticles) 
    {
        
        DoubleFactory2D factory2D = DoubleFactory2D.dense;
        DoubleMatrix2D coalRates = factory2D.make(states, states);
        //DoubleMatrix2D otherCoalRates = factory2D.make(states, states);
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        int coalNode; int daughterLineage1Index; int daughterLineage2Index;
        int mapIndexDaughter1; int mapIndexDaughter2; int mapIndexParent;
        double lambda; double lambdaSum;
        double totalProbStateK;
        
        //Find coalescent events at this time
        ArrayList<Integer> coalNodesAtEvent = new ArrayList<Integer>();
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event)==1) {
                coalNode = dataZ.nodePointers.get(zLoc).get(event);
                coalNodesAtEvent.add(coalNode);
            }
        }
        
        int coalEvents = coalNodesAtEvent.size();
        for (int event = 0; event < coalEvents; event++) {
                
            coalNode = coalNodesAtEvent.get(event);
            daughterLineage1Index = tree[coalNode].childNodes[0] - 1;
            daughterLineage2Index = tree[coalNode].childNodes[1] - 1;
    
            //mapIndexDaughter1 = lineageMap.get(time).indexOf(daughterLineage1Index);
            //mapIndexDaughter2 = lineageMap.get(time).indexOf(daughterLineage2Index);
            //mapIndexParent = lineageMap.get(time).indexOf(coalNode);
            
            ////Get daughter and parent indexes in activeLineages
            mapIndexDaughter1 = activeLineages.indexOf(daughterLineage1Index);
            mapIndexDaughter2 = activeLineages.indexOf(daughterLineage2Index);
            //mapIndexParent = activeLineages.indexOf(coalNode);
            
            //For each particle
            int pIndex;
            //int nextX;
            int particlesAtTime = uniqueParticles.size();
            for (int x = 0; x < particlesAtTime; x++) {
                
                //For stateProbs at t: use x
                //For xMatrix: use particle index
                pIndex = ancestors.get(uniqueParticles.get(x));
                
                coal.updateF(xCurr.viewColumn(pIndex), theta, currTime); //xCurr from time-1
                coal.updateG(xCurr.viewColumn(pIndex), theta);
                coal.updateY(xCurr.viewColumn(pIndex));
                
                ArrayList<Double> daughter1Probs = particleMiniProbs.get(x).get(mapIndexDaughter1);
                ArrayList<Double> daughter2Probs = particleMiniProbs.get(x).get(mapIndexDaughter2);
            
                lambdaSum = 0.0;
                for (int k = 0; k < states; k++) {
                    for (int l = 0; l < states; l++) {

                        //New way for dealing with Y's = 0
                        if (coal.Y.getQuick(k) == 0.0 || coal.Y.getQuick(l) == 0.0) {
                            lambda = 0.0;
                        } else {
                            //THERE WAS A BUG HERE WITH HOW LAMBDA WAS BEING COMPUTED WITH THE INDEXING OF K AND L IN DAUGTHER PROBS
                            //lambda = (coal.F.getQuick(k,l)/(coal.Y.getQuick(k)*coal.Y.getQuick(l))) * (daughter1Probs.get(k)*daughter2Probs.get(l) + daughter2Probs.get(l)*daughter1Probs.get(k));
                            lambda = (coal.F.getQuick(k,l)/(coal.Y.getQuick(k)*coal.Y.getQuick(l))) * ((daughter1Probs.get(k)*daughter2Probs.get(l)) + (daughter1Probs.get(l)*daughter2Probs.get(k)));
                        }

                        //Old way
                        //lambda = (coal.F.getQuick(k,l)/(coal.Y.getQuick(k)*coal.Y.getQuick(l))) * (daughter1Probs.get(k)*daughter2Probs.get(l) + daughter2Probs.get(l)*daughter1Probs.get(k));

                        coalRates.setQuick(k,l,lambda);
                        lambdaSum += lambda;
                    }
                }
                
                //Remove arrays for daughter lineages --- removing larger indexed lineage first to retain correct index of the smaller indexed lineage
                if (mapIndexDaughter1 > mapIndexDaughter2) {
                    particleMiniProbs.get(x).remove(mapIndexDaughter1);
                    particleMiniProbs.get(x).remove(mapIndexDaughter2);
                } else {
                    particleMiniProbs.get(x).remove(mapIndexDaughter2);
                    particleMiniProbs.get(x).remove(mapIndexDaughter1);
                }
                

                //Compute prob of parent of being in each state k after coal event
                particleMiniProbs.get(x).add(new ArrayList<Double> ());
                for (int k = 0; k < states; k++) {
                    totalProbStateK = 0.0;
                    for (int l = 0; l < states; l++) {
                        totalProbStateK += coalRates.getQuick(k,l);
                    }
                    totalProbStateK = totalProbStateK / lambdaSum;
                    particleMiniProbs.get(x).get(activeLineages.size() - 2).add(totalProbStateK);
                }
            
            } //end FOR particles
            
            //Remove daughter lineages from and add parent lineage to activeLineages
            int listIndex1 = activeLineages.indexOf(daughterLineage1Index);
            if (listIndex1 != -1) {
                activeLineages.remove(listIndex1);
            } else {
                System.out.println("Cannot remove daughter lineage");
            }

            int listIndex2 = activeLineages.indexOf(daughterLineage2Index);
            if (listIndex2 != -1) {
                activeLineages.remove(listIndex2);
            } else {
                System.out.println("Cannot remove daughter lineage");
            }
            
            if (coalNode != rootNode) { //Do not add root lineage
                activeLineages.add(coalNode);
            }

        }
    }
    
    public void updateParticleStateSums(int time, int particle)
    {
        for (int st = 0; st < states; st++) {
            stateSums.setQuick(st, 0.0); 
        }
        double newSum; int currLineage; int currLineageIndex;
        for (int lin = 0; lin < activeLineages.size(); lin++) {
            currLineage = activeLineages.get(lin);
            currLineageIndex = lineageMap.get(time).indexOf(currLineage);
            for (int st = 0; st < states; st++) {
                newSum = stateSums.getQuick(st) + particleProbs.get(time).get(particle).get(currLineageIndex).get(st);
                stateSums.setQuick(st, newSum);
            }
        }
        
    }
    
    public void updateMiniParticleStateSums(int particle)
    {
        
        double currStateSum;
        for (int st = 0; st < states; st++) {
            currStateSum = 0.0;
            for (int lin = 0; lin < activeLineages.size(); lin++) {
                currStateSum += particleMiniProbs.get(particle).get(lin).get(st);
            }
            stateSums.setQuick(st, currStateSum);
        }
        
    }
    
    public DoubleArrayList getParticleLineagesByState(int time, int particle)
    {
        DoubleArrayList linesAk = new DoubleArrayList();
        int numLines = lineageMap.get(time).size();
        for (int st = 0; st < states; st++) {
            double stateSum = 0.0;
            for (int lin = 0; lin < numLines; lin++) {
                stateSum += particleProbs.get(time).get(particle).get(lin).get(st);
            }
            linesAk.add(stateSum);
        }
        return linesAk;
        
    }
    
    public DoubleArrayList getMiniParticleLineagesByState(int particle)
    {
        DoubleArrayList linesAk = new DoubleArrayList();
        int numLines = activeLineages.size();
        double currStateSum;
        for (int st = 0; st < states; st++) {
            currStateSum = 0.0;
            for (int lin = 0; lin < numLines; lin++) {
                currStateSum += particleMiniProbs.get(particle).get(lin).get(st);
            }
            linesAk.add(currStateSum);
        }
        return linesAk;
        
    }
    
}
