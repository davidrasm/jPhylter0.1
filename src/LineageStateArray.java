import cern.colt.list.*;
import cern.colt.matrix.*;
import java.util.ArrayList;

/**
 *
 * @author David
 */
public class LineageStateArray {
    
    int particles;
    int states;
    int lineages;
    ArrayList<ArrayList<Integer>> array = new ArrayList<ArrayList<Integer>>();
    ArrayList<Integer> activeLineages = new ArrayList<Integer>();
    ArrayList<Integer> lineagesRemoved = new ArrayList<Integer>();
    ArrayList<Integer> lineagesAdded = new ArrayList<Integer>();
    
    public void getArray(int jParticles, int lines, int sts) {
        
        states = sts;
        particles = jParticles;
        lineages = lines;
        //Array row correspond to lineages, columns to particles
        for (int j = 0; j < lines; j++) {
            array.add(new ArrayList<Integer>());
            for (int k = 0; k < particles; k++) {
                array.get(j).add(-1);
            }
        }
        
    }
    
    /**
    public void addCoalLineages(ZVectors dataZ, int zLoc, TreeNode[] tree, DoubleMatrix2D xCurr, DoubleArrayList theta, StructCoalModel coal)
    {
    
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        double currTime = dataZ.absoluteTimes.get(zLoc);
        for (int event = 0; event < totalEvents; event++) {
            
            if (dataZ.omegaEvents.get(zLoc).get(event) == 1) { //if coalescent event
                
                int coalLineage = dataZ.nodePointers.get(zLoc).get(event);
                int daughterLineage1 = tree[coalLineage].childNodes[0] - 1;
                int daughterLineage2 = tree[coalLineage].childNodes[1] - 1;
                //Add daughter lineages to activeLineages and remove parent
                activeLineages.add(daughterLineage1);
                activeLineages.add(daughterLineage2);
                int listIndexCoalLineage = activeLineages.indexOf(coalLineage);
                if (listIndexCoalLineage == -1) {
                    //DoubleArrayList coalTimes = TreeUtils.getInternalNodeTimes(tree);
                    System.out.println("Cannot remove lineage");
                } else {
                    activeLineages.remove(listIndexCoalLineage);
                }
                
                for (int x = 0; x < particles; x++) {
                
                    coal.updateF(xCurr.viewColumn(x), theta, currTime);
                    coal.updateG(xCurr.viewColumn(x), theta);
                    coal.updateY(xCurr.viewColumn(x));
                    double totalRateSum = 0.0;
                    for (int stateK = 0; stateK < states; stateK++) {
                        for (int stateL = 0; stateL < states; stateL++) {
                            double thisRate = (coal.F.get(stateL, stateK) / coal.Y.get(stateL)) * matrix.get(x, coalLineage, stateL);
                            totalRateSum += thisRate; 
                        }
                    }

                    for (int stateK = 0; stateK < states; stateK++) {
                        
                        double nonTransProbInStateK = 0.0;
                        for (int stateL = 0; stateL < states; stateL++) {
                            double thisRate = ((coal.F.get(stateL, stateK) / coal.Y.get(stateL)) * matrix.get(x, coalLineage, stateL));
                            nonTransProbInStateK += thisRate;
                        }
                        
                        double transProbInStateK = 0.0;
                        for (int stateL = 0; stateL < states; stateL++) {
                            double thisRate = ((coal.F.get(stateK, stateL) / coal.Y.get(stateK)) * matrix.get(x, coalLineage, stateK));
                            transProbInStateK += thisRate;
                        }
                        
                        double probStateK = (nonTransProbInStateK + transProbInStateK) / (2*totalRateSum);
                        
                        //Assign probabilities to daughter lineages
                        matrix.set(x, daughterLineage1, stateK, probStateK);
                        matrix.set(x, daughterLineage2, stateK, probStateK);
                    }
                    //System.out.println();
                } //END FOR particles    
            }
        }    
    }
        
    public void removeSamples(ZVectors dataZ, int zLoc) 
    {
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        double eventTime = dataZ.absoluteTimes.get(zLoc);
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                
                int sampledLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                int sampleListIndex = activeLineages.indexOf(sampledLineage);
                activeLineages.remove(sampleListIndex);
                
            }
        }      
    }*/
    
    public void addSamples(ZVectors dataZ, int zLoc, TreeNode[] tree) 
    {
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        //double eventTime = dataZ.absoluteTimes.get(zLoc);
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                
                int newLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                activeLineages.add(newLineage);
                lineagesAdded.add(newLineage);
                int newState = tree[newLineage].nodeLabel - 1; //nodeLabels are indexed from one
                for (int j = 0; j < particles; j++) {
                    if (newState == -1) {
                        System.out.println("Hit unknown lineage states");
                    }
                    //matrix.set(j, newLineage, newState, 1.0); //set new state prob to 1.0
                    array.get(newLineage).set(j, newState);
                }
            }
        }
        
    }
    
    public void updateStates(DoubleMatrix2D xCurr, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime) 
    {
        
        int numLineages = activeLineages.size();
        int i; int currState; int nextState; double birthRateOut; double migrationRateOut; int linesInLNotInA;
        DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
	DoubleMatrix2D transMatrix = factory2D.make(states,states);
        DoubleMatrix2D normCumTransMatrix = factory2D.make(states,states);
        
        //For each particle
        for (int x = 0; x < particles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j
        
            coal.updateF(xCurr.viewColumn(x), theta, currTime);
            coal.updateG(xCurr.viewColumn(x), theta);
            coal.updateY(xCurr.viewColumn(x));
            
            //Precompute A(l) terms by summing probs in matrix
            ArrayList<Integer> lineagesInAByState = this.getLineagesByState(x);
            
            //Can precompute lineage transition probabilties here
            double rowSum; double pMove;
            for (int k = 0; k < states; k++) {
                rowSum = 0.0;
                for (int l = 0; l < states; l++) {
                    if (l != k) {
                        if (coal.Y.getQuick(k) > 0.0 & coal.Y.getQuick(l) > 0.0) {
                            linesInLNotInA = (int) coal.Y.getQuick(l) - lineagesInAByState.get(l);
                            if (linesInLNotInA > 0) {
                                birthRateOut = (coal.F.getQuick(l,k)/coal.Y.getQuick(k))*((linesInLNotInA) / coal.Y.getQuick(l));
                            } else {
                                birthRateOut = 0.0;
                            }
                        } else {
                            birthRateOut = 0.0;
                        }
                        if (coal.Y.getQuick(k) > 0.0) {
                            migrationRateOut = coal.G.getQuick(l,k)/coal.Y.getQuick(k);
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
                double cumSum = 0;
                for (int l = 0; l < states; l++) {
                    cumSum += transMatrix.getQuick(k,l);
                    normCumTransMatrix.setQuick(k,l,cumSum);
                }
            }
            //System.out.println();
            
            for (int lin = 0; lin < numLineages; lin++) { //for all lineages i
                
                i = activeLineages.get(lin);
                currState = array.get(i).get(x);
                nextState = this.getNextState(currState, normCumTransMatrix);
                array.get(i).set(x,nextState);
                
            }
        }
    }
    
    public void updateStatesAfterCoal(int x, StructCoalModel coal, int coalNode, int leftDaughterState, int rightDaughterState) 
    {
        
        //Probabilistically choose a state for the parent of the two daughter lineages
        double fKtoL = coal.F.get(leftDaughterState, rightDaughterState);
        double fLtoK = coal.F.get(rightDaughterState, leftDaughterState);
        double pKtoL = fKtoL / (fKtoL + fLtoK);  //probability left daughter transmitted
        double uniRand = Math.random();
        int coalNodeState = -1;
        if (uniRand <= pKtoL) {   
            coalNodeState = leftDaughterState;   
        } else {
            coalNodeState = rightDaughterState;    
        }
        array.get(coalNode).set(x, coalNodeState);
        //System.out.println(); 
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
    
    public ArrayList<Integer> getLineagesByState(int x) 
    {
        ArrayList<Integer> lineagesByState = new ArrayList<Integer>();
        int numLineages = activeLineages.size();
        int i; int stateSum;
        for (int k = 0; k < states; k++) {
            stateSum = 0;
            for (int lin = 0; lin < numLineages; lin++) {
                i = activeLineages.get(lin);
                if (array.get(i).get(x) == k) {
                    stateSum++;
                }
            }
            lineagesByState.add(stateSum);
        }
        return lineagesByState;
    }
    
    
    
    
}