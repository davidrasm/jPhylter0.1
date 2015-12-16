import cern.colt.matrix.*;
import cern.colt.list.*;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Superclass for the structured coalescent model
 * @author David
 */
public class StructCoalModel {
    
    DoubleMatrix2D F;
    DoubleMatrix2D G;
    DoubleMatrix1D Y;
    int states;
    int transCount;
    int nonTransCount;
    ArrayList<ArrayList<Integer>> transCountArray = new ArrayList<ArrayList<Integer>>();
    double currTime;
    DoubleArrayList recordedTimes = new DoubleArrayList();

    
    public void make() 
    {   
        
        /**
         * 
         * NEED TO UPDATE THIS FOR EACH MODEL!!
         * 
         */
        
    }
    

    public void updateF(DoubleMatrix1D xCurr, DoubleArrayList theta, double currTime) 
    {
        /**
         * 
         * NEED TO UPDATE THIS FOR EACH MODEL!!
         * 
         */
        
    }
    
    public void updateG(DoubleMatrix1D xCurr, DoubleArrayList theta)
    {
         /**
         * 
         * NEED TO UPDATE THIS FOR EACH MODEL!!
         * 
         */
    }
    
    public void updateY(DoubleMatrix1D xCurr) 
    {
         /**
         * 
         * NEED TO UPDATE THIS FOR EACH MODEL!!
         * 
         */
    }
    
    public double computeCoalRate(int k, int l) 
    {
        double rate = (F.getQuick(k,l) + F.getQuick(l,k)) / (Y.getQuick(k)*Y.getQuick(l));
//        double rate = 0.0;
//        if (k == 0 & l == 0) {
//            rate = F.getQuick(k,l);
//        } else {
//            rate = (F.getQuick(k,l) + F.getQuick(l,k)) / (Y.getQuick(k)*Y.getQuick(l));
//        }
        return rate; 
    }
    
    public double computeCoalRateOneDirection(int k, int l) 
    {
        double rate = F.getQuick(k,l) / (Y.getQuick(k)*Y.getQuick(l));
        return rate;
        
    }
    
    public double computeMoveRate(int k, int l, int[] lineagesByState)
    {
        
        double rate = (G.getQuick(l,k)/Y.getQuick(k)) + (F.getQuick(l,k)/Y.getQuick(k))*((Y.getQuick(l)-lineagesByState[l])/Y.getQuick(l));
        return rate;
    }
    
    public double computeTransProbs(DoubleMatrix1D xCurr, ParticleLineageStateArray stateArray, DoubleArrayList theta, int time, double dtTime, double currTime, ArrayList<Integer> eventTypes)
    {

        int numLineages = stateArray.arrayMap.get(time).size(); //actual number of lineages present between time t and t-1 might be less than this
        DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
        DoubleMatrix2D birthRateMatrix = factory2D.make(states,states);
        DoubleMatrix2D migrationRateMatrix = factory2D.make(states,states);
        DoubleMatrix2D transMatrix = factory2D.make(states, states);
        double birthRateOut; double migrationRateOut; double rowSum; double moveRate; int currState; int nextState;
        ArrayList<Integer> currTimeLineageStates = new ArrayList<Integer>();
        ArrayList<Integer> nextTimeLineageStates = new ArrayList<Integer>();
        int activeLineages = 0;
        ArrayList<Integer> lineagesInAByState = new ArrayList<Integer>(); int linesInLNotInA;
        
        double transProb = 1.0;
        if (numLineages >= 1) {
            
            this.updateF(xCurr, theta, currTime);
            this.updateG(xCurr, theta);
            this.updateY(xCurr);
            for (int k = 0; k < states; k++) {
                for (int l = 0; l < states; l++) {
                    if (l != k) {
                        if (Y.getQuick(k) == 0.0 || Y.getQuick(l) == 0.0) {
                            birthRateMatrix.setQuick(k,l,0.0);
                            migrationRateMatrix.setQuick(k,l,0.0);
                        } else {
                            birthRateOut = (F.getQuick(l,k)/Y.getQuick(k))*(1.0 / Y.getQuick(l));
                            birthRateMatrix.setQuick(k,l,birthRateOut);
                            migrationRateOut = G.getQuick(l,k) / Y.getQuick(k); //could argue that in general only Y(k) has to be > 0.0 
                            migrationRateMatrix.setQuick(k,l,migrationRateOut);
                        }
                    }
                }
            }
            
            if (eventTypes.contains(1)) { //if a coalescent event occured
                //Coalescent event occured so some lineages present in stateArray at time t will NOT be present at time t-1
                int currTimeLineageIndex; int indexOfLineageAtNextTime;
                for (int lin = 0; lin < numLineages; lin++) { //Find lineages present at time t and time t-1 and get their states
                    currTimeLineageIndex = stateArray.arrayMap.get(time).get(lin);
                    if (stateArray.arrayMap.get(time-1).contains(currTimeLineageIndex)) {
                        indexOfLineageAtNextTime = stateArray.arrayMap.get(time-1).indexOf(currTimeLineageIndex);
                        currTimeLineageStates.add(stateArray.sampleArray.get(time).get(lin));
                        nextTimeLineageStates.add(stateArray.sampleArray.get(time-1).get(indexOfLineageAtNextTime));
                    }
                }
                activeLineages = currTimeLineageStates.size();
                //Compute lineagesInAByState for continuing lineages
                int stateCount;
                for (int st = 0; st < states; st++) {
                    stateCount = 0;
                    for (int lin = 0; lin < activeLineages; lin++) {
                        if (currTimeLineageStates.get(lin) == st) {
                            stateCount++;
                        }
                    }
                    lineagesInAByState.add(stateCount);
                }
                
            } else { //if no coalescent event occured
                lineagesInAByState = stateArray.sampleArrayLineageStates.get(time);
                activeLineages = stateArray.arrayMap.get(time).size();
                currTimeLineageStates = stateArray.sampleArray.get(time);
                nextTimeLineageStates = stateArray.sampleArray.get(time-1);
                //if (currTimeLineageStates.size() != nextTimeLineageStates.size()) {
                    //System.out.println("Found lineage state index mismatch");
                //}
                
            }
            
            //Check if population size is compatible with lineage state mapping
            boolean sizeFlag = false;
            for (int st = 0; st < states; st++) {
                if (lineagesInAByState.get(st) > Y.getQuick(st)) {
                    sizeFlag = true;
                }
            }

            if (sizeFlag) {

                transProb = 0.0; //System.out.println("Ak is larger than Yk at time " + currTime);

            } else {
            
                for (int k = 0; k < states; k++) {
                    rowSum = 0.0;
                    for (int l = 0; l < states; l++) {
                        if (l != k) {
                            linesInLNotInA = (int) Y.getQuick(l) - lineagesInAByState.get(l);
                            if (linesInLNotInA > 0) {
                                birthRateOut = birthRateMatrix.getQuick(k,l) * linesInLNotInA;
                            } else {
                                //System.out.println("More lineages in state " + l + " than in pop");
                                birthRateOut = 0.0;
                            }
                            migrationRateOut = migrationRateMatrix.getQuick(k,l);
                            moveRate = birthRateOut + migrationRateOut;
                            rowSum += moveRate;
                            transMatrix.setQuick(k, l, moveRate);
                        }   
                    }
                    transMatrix.setQuick(k, k, (rowSum));
                }

                for (int ln = 0; ln < activeLineages; ln++) { //for all lineages i
                    currState = currTimeLineageStates.get(ln); //stateArray.sampleArray.get(time).get(prevLineageIndexes.get(ln));
                    nextState = nextTimeLineageStates.get(ln); //stateArray.sampleArray.get(time-1).get(ln);
                    transCountArray.get(nextState).set(currState, (transCountArray.get(nextState).get(currState) + 1));
                    if (currState == nextState) {
                        transProb *= Math.exp(-transMatrix.getQuick(currState,currState) * dtTime); //CHANGED FROM Math.exp(-(1-lambda*dt) 
                        //nonTransCount++;
                    } else {
                        transProb *= 1 - Math.exp(-transMatrix.getQuick(currState,nextState) * dtTime);
                        transCount++;
                    }
                    if (nextState == 0 & currState == 1) {
                        //transCount++;
                        recordedTimes.add(currTime);
                    }
                }
            }
        }
        if (transProb <= 0.0) {
            transProb = Double.MIN_VALUE;
            //System.out.println("transProb zero at non-coalescent event");
        }
        if (Double.isNaN(transProb)) {
            transProb = Double.MIN_VALUE;
            //System.out.println("transProb NaN at non-coalescent event");
        }
        return transProb;                
    }
    
    public double computeTransProbsNew(DoubleMatrix1D xCurr, ParticleLineageStateArray stateArray, DoubleArrayList theta, int time, double dtTime, double currTime, ArrayList<Integer> eventTypes)
    {

        int numLineages = stateArray.arrayMap.get(time).size(); //actual number of lineages present between time t and t-1 might be less than this
        DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
        DoubleMatrix2D birthRateMatrix = factory2D.make(states,states);
        DoubleMatrix2D migrationRateMatrix = factory2D.make(states,states);
        DoubleMatrix2D transRates = factory2D.make(states, states);
        DoubleMatrix2D transProbs = factory2D.make(states, states);
        double birthRateOut; double migrationRateOut; double rowSum; double moveRate; int currState; int nextState;
        ArrayList<Integer> currTimeLineageStates = new ArrayList<Integer>();
        ArrayList<Integer> nextTimeLineageStates = new ArrayList<Integer>();
        int activeLineages = 0;
        ArrayList<Integer> lineagesInAByState = new ArrayList<Integer>(); int linesInLNotInA;
        
        double transProb = 1.0;
        if (numLineages >= 1) {
            
            this.updateF(xCurr, theta, currTime);
            this.updateG(xCurr, theta);
            this.updateY(xCurr);
            for (int k = 0; k < states; k++) {
                for (int l = 0; l < states; l++) {
                    //if (l != k) {
                        if (Y.getQuick(k) == 0.0 || Y.getQuick(l) == 0.0) {
                            birthRateMatrix.setQuick(k,l,0.0);
                            migrationRateMatrix.setQuick(k,l,0.0);
                        } else {
                            birthRateOut = (F.getQuick(l,k)/Y.getQuick(k))*(1.0 / Y.getQuick(l));
                            birthRateMatrix.setQuick(k,l,birthRateOut);
                            migrationRateOut = G.getQuick(l,k) / Y.getQuick(k); //could argue that in general only Y(k) has to be > 0.0 
                            migrationRateMatrix.setQuick(k,l,migrationRateOut);
                        }
                    //}
                }
            }
            
            if (eventTypes.contains(1)) { //if a coalescent event occured
                //Coalescent event occured so some lineages present in stateArray at time t will NOT be present at time t-1
                int currTimeLineageIndex; int indexOfLineageAtNextTime;
                for (int lin = 0; lin < numLineages; lin++) { //Find lineages present at time t and time t-1 and get their states
                    currTimeLineageIndex = stateArray.arrayMap.get(time).get(lin);
                    if (stateArray.arrayMap.get(time-1).contains(currTimeLineageIndex)) {
                        indexOfLineageAtNextTime = stateArray.arrayMap.get(time-1).indexOf(currTimeLineageIndex);
                        currTimeLineageStates.add(stateArray.sampleArray.get(time).get(lin));
                        nextTimeLineageStates.add(stateArray.sampleArray.get(time-1).get(indexOfLineageAtNextTime));
                    }
                }
                activeLineages = currTimeLineageStates.size();
                //Compute lineagesInAByState for continuing lineages
                int stateCount;
                for (int st = 0; st < states; st++) {
                    stateCount = 0;
                    for (int lin = 0; lin < activeLineages; lin++) {
                        if (currTimeLineageStates.get(lin) == st) {
                            stateCount++;
                        }
                    }
                    lineagesInAByState.add(stateCount);
                }
                
            } else { //if no coalescent event occured
                lineagesInAByState = stateArray.sampleArrayLineageStates.get(time);
                activeLineages = stateArray.arrayMap.get(time).size();
                currTimeLineageStates = stateArray.sampleArray.get(time);
                nextTimeLineageStates = stateArray.sampleArray.get(time-1);
                //if (currTimeLineageStates.size() != nextTimeLineageStates.size()) {
                    //System.out.println("Found lineage state index mismatch");
                //}
                
            }
            
            //Check if population size is compatible with lineage state mapping
            boolean sizeFlag = false;
            for (int st = 0; st < states; st++) {
                if (lineagesInAByState.get(st) > Y.getQuick(st)) {
                    sizeFlag = true;
                }
            }

            if (sizeFlag) {

                transProb = 0.0; System.out.println("Ak is larger than Yk at time " + currTime);

            } else {
                
                double probOut;
                for (int k = 0; k < states; k++) {
                    rowSum = 0.0;
                    for (int l = 0; l < states; l++) {
                        //if (l != k) {
                            linesInLNotInA = (int) Y.getQuick(l) - lineagesInAByState.get(l);
                            if (linesInLNotInA > 0) {
                                birthRateOut = birthRateMatrix.getQuick(k,l) * linesInLNotInA;
                            } else {
                                //System.out.println("More lineages in state " + l + " than in pop");
                                birthRateOut = 0.0;
                            }
                            migrationRateOut = migrationRateMatrix.getQuick(k,l);
                            moveRate = birthRateOut + migrationRateOut;
                            rowSum += moveRate;
                            transRates.setQuick(k, l, moveRate);
                        //}   
                    }
                    double probEvent = 1 - Math.exp(-rowSum*dtTime);
                    double probNoEvent = Math.exp(-rowSum*dtTime);
                    for (int l = 0; l < states; ++l) {
                        if (l == k) {
                            probOut = probNoEvent + (probEvent*(transRates.getQuick(k,l)/rowSum)); //don't actually need to compute this
                        } else {
                            probOut = probEvent * (transRates.getQuick(k,l)/rowSum);
                        }
                        transProbs.setQuick(k, l, probOut);
                    }
                    //transMatrix.setQuick(k, k, (rowSum));
                }
                
                //System.out.println();
                for (int ln = 0; ln < activeLineages; ln++) { //for all lineages i
                    currState = currTimeLineageStates.get(ln); //stateArray.sampleArray.get(time).get(prevLineageIndexes.get(ln));
                    nextState = nextTimeLineageStates.get(ln); //stateArray.sampleArray.get(time-1).get(ln);
                    transProb *= transProbs.getQuick(currState,nextState);
                    //if (currState == nextState) {
                        //transProb *= Math.exp(-transMatrix.getQuick(currState,currState) * dtTime); //CHANGED FROM Math.exp(-(1-lambda*dt)      
                    //} else {
                        //transProb *= 1 - Math.exp(-transMatrix.getQuick(currState,nextState) * dtTime);
                        //transCount++;
                    //}
                    if (nextState == 0 && currState == 0) {
                        transCount++;
                    }
                }
            }
        }
        if (transProb <= 0.0) {
            transProb = Double.MIN_VALUE;
            //System.out.println("transProb zero at non-coalescent event");
        }
        if (Double.isNaN(transProb)) {
            transProb = Double.MIN_VALUE;
            //System.out.println("transProb NaN at non-coalescent event");
        }
        return transProb;                
    }
    
    public double computeTransCoalProbs(ZVectors dataZ, int zLoc, TreeNode[] tree, DoubleMatrix1D xCurr, ParticleLineageStateArray stateArray, DoubleArrayList theta, int time, double currTime) 
    {
    
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        int lineageMapIndex;
        
        //Get coal nodes at this time
        ArrayList<Integer> coalNodesAtEvent = new ArrayList<Integer>();
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event)==1) {
                int coalNode = dataZ.nodePointers.get(zLoc).get(event);
                coalNodesAtEvent.add(coalNode);
            }
        }

        int coalEvents = coalNodesAtEvent.size();
        
        double transProb = 1.0;
        for (int event = 0; event < coalEvents; event++) {
                
            int coalNode = coalNodesAtEvent.get(event);
            int daughterLineage1Index = tree[coalNode].childNodes[0] - 1;
            int daughterLineage2Index = tree[coalNode].childNodes[1] - 1;
            
            this.updateF(xCurr, theta, currTime);
            this.updateG(xCurr, theta);
            this.updateY(xCurr);

            lineageMapIndex = stateArray.arrayMap.get(time-1).indexOf(daughterLineage1Index);
            int leftDaughterState = stateArray.sampleArray.get(time-1).get(lineageMapIndex);
            lineageMapIndex = stateArray.arrayMap.get(time-1).indexOf(daughterLineage2Index);
            int rightDaughterState = stateArray.sampleArray.get(time-1).get(lineageMapIndex);
            lineageMapIndex = stateArray.arrayMap.get(time-1).indexOf(coalNode);
            int coalNodeState = stateArray.sampleArray.get(time-1).get(lineageMapIndex);
            
            //Probabilistically choose a state for the parent of the two daughter lineages
            double fKtoL = F.get(leftDaughterState, rightDaughterState);
            double fLtoK = F.get(rightDaughterState, leftDaughterState);
            if (leftDaughterState == rightDaughterState) {
                transProb *= 1.0;
                //transCountArray.get(coalNodeState).set(coalNodeState, (transCountArray.get(coalNodeState).get(coalNodeState) + 1));
            } else {
                if (leftDaughterState == coalNodeState) {
                    transProb *= fKtoL / (fKtoL + fLtoK);
                    //transCountArray.get(coalNodeState).set(rightDaughterState, (transCountArray.get(coalNodeState).get(rightDaughterState) + 1));
                }
                if (rightDaughterState == coalNodeState) {
                    transProb *= fLtoK / (fKtoL + fLtoK);
                    //transCountArray.get(coalNodeState).set(leftDaughterState, (transCountArray.get(coalNodeState).get(leftDaughterState) + 1));
                }
                if (leftDaughterState != coalNodeState & rightDaughterState != coalNodeState) {
                    System.out.println("WARNING: Illegitimate coalescent event!!!"); //This seems to only happen when two coalescent events occur at the same event time
                }
            }
            
        }   
        if (transProb <= 0.0) {
            transProb = Double.MIN_VALUE;
            //System.out.println("transProb zero at coalescent event");
        }
        if (Double.isNaN(transProb)) {
            transProb = Double.MIN_VALUE;
            //System.out.println("transProb NaN at coalescent event");
        }
        return transProb;
    }
    
    public double computeGLikeNoCoalInterval(DoubleMatrix1D xCurr, ParticleLineageStateArray stateArray, DoubleArrayList theta, int time, double dtTime, double currTime, ArrayList<Integer> eventTypes) 
    {
        double pGInterval = 1.0;
        int numLineages = stateArray.arrayMap.get(time).size();
        double Ak; double Al; //double rateAkAk; double rateAkAl; double rateAlAl; //double lambdaSumAg; double absDiffLambda; double errTolerance = 0.001;
        ArrayList<Integer> linesAk = new ArrayList<Integer>();
        
        if (numLineages > 1) {
            
            this.updateF(xCurr, theta, currTime);
            this.updateG(xCurr, theta);
            this.updateY(xCurr);
            
            if (eventTypes.contains(1)) { //if a coalescent event occured
                
                //Coalescent event occured so some lineages present in stateArray at time t will NOT be present at time t-1
                
                //Find lineages present at time t and time t-1 and get their states
                ArrayList<Integer> currTimeLineageStates = new ArrayList<Integer>();
                ArrayList<Integer> nextTimeLineageStates = new ArrayList<Integer>();
                int currLineage; int indexOfLineageAtNextTime;
                for (int lin = 0; lin < numLineages; lin++) {
                    currLineage = stateArray.arrayMap.get(time).get(lin);
                    if (stateArray.arrayMap.get(time-1).contains(currLineage)) {
                        indexOfLineageAtNextTime = stateArray.arrayMap.get(time-1).indexOf(currLineage);
                        //currTimeLineageStates.add(stateArray.trueArray.get(time).get(lin));
                        //nextTimeLineageStates.add(stateArray.trueArray.get(time-1).get(indexOfLineageAtNextTime));
                        currTimeLineageStates.add(stateArray.sampleArray.get(time).get(lin));
                        nextTimeLineageStates.add(stateArray.sampleArray.get(time-1).get(indexOfLineageAtNextTime));
                    }
                }
                int activeLineages = currTimeLineageStates.size();
                
                //Compute lineagesInAByState for continuing lineages
                int stateCount;
                for (int st = 0; st < states; st++) {
                    stateCount = 0;
                    for (int lin = 0; lin < activeLineages; lin++) {
                        if (currTimeLineageStates.get(lin) == st) {
                            stateCount++;
                        }
                    }
                    linesAk.add(stateCount);
                }
            } else {
                
                //No coalescent event occured so all lineages present in stateArray at time t will be present at time t-1
                linesAk = stateArray.sampleArrayLineageStates.get(time);
                //linesAk = stateArray.trueArrayLineageStates.get(time);
            }
            
            //Check for incompatible population sizes
            boolean sizeFlag = false;
            for (int st = 0; st < states; st++) {
                if (linesAk.get(st) > Y.getQuick(st)) {
                    sizeFlag = true;
                }
            }

            if (sizeFlag) {
                
                System.out.println("Hit size flag: " + currTime);
                pGInterval = 0.0;
                
            } else {
                
                //General case for m > 1
                double lambdaSum = 0.0;
                for (int k = 0; k < states; k++) {
                    if (Y.getQuick(k) > 0.0) {
                        Ak = linesAk.get(k);
                        if (Ak > 0) {
                            lambdaSum += ((Ak*(Ak-1))/2.0) * (2.0 * F.getQuick(k,k) / (Y.getQuick(k)*Y.getQuick(k)));
                        }
                    }
                }

                for (int k = 0; k < states; k++) {
                    for (int l = 0; l < k; l++) { //l < k so not double counting here
                        if (Y.getQuick(k) > 0.0 & Y.getQuick(l) > 0.0) {
                            Ak = linesAk.get(k);
                            Al = linesAk.get(l);
                            lambdaSum += (Ak * Al) * ((F.getQuick(k,l) + F.getQuick(l,k)) / (Y.getQuick(k)*Y.getQuick(l)));
                        }
                    }
                }
            
                //For specific case of m = 2    
//                Ak = linesAk.get(0);
//                Al = linesAk.get(1);
//                if (Y.getQuick(0) > 0 & Ak > 1) {
//                    rateAkAk = ((Ak*(Ak-1))/2.0) * (2.0 * F.getQuick(0,0) / (Y.getQuick(0)*Y.getQuick(0)));
//                } else {
//                    rateAkAk = 0.0;
//                }
//                if (Y.getQuick(0) > 0 & Y.getQuick(1) > 0) {
//                    rateAkAl = (Ak * Al) * ((F.getQuick(1,0) + F.getQuick(0,1)) / (Y.getQuick(1)*Y.getQuick(0)));
//                } else {
//                    rateAkAl = 0.0;
//                }
//                if (Y.getQuick(1) > 0 & Al > 1) {
//                    rateAlAl = ((Al*(Al-1))/2.0) * (2.0 * F.getQuick(1,1) / (Y.getQuick(1)*Y.getQuick(1)));
//                } else {
//                    rateAlAl = 0.0;
//                }
//                lambdaSum = rateAkAk + rateAkAl + rateAlAl;

                pGInterval = Math.exp(-lambdaSum * dtTime);
            }
            
            if (Double.isNaN(pGInterval)) {
                System.out.println("pGInterval returned NaN");
            }
            if (pGInterval <= 0.0) {
                pGInterval = Double.MIN_VALUE;
                //System.out.println("WARNING: zero weights found");
            }
            //if (pGInterval > 1.0) { //this is a likelihood so its actually ok if this is greater than 1.0
                //pGInterval = 1.0;
                //System.out.println("pGInterval returned greater than one");
            //}
        }
        return pGInterval;
        
    }
    
    public double computeGMarginalLikeNoCoalInterval(DoubleMatrix1D xCurr, LineageStateProbs stateProbs, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime)
    {
        
        //THIS VERSION COMPUTES THE MARGINAL LIKELIHOOD OF THE GENEALOGY WHILE INTEGRATING OVER LINEAGE STATES

        //For each particle
        //int states = coal.Y.size();
        int numLineages = stateProbs.activeLineages.size();
        double lambdaSum; int i; int j; double popCoalRate; double pIIsInK; double pIIsInL; double pJIsInK; double pJIsInL; double linPairCoalRate;
        double Ak; double Al; double rateAkAk; double rateAkAl; double rateAlAl; double lambdaSumAg;
        DoubleArrayList linesAk = new DoubleArrayList();
        
        //int jParticles = 1;
        double pWeight = 1.0;
        
        if (numLineages > 1) {

            //for (int x = 0; x < jParticles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j

                coal.updateF(xCurr, theta, currTime);
                coal.updateG(xCurr, theta);
                coal.updateY(xCurr);
                
                linesAk = stateProbs.getLineagesByState(0); //only one particle
                
                //Check for incompatible population sizes
                boolean sizeFlag = false;
                boolean lineFlag = false;
                for (int st = 0; st < states; st++) {
                    if (linesAk.get(st) > coal.Y.getQuick(st)) {
                        sizeFlag = true;
                    }
                    if (linesAk.get(st) < 1.0) {
                        lineFlag = true;
                    }
                }
                
                if (sizeFlag) {
                    //System.out.println("Yk smaller than Ak");
                    pWeight = 0.0;
                } else {
                
                    lambdaSum = 0.0; //total rate of coalescence across all pairs of lineages

                    if (lineFlag == false)  {
                        
                        //General case for m > 1
                        for (int k = 0; k < states; k++) {
                            if (coal.Y.getQuick(k) > 0.0) {
                                Ak = linesAk.getQuick(k);
                                lambdaSum += ((Ak*(Ak-1))/2) * (2 * coal.F.getQuick(k,k) / (coal.Y.getQuick(k)*coal.Y.getQuick(k)));
                            }
                        }
                        
                        for (int k = 0; k < states; k++) {
                            for (int l = 0; l < k; l++) { //l < k so not double counting here
                                if (coal.Y.getQuick(k) > 0.0 & coal.Y.getQuick(l) > 0.0) {
                                    Ak = linesAk.getQuick(k);
                                    Al = linesAk.getQuick(l);
                                    lambdaSum += (Ak * Al) * ((coal.F.getQuick(k,l) + coal.F.getQuick(l,k)) / (coal.Y.getQuick(k)*coal.Y.getQuick(l)));
                                }
                            }
                        }
                        
                    } else {
                        //lambdaSumAg = Double.NaN;
                        for (int linI = 0; linI < numLineages; linI++) { //for all lineages i

                            i = stateProbs.activeLineages.get(linI);

                            for (int linJ = linI+1; linJ < numLineages; linJ++) { //for all unique pairs of lineages i and j

                                j = stateProbs.activeLineages.get(linJ);

                                //Compute probability of lineages i and j NOT coalescing
                                for (int k = 0; k < states; k++) { //for all states k
                                    if (coal.Y.getQuick(k) <= 0.0) {
                                        //lambdaSum can not increase
                                    } else {
                                        for (int l = 0; l < states; l++) {//for all state l
                                            if (coal.Y.getQuick(l) <= 0.0) {
                                                //lambdaSum can not increase
                                            } else {
                                                //popCoalRate = coal.computeCoalRateOneDirection(k, l);
                                                popCoalRate = coal.F.getQuick(k,l) / (coal.Y.getQuick(k)*coal.Y.getQuick(l));
                                                pIIsInK = stateProbs.matrix.getQuick(0, i, k); //(particle, lineage, state)
                                                pIIsInL = stateProbs.matrix.getQuick(0, i, l);
                                                pJIsInK = stateProbs.matrix.getQuick(0, j, k);
                                                pJIsInL = stateProbs.matrix.getQuick(0, j, l);
                                                linPairCoalRate = popCoalRate * (pIIsInK*pJIsInL + pIIsInL*pJIsInK);
                                                if (Double.isNaN(linPairCoalRate)) { //Checking for NaN's here might not be needed
                                                    linPairCoalRate = 0.0;
                                                }
                                                lambdaSum += linPairCoalRate;

                                            }
                                        }
                                    }
                                }
                            } 
                        }
                    }

                    //System.out.println("LambdaSum = " + lambdaSum + " LambdaSumAg = " + lambdaSumAg);

                    pWeight = Math.exp(-lambdaSum * dtTime);
                }
                if (Double.isNaN(pWeight)) {
                    System.out.println("WARNING: Some particle weights returned NaN");
                    pWeight = Double.MIN_VALUE;
                }
                if (pWeight <= 0.0) {
                    pWeight = Double.MIN_VALUE;
                    //System.out.println("WARNING: Zero or negative particle weights found");
                }
                //if (pWeight > 1.0) {
                    //pWeight = 1.0;
                    //System.out.println("WARNING: Particle weights greater than one found");
                //} 
                //matrix.setQuick(x,time-1,pWeight);
                
            //}
        } else { //less than two lineages
            //for (int x = 0; x < jParticles; x++) {
                //matrix.setQuick(x,time-1,1.0); //assign weight = 1.0 if less than one lineage present
            //}
        }
        
        return pWeight;

    }
    
    public double computeGLikeCoalInterval(ZVectors dataZ, int zLoc, TreeNode[] tree, DoubleMatrix1D xCurr, ParticleLineageStateArray stateArray, DoubleArrayList theta, int time, double currTime) 
    {
        double pGInterval = 1.0;
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        int lineageMapIndex;
        
        //Get coal nodes at this time
        ArrayList<Integer> coalNodesAtEvent = new ArrayList<Integer>();
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event)==1) {
                int coalNode = dataZ.nodePointers.get(zLoc).get(event);
                coalNodesAtEvent.add(coalNode);
            }
        }
        
        this.updateF(xCurr, theta, currTime);
        this.updateG(xCurr, theta);
        this.updateY(xCurr);
        
        int coalEvents = coalNodesAtEvent.size();
        for (int event = 0; event < coalEvents; event++) {
                
            int coalNode = coalNodesAtEvent.get(event);
            int daughterLineage1Index = tree[coalNode].childNodes[0] - 1;
            int daughterLineage2Index = tree[coalNode].childNodes[1] - 1;

            lineageMapIndex = stateArray.arrayMap.get(time-1).indexOf(daughterLineage1Index);
            int daughterLineage1State = stateArray.sampleArray.get(time-1).get(lineageMapIndex);
            lineageMapIndex = stateArray.arrayMap.get(time-1).indexOf(daughterLineage2Index);
            int daughterLineage2State = stateArray.sampleArray.get(time-1).get(lineageMapIndex);
            lineageMapIndex = stateArray.arrayMap.get(time-1).indexOf(coalNode);
            int parentLineageState = stateArray.sampleArray.get(time-1).get(lineageMapIndex);
            
            double probOfEvent = 1.0;
            if (Y.getQuick(daughterLineage1State) == 0.0 || Y.getQuick(daughterLineage2State) == 0.0) {
                probOfEvent = 0.0;
            } else {
                double probParentType = 1.0;
                if (daughterLineage1State == parentLineageState & daughterLineage2State != parentLineageState) {
                    probParentType = F.getQuick(daughterLineage1State,daughterLineage2State) / (F.getQuick(daughterLineage1State,daughterLineage2State) + F.getQuick(daughterLineage2State,daughterLineage1State));
                }
                if (daughterLineage2State == parentLineageState & daughterLineage1State != parentLineageState) {
                    probParentType = F.getQuick(daughterLineage2State,daughterLineage1State) / (F.getQuick(daughterLineage1State,daughterLineage2State) + F.getQuick(daughterLineage2State,daughterLineage1State));
                }
                probOfEvent = probParentType * this.computeCoalRate(daughterLineage1State, daughterLineage2State);

            }
            
            //System.out.println("probParentType" + probParentType);
            //System.out.println("lambda" + lambda);
            
            pGInterval *= probOfEvent;
            if (Double.isNaN(pGInterval)) {
                System.out.println("Some weights returned NaN");
            }
            if (pGInterval <= 0.0) {
                pGInterval = Double.MIN_VALUE;
                    //System.out.println("WARNING: zero weights found");
            }
//            if (pGInterval > 1.0) { //this is a likelihood so its actually ok if this is greater than 1.0
//                pGInterval = 1.0;
//                System.out.println("pGInterval returned greater than one");
//            } 
        }
        return pGInterval;
    }
    
    public double computeGMarginalLikeCoalInterval(ZVectors dataZ, int zLoc, TreeNode[] tree, DoubleMatrix1D xCurr, LineageStateProbs stateProbs, DoubleArrayList theta, int time, StructCoalModel coal, double currTime) 
    {
    
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        //double eventTime = dataZ.absoluteTimes.get(zLoc);
        //int states = coal.Y.size();
        
        //Establish precedence among coalescing lineages
        ArrayList<Integer> coalNodesAtEvent = new ArrayList<Integer>();
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event)==1) {
                int coalNode = dataZ.nodePointers.get(zLoc).get(event);
                coalNodesAtEvent.add(coalNode);
            }
        }
        Collections.sort(coalNodesAtEvent);
        Collections.reverse(coalNodesAtEvent);
        double pWeightNew = 1.0;
        
        int coalEvents = coalNodesAtEvent.size();
        for (int event = 0; event < coalEvents; event++) {
                
            int coalNode = coalNodesAtEvent.get(event);
            int daughterLineage1Index = tree[coalNode].childNodes[0] - 1;
            int daughterLineage2Index = tree[coalNode].childNodes[1] - 1;

            //for (int x = 0; x < jParticles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j

                coal.updateF(xCurr, theta, currTime);
                coal.updateG(xCurr, theta);
                coal.updateY(xCurr);

                double lambdaSum = 0.0;
                DoubleArrayList lambdaSumForEachK = new DoubleArrayList();
                for (int k = 0; k < states; k++) { //for all states k
                    double lambdaSumK = 0.0;
                    if (coal.Y.getQuick(k) <= 0.0) {
                        //lambdaSumK does not increase
                    } else {
                        for (int l = 0; l < states; l++) {//for all state l
                            if (coal.Y.getQuick(l) <= 0.0) {
                                //lambdaSum can not increase
                            } else {
                                double popCoalRate = coal.computeCoalRateOneDirection(k, l);
                                double pIIsInK = stateProbs.matrix.getQuick(0, daughterLineage1Index, k); //(particle, lineage, state)
                                double pIIsInL = stateProbs.matrix.getQuick(0, daughterLineage1Index, l);
                                double pJIsInK = stateProbs.matrix.getQuick(0, daughterLineage2Index, k);
                                double pJIsInL = stateProbs.matrix.getQuick(0, daughterLineage2Index, l);
                                double linPairCoalRate = popCoalRate * (pIIsInK*pJIsInL + pIIsInL*pJIsInK);
                                if (Double.isNaN(linPairCoalRate)) { //checking for NaN's here might not be needed
                                    linPairCoalRate = 0.0;
                                }
                                if (Double.isInfinite(linPairCoalRate)) {
                                    linPairCoalRate = 0.0;
                                }
                                lambdaSum += linPairCoalRate;
                                lambdaSumK += linPairCoalRate;
                            }
                        }
                    }
                    lambdaSumForEachK.add(lambdaSumK);
                }
                
                if (lambdaSum == 0) {
                    //System.out.println("WARNING: Found zero probability of a coalescence event");
                }

                //Update weight to reflect coalescent event}
                pWeightNew *= lambdaSum;
                if (Double.isNaN(pWeightNew)) {
                    System.out.println("WARNING: Some particle weights returned NaN");
                }
                if (pWeightNew <= 0.0) {
                    pWeightNew = Double.MIN_VALUE;
                    //System.out.println("WARNING: Zero  or negative particle weights found");
                }
                //if (pWeightNew > 1.0) {
                    //pWeightNew = 1.0;
                //} 
                //matrix.setQuick(x,time-1,pWeightNew);

                //Update LineageStateProbs
                stateProbs.updateProbsAfterCoal(0, coalNode, daughterLineage1Index, daughterLineage2Index, lambdaSum, lambdaSumForEachK);
            //}

            //Remove daughter lineages from and add parent lineage to activeLineages
            int listIndex1 = stateProbs.activeLineages.indexOf(daughterLineage1Index);
            if (listIndex1 != -1) {
                stateProbs.activeLineages.remove(listIndex1);
            } else {
                System.out.println("Cannot remove daughter lineage");
            }

            int listIndex2 = stateProbs.activeLineages.indexOf(daughterLineage2Index);
            if (listIndex2 != -1) {
                stateProbs.activeLineages.remove(listIndex2);
            } else {
                System.out.println("Cannot remove daughter lineage");
            }

            stateProbs.lineagesRemoved.add(daughterLineage1Index);
            stateProbs.lineagesRemoved.add(daughterLineage2Index);

            stateProbs.activeLineages.add(coalNode);
            stateProbs.lineagesAdded.add(coalNode);     
            
        }
        return pWeightNew;
    }
    
    
    
}
