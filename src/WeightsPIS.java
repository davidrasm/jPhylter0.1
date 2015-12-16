import cern.colt.list.*;
import cern.colt.matrix.*;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Update particle weights and compute marginal likelihoods
 * @author David
 */
public class WeightsPIS {
    
    DoubleMatrix2D matrix;
    DoubleMatrix2D cMatrix;
    DoubleArrayList cumulativeWeights;
    DoubleArrayList smoothingWeights;
    int jParticles;
    int times;
    int coalEventsSeen;
    
    //These are here just to speed things up
    DoubleFactory2D factory2D = DoubleFactory2D.dense;
    //DoubleMatrix2D birthRateMatrix; DoubleMatrix2D migrationRateMatrix; DoubleMatrix2D transMatrix;
    
    public void getMatrix(int jNum, int ts, DoubleFactory2D factory2D)
    {
        jParticles = jNum;
        times = ts- 1;
	matrix = factory2D.make(jParticles, (times));
        matrix.assign(1.0);
        coalEventsSeen = 0;
        
    }
    
    public void setCumulativeWeights()
    {
        cumulativeWeights = new DoubleArrayList();
        for (int j = 0; j < jParticles; j++) {
            cumulativeWeights.add(0.0); //log values
        }
    }
    
    public void updateCumulativeWeights(int time)
    {
        //Track cumulative weights on log scale to prevent underflow
        double newWeight;
        for (int j = 0; j < jParticles; j++) {
            newWeight = cumulativeWeights.getQuick(j) + Math.log(matrix.getQuick(j,time));
            cumulativeWeights.setQuick(j, newWeight);
        }
    }
    
    public void resetCumulativeWeights()
    {
        cumulativeWeights.fillFromToWith(0, (jParticles-1), 0.0); //log values
    }
    
    public double computeEffectiveSampleSize()
    {
        //Compute ESS as the inverse of the sum of the squared normalized wieghts
        double sumOfWeights = 0.0;
        for (int j = 0; j < jParticles; j++) {
            sumOfWeights += cumulativeWeights.getQuick(j);    
        }
        double sumOfSquaredNormWeights = 0.0;
        for (int j = 0; j < jParticles; j++) {
            sumOfSquaredNormWeights += Math.pow((cumulativeWeights.getQuick(j) / sumOfWeights), 2.0);    
        }
        double eSS = 1.0 / sumOfSquaredNormWeights;
        return eSS;
        
    }
    
    public void updateWeightsNoCoal(DoubleMatrix2D xCurr, LineageStateProbsPIS stateProbs, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime)
    {

        //For each particle
        int states = coal.Y.size();

        int numLineages = stateProbs.lineageMap.get(time).size();
        double lambdaSum; double popCoalRate; double pIIsInK; double pIIsInL; double pJIsInK; double pJIsInL; double linPairCoalRate;
        double Ak; double Al;
        DoubleArrayList linesAk;
        /**
         * This does not take into consideration that some of the lineages may have coalesced at time t+1 !!!
         */
        linesAk = stateProbs.getLineagesByState(time);
        
        if (numLineages > 1) {

            for (int x = 0; x < jParticles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j

                double pWeight = 1.0;
                
                coal.updateF(xCurr.viewColumn(x), theta, currTime);
                coal.updateG(xCurr.viewColumn(x), theta);
                coal.updateY(xCurr.viewColumn(x));
                
                
                //Check for incompatible population sizes
                boolean sizeFlag = false;
                boolean lineFlag = false;
                for (int st = 0; st < states; st++) {
                    if (linesAk.get(st) > coal.Y.getQuick(st)) {
                        //linesAk.set(st, coal.Y.getQuick(st));
                        sizeFlag = false; //Change this back!
                        
                    }
                    if (linesAk.get(st) < 1.0) {
                        lineFlag = false; //CHANGE THIS BACK!
                    }
                }
                
                if (sizeFlag) {
                    System.out.println("Yk smaller than Ak: " + currTime);
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

                        for (int linI = 0; linI < numLineages; linI++) { //for all lineages i

                            for (int linJ = linI+1; linJ < numLineages; linJ++) { //for all unique pairs of lineages i and j

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
                                                pIIsInK = stateProbs.probs.get(time).get(linI).get(k); //(particle, lineage, state)
                                                pIIsInL = stateProbs.probs.get(time).get(linI).get(l);
                                                pJIsInK = stateProbs.probs.get(time).get(linJ).get(k);
                                                pJIsInL = stateProbs.probs.get(time).get(linJ).get(l);
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

                    pWeight = Math.exp(-lambdaSum * dtTime);
                }
                if (Double.isNaN(pWeight)) {
                    System.out.println("WARNING: Some particle weights returned NaN in step 3");
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
                double newWeight = pWeight * matrix.getQuick(x,time-1);
                //System.out.println("Weight no coal = " + newWeight);
                matrix.setQuick(x,time-1,newWeight);
                
            }
        } //else { //less than two lineages
            //for (int x = 0; x < jParticles; x++) {
                //matrix.setQuick(x,time-1,1.0); //assign weight = 1.0 if less than one lineage present
            //}
        //}

    }
    
    public void updateWeightsCoal(ZVectors dataZ, int zLoc, TreeNode[] tree, DoubleMatrix2D xCurr, LineageStateProbsPIS stateProbs, DoubleArrayList theta, int time, StructCoalModel coal, double currTime) 
    {
    
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        int states = coal.Y.size();
        
        ArrayList<Integer> coalNodesAtEvent = new ArrayList<Integer>();
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event)==1) {
                int coalNode = dataZ.nodePointers.get(zLoc).get(event);
                coalNodesAtEvent.add(coalNode);
            }
        }
        //Not necesary:
        Collections.sort(coalNodesAtEvent);
        Collections.reverse(coalNodesAtEvent);
        
        int coalEvents = coalNodesAtEvent.size();
        double lambdaSum; double popCoalRate; double pIIsInK; double pIIsInL; double pJIsInK; double pJIsInL; double linPairCoalRate; 
        
        for (int event = 0; event < coalEvents; event++) {
                
            int coalNode = coalNodesAtEvent.get(event);
            int daughterLineage1Index = tree[coalNode].childNodes[0] - 1;
            int daughterLineage2Index = tree[coalNode].childNodes[1] - 1;
            int mapIndexDaughter1 = stateProbs.lineageMap.get(time).indexOf(daughterLineage1Index);
            int mapIndexDaughter2 = stateProbs.lineageMap.get(time).indexOf(daughterLineage2Index);

            for (int x = 0; x < jParticles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j

                coal.updateF(xCurr.viewColumn(x), theta, currTime);
                coal.updateG(xCurr.viewColumn(x), theta);
                coal.updateY(xCurr.viewColumn(x));

                lambdaSum = 0.0;
                //DoubleArrayList lambdaSumForEachK = new DoubleArrayList();
                for (int k = 0; k < states; k++) { //for all states k
                    //double lambdaSumK = 0.0;
                    if (coal.Y.getQuick(k) <= 0.0) {
                        //lambdaSumK does not increase
                    } else {
                        for (int l = 0; l < states; l++) {//for all state l
                            if (coal.Y.getQuick(l) <= 0.0) {
                                //lambdaSum can not increase
                            } else {
                                popCoalRate = coal.computeCoalRateOneDirection(k, l);
                                pIIsInK = stateProbs.probs.get(time).get(mapIndexDaughter1).get(k); //(particle, lineage, state)
                                pIIsInL = stateProbs.probs.get(time).get(mapIndexDaughter1).get(l);
                                pJIsInK = stateProbs.probs.get(time).get(mapIndexDaughter2).get(k);
                                pJIsInL = stateProbs.probs.get(time).get(mapIndexDaughter2).get(l);
                                linPairCoalRate = popCoalRate * (pIIsInK*pJIsInL + pIIsInL*pJIsInK);
                                if (Double.isNaN(linPairCoalRate)) { //checking for NaN's here might not be needed
                                    linPairCoalRate = 0.0;
                                }
                                if (Double.isInfinite(linPairCoalRate)) {
                                    linPairCoalRate = 0.0;
                                }
                                lambdaSum += linPairCoalRate;
                                //lambdaSumK += linPairCoalRate;
                            }
                        }
                    }
                    //lambdaSumForEachK.add(lambdaSumK);
                }

                //Update weight to reflect coalescent event}
                double pWeightNew = matrix.get(x,time-1) * lambdaSum;
                if (Double.isNaN(pWeightNew)) {
                    System.out.println("WARNING: Some particle weights returned NaN in step 3 at coal event");
                }
                if (pWeightNew <= 0.0) {
                    pWeightNew = Double.MIN_VALUE;
                    //System.out.println("WARNING: Zero  or negative particle weights found");
                }
                //if (pWeightNew > 1.0) {
                    //pWeightNew = 1.0;
                //}
                //System.out.println("Weight coal = " + pWeightNew);
                matrix.setQuick(x,time-1,pWeightNew);

            }    
            
        }
    }
    
//    public void updateParticleWeightsNoCoal(DoubleMatrix2D xCurr, LineageStateProbsPIS stateProbs, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime)
//    {
//
//        //For each particle
//        int states = coal.Y.size();
//
//        int numLineages = stateProbs.lineageMap.get(time).size();
//        double lambdaSum; double popCoalRate; double pIIsInK; double pIIsInL; double pJIsInK; double pJIsInL; double linPairCoalRate;
//        double Ak; double Al;
//        DoubleArrayList linesAk;
//        
//        if (numLineages > 1) {
//
//            for (int x = 0; x < jParticles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j
//
//                double pWeight = 1.0;
//                
//                 /**
//                 * This does not take into consideration that some of the lineages may have coalesced at time t+1 !!!
//                 */
//                linesAk = stateProbs.getParticleLineagesByState(time, x);
//                
//                coal.updateF(xCurr.viewColumn(x), theta, currTime);
//                coal.updateG(xCurr.viewColumn(x), theta);
//                coal.updateY(xCurr.viewColumn(x));
//                
//                //Check for incompatible population sizes
//                boolean sizeFlag = false;
//                boolean lineFlag = false;
//                for (int st = 0; st < states; st++) {
//                    if (linesAk.get(st) > coal.Y.getQuick(st)) {
//                        sizeFlag = true;
//                        
//                    }
//                    if (linesAk.get(st) < 1.0) {
//                        lineFlag = false; //CHANGE THIS BACK!
//                    }
//                }
//                
//                if (sizeFlag) {
//                    //System.out.println("Yk smaller than Ak: " + currTime);
//                    pWeight = 0.0;
//                    
//                } else {
//                
//                    lambdaSum = 0.0; //total rate of coalescence across all pairs of lineages
//
//                    if (lineFlag == false)  {
//                        
//                        //General case for m > 1
//                        for (int k = 0; k < states; k++) {
//                            if (coal.Y.getQuick(k) > 0.0) {
//                                Ak = linesAk.getQuick(k);
//                                lambdaSum += ((Ak*(Ak-1))/2) * (2 * coal.F.getQuick(k,k) / (coal.Y.getQuick(k)*coal.Y.getQuick(k)));
//                            }
//                        }
//                        
//                        for (int k = 0; k < states; k++) {
//                            for (int l = 0; l < k; l++) { //l < k so not double counting here
//                                if (coal.Y.getQuick(k) > 0.0 & coal.Y.getQuick(l) > 0.0) {
//                                    Ak = linesAk.getQuick(k);
//                                    Al = linesAk.getQuick(l);
//                                    lambdaSum += (Ak * Al) * ((coal.F.getQuick(k,l) + coal.F.getQuick(l,k)) / (coal.Y.getQuick(k)*coal.Y.getQuick(l)));
//                                }
//                            }
//                        }
//                        
//                    } else {
//                        //lambdaSumAg = Double.NaN;
//                        for (int linI = 0; linI < numLineages; linI++) { //for all lineages i
//
//                            //i = stateProbs.activeLineages.get(linI);
//
//                            for (int linJ = linI+1; linJ < numLineages; linJ++) { //for all unique pairs of lineages i and j
//
//                                //j = stateProbs.activeLineages.get(linJ);
//
//                                //Compute probability of lineages i and j NOT coalescing
//                                for (int k = 0; k < states; k++) { //for all states k
//                                    if (coal.Y.getQuick(k) <= 0.0) {
//                                        //lambdaSum can not increase
//                                    } else {
//                                        for (int l = 0; l < states; l++) {//for all state l
//                                            if (coal.Y.getQuick(l) <= 0.0) {
//                                                //lambdaSum can not increase
//                                            } else {
//                                                //popCoalRate = coal.computeCoalRateOneDirection(k, l);
//                                                popCoalRate = coal.F.getQuick(k,l) / (coal.Y.getQuick(k)*coal.Y.getQuick(l));
//                                                pIIsInK = stateProbs.particleProbs.get(x).get(time).get(linI).get(k); //(particle, lineage, state)
//                                                pIIsInL = stateProbs.particleProbs.get(x).get(time).get(linI).get(l);
//                                                pJIsInK = stateProbs.particleProbs.get(x).get(time).get(linJ).get(k);
//                                                pJIsInL = stateProbs.particleProbs.get(x).get(time).get(linJ).get(l);
//                                                linPairCoalRate = popCoalRate * (pIIsInK*pJIsInL + pIIsInL*pJIsInK);
//                                                if (Double.isNaN(linPairCoalRate)) { //Checking for NaN's here might not be needed
//                                                    linPairCoalRate = 0.0;
//                                                }
//                                                lambdaSum += linPairCoalRate;
//                                            }
//                                        }
//                                    }
//                                }
//                            } 
//                        }
//                    }
//
//                    pWeight = Math.exp(-lambdaSum * dtTime);
//                }
//                if (Double.isNaN(pWeight)) {
//                    System.out.println("WARNING: Some particle weights returned NaN");
//                    pWeight = Double.MIN_VALUE;
//                }
//                if (pWeight <= 0.0) {
//                    pWeight = Double.MIN_VALUE;
//                    //System.out.println("WARNING: Zero or negative particle weights found");
//                }
//                if (pWeight > 1.0) {
//                    //pWeight = 1.0;
//                    //System.out.println("WARNING: Particle weights greater than one found");
//                } 
//                matrix.setQuick(x,time-1,pWeight);
//                
//            }
//        } else { //less than two lineages
//            for (int x = 0; x < jParticles; x++) {
//                matrix.setQuick(x,time-1,1.0); //assign weight = 1.0 if less than one lineage present
//            }
//        }
//
//    }
    
    public void updateParticleWeightsNoCoal(DoubleMatrix2D xCurr, LineageStateProbsPIS stateProbs, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime, ArrayList<Integer> ancestors, ArrayList<Integer> finalParticles, ArrayList<Integer> uniqueParticles)
    {

        //For each particle
        int states = coal.Y.size();
        int particlesAtTime = uniqueParticles.size(); int pIndex; int uniqueIndex;
        //int numLineages = stateProbs.lineageMap.get(time).size();
        int numLineages = stateProbs.activeLineages.size();
        double lambdaSum; double popCoalRate; double pIIsInK; double pIIsInL; double pJIsInK; double pJIsInL; double linPairCoalRate;
        double Ak; double Al;
        DoubleArrayList linesAk;
        DoubleArrayList uniqueParticleWeights = new DoubleArrayList();
        
        if (numLineages > 1) {

            for (int x = 0; x < particlesAtTime; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j
                
                //For stateProbs: use x
                //For xMatrix: use particle index
                pIndex = ancestors.get(uniqueParticles.get(x));
                double pWeight = 1.0;
                
                 /**
                 * This does not take into consideration that some of the lineages may have coalesced at time t+1 !!!
                 */
                linesAk = stateProbs.getMiniParticleLineagesByState(x);
                
                coal.updateF(xCurr.viewColumn(pIndex), theta, currTime);
                coal.updateG(xCurr.viewColumn(pIndex), theta);
                coal.updateY(xCurr.viewColumn(pIndex));
                
                //Check for incompatible population sizes
                boolean sizeFlag = false;
                boolean lineFlag = false;
                for (int st = 0; st < states; st++) {
                    if (linesAk.get(st) > coal.Y.getQuick(st)) {
                        sizeFlag = false; //change this back!
                        
                    }
                    if (linesAk.get(st) < 1.0) {
                        lineFlag = false; //CHANGE THIS BACK!
                    }
                }
                
                if (sizeFlag) {
                    //System.out.println("Yk smaller than Ak: " + currTime);
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

                            //i = stateProbs.activeLineages.get(linI);

                            for (int linJ = linI+1; linJ < numLineages; linJ++) { //for all unique pairs of lineages i and j

                                //j = stateProbs.activeLineages.get(linJ);

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
                                                //pIIsInK = stateProbs.particleProbs.get(time).get(x).get(linI).get(k); //(particle, lineage, state)
                                                //pIIsInL = stateProbs.particleProbs.get(time).get(x).get(linI).get(l);
                                                //pJIsInK = stateProbs.particleProbs.get(time).get(x).get(linJ).get(k);
                                                //pJIsInL = stateProbs.particleProbs.get(time).get(x).get(linJ).get(l);
                                                pIIsInK = stateProbs.particleMiniProbs.get(x).get(linI).get(k); //(particle, lineage, state)
                                                pIIsInL = stateProbs.particleMiniProbs.get(x).get(linI).get(l);
                                                pJIsInK = stateProbs.particleMiniProbs.get(x).get(linJ).get(k);
                                                pJIsInL = stateProbs.particleMiniProbs.get(x).get(linJ).get(l);
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

                    pWeight = Math.exp(-lambdaSum * dtTime);
                }
                if (Double.isNaN(pWeight)) {
                    //System.out.println("WARNING: Some particle weights returned NaN in step4");
                    pWeight = Double.MIN_VALUE;
                }
                if (pWeight <= 0.0) {
                    pWeight = Double.MIN_VALUE;
                    //System.out.println("WARNING: Zero or negative particle weights found");
                }
                if (pWeight > 1.0) {
                    //pWeight = 1.0;
                    //System.out.println("WARNING: Particle weights greater than one found");
                }
                //double newWeight = matrix.getQuick
                uniqueParticleWeights.add(pWeight);
                //matrix.setQuick(x,time-1,pWeight);
                
            }
            
            for (int x = 0; x < jParticles; x++) {
                pIndex = finalParticles.get(x);
                uniqueIndex = uniqueParticles.indexOf(pIndex);
                //pIndex = ancestors.get(x);
                //uniqueIndex = uniqueParticles.indexOf(pIndex);
                double newWeight = matrix.getQuick(x,time-1) * uniqueParticleWeights.getQuick(uniqueIndex);
                matrix.setQuick(x,time-1,newWeight);
            }
            
        } //else { //less than two lineages
            //for (int x = 0; x < jParticles; x++) {
                //matrix.setQuick(x,time-1,1.0); //assign weight = 1.0 if less than one lineage present
            //}
        //}

    }
    
//    public void updateParticleWeightsCoal(ZVectors dataZ, int zLoc, TreeNode[] tree, DoubleMatrix2D xCurr, LineageStateProbsPIS stateProbs, DoubleArrayList theta, int time, StructCoalModel coal, double currTime) 
//    {
//    
//        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
//        //double eventTime = dataZ.absoluteTimes.get(zLoc);
//        int states = coal.Y.size();
//        
//        ArrayList<Integer> coalNodesAtEvent = new ArrayList<Integer>();
//        for (int event = 0; event < totalEvents; event++) {
//            if (dataZ.omegaEvents.get(zLoc).get(event)==1) {
//                int coalNode = dataZ.nodePointers.get(zLoc).get(event);
//                coalNodesAtEvent.add(coalNode);
//            }
//        }
//        //Not necesary:
//        Collections.sort(coalNodesAtEvent);
//        Collections.reverse(coalNodesAtEvent);
//        
//        int coalEvents = coalNodesAtEvent.size();
//        for (int event = 0; event < coalEvents; event++) {
//                
//            int coalNode = coalNodesAtEvent.get(event);
//            int daughterLineage1Index = tree[coalNode].childNodes[0] - 1;
//            int daughterLineage2Index = tree[coalNode].childNodes[1] - 1;
//            int mapIndexDaughter1 = stateProbs.lineageMap.get(time).indexOf(daughterLineage1Index);
//            int mapIndexDaughter2 = stateProbs.lineageMap.get(time).indexOf(daughterLineage2Index);
//
//            for (int x = 0; x < jParticles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j
//
//                coal.updateF(xCurr.viewColumn(x), theta, currTime);
//                coal.updateG(xCurr.viewColumn(x), theta);
//                coal.updateY(xCurr.viewColumn(x));
//
//                double lambdaSum = 0.0;
//                //DoubleArrayList lambdaSumForEachK = new DoubleArrayList();
//                for (int k = 0; k < states; k++) { //for all states k
//                    //double lambdaSumK = 0.0;
//                    if (coal.Y.getQuick(k) <= 0.0) {
//                        //lambdaSumK does not increase
//                    } else {
//                        for (int l = 0; l < states; l++) {//for all state l
//                            if (coal.Y.getQuick(l) <= 0.0) {
//                                //lambdaSum can not increase
//                            } else {
//                                double popCoalRate = coal.computeCoalRateOneDirection(k, l);
//                                double pIIsInK = stateProbs.particleProbs.get(x).get(time).get(mapIndexDaughter1).get(k); //(particle, lineage, state)
//                                double pIIsInL = stateProbs.particleProbs.get(x).get(time).get(mapIndexDaughter1).get(l);
//                                double pJIsInK = stateProbs.particleProbs.get(x).get(time).get(mapIndexDaughter2).get(k);
//                                double pJIsInL = stateProbs.particleProbs.get(x).get(time).get(mapIndexDaughter2).get(l);
//                                double linPairCoalRate = popCoalRate * (pIIsInK*pJIsInL + pIIsInL*pJIsInK);
//                                if (Double.isNaN(linPairCoalRate)) { //checking for NaN's here might not be needed
//                                    linPairCoalRate = 0.0;
//                                }
//                                if (Double.isInfinite(linPairCoalRate)) {
//                                    linPairCoalRate = 0.0;
//                                }
//                                lambdaSum += linPairCoalRate;
//                                //lambdaSumK += linPairCoalRate;
//                            }
//                        }
//                    }
//                    //lambdaSumForEachK.add(lambdaSumK);
//                }
//
//                //Update weight to reflect coalescent event}
//                double pWeightNew = matrix.get(x,time-1) * lambdaSum;
//                if (Double.isNaN(pWeightNew)) {
//                    System.out.println("WARNING: Some particle weights returned NaN");
//                }
//                if (pWeightNew <= 0.0) {
//                    pWeightNew = Double.MIN_VALUE;
//                    //System.out.println("WARNING: Zero  or negative particle weights found");
//                }
//                //if (pWeightNew > 1.0) {
//                    //pWeightNew = 1.0;
//                //} 
//                matrix.setQuick(x,time-1,pWeightNew);
//
//            }    
//            
//        }
//    }
    
    public void updateParticleWeightsCoal(ZVectors dataZ, int zLoc, TreeNode[] tree, DoubleMatrix2D xCurr, LineageStateProbsPIS stateProbs, DoubleArrayList theta, int time, StructCoalModel coal, double currTime, ArrayList<Integer> ancestors, ArrayList<Integer> finalParticles, ArrayList<Integer> uniqueParticles) 
    {
    
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        //double eventTime = dataZ.absoluteTimes.get(zLoc);
        int states = coal.Y.size();
        
        ArrayList<Integer> coalNodesAtEvent = new ArrayList<Integer>();
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event)==1) {
                int coalNode = dataZ.nodePointers.get(zLoc).get(event);
                coalNodesAtEvent.add(coalNode);
            }
        }
        //Not necesary:
        Collections.sort(coalNodesAtEvent);
        Collections.reverse(coalNodesAtEvent);
        
        int particlesAtTime = uniqueParticles.size();
        int pIndex; int uniqueIndex;
        DoubleArrayList uniqueParticleWeights = new DoubleArrayList();
        
        int coalEvents = coalNodesAtEvent.size();
        for (int event = 0; event < coalEvents; event++) {
            
            /**
             * Time index for state probs should be time - 1??? No.
             */
            
            int coalNode = coalNodesAtEvent.get(event);
            int daughterLineage1Index = tree[coalNode].childNodes[0] - 1;
            int daughterLineage2Index = tree[coalNode].childNodes[1] - 1;
            //int mapIndexDaughter1 = stateProbs.lineageMap.get(time).indexOf(daughterLineage1Index);
            //int mapIndexDaughter2 = stateProbs.lineageMap.get(time).indexOf(daughterLineage2Index);
            int mapIndexDaughter1 = stateProbs.activeLineages.indexOf(daughterLineage1Index);
            int mapIndexDaughter2 = stateProbs.activeLineages.indexOf(daughterLineage2Index);
            if (mapIndexDaughter1 == -1 || mapIndexDaughter2 == -1) {
                System.out.println("List of activeLineages is off!!");
            }

            for (int x = 0; x < particlesAtTime; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j

                //For stateProbs: use x
                //For xMatrix: use particle index
                pIndex = ancestors.get(uniqueParticles.get(x));
                
                coal.updateF(xCurr.viewColumn(pIndex), theta, currTime);
                coal.updateG(xCurr.viewColumn(pIndex), theta);
                coal.updateY(xCurr.viewColumn(pIndex));

                double lambdaSum = 0.0;
                //DoubleArrayList lambdaSumForEachK = new DoubleArrayList();
                for (int k = 0; k < states; k++) { //for all states k
                    //double lambdaSumK = 0.0;
                    if (coal.Y.getQuick(k) <= 0.0) {
                        //lambdaSumK does not increase
                    } else {
                        for (int l = 0; l < states; l++) {//for all state l
                            if (coal.Y.getQuick(l) <= 0.0) {
                                //lambdaSum can not increase
                            } else {
                                double popCoalRate = coal.computeCoalRateOneDirection(k, l);
                                //double pIIsInK = stateProbs.particleProbs.get(time).get(x).get(mapIndexDaughter1).get(k); //(particle, lineage, state)
                                //double pIIsInL = stateProbs.particleProbs.get(time).get(x).get(mapIndexDaughter1).get(l);
                                //double pJIsInK = stateProbs.particleProbs.get(time).get(x).get(mapIndexDaughter2).get(k);
                                //double pJIsInL = stateProbs.particleProbs.get(time).get(x).get(mapIndexDaughter2).get(l);
                                double pIIsInK = stateProbs.particleMiniProbs.get(x).get(mapIndexDaughter1).get(k); //(particle, lineage, state)
                                double pIIsInL = stateProbs.particleMiniProbs.get(x).get(mapIndexDaughter1).get(l);
                                double pJIsInK = stateProbs.particleMiniProbs.get(x).get(mapIndexDaughter2).get(k);
                                double pJIsInL = stateProbs.particleMiniProbs.get(x).get(mapIndexDaughter2).get(l);
                                double linPairCoalRate = popCoalRate * (pIIsInK*pJIsInL + pIIsInL*pJIsInK);
                                if (Double.isNaN(linPairCoalRate)) { //checking for NaN's here might not be needed
                                    linPairCoalRate = 0.0;
                                }
                                if (Double.isInfinite(linPairCoalRate)) {
                                    linPairCoalRate = 0.0;
                                }
                                lambdaSum += linPairCoalRate;
                                //lambdaSumK += linPairCoalRate;
                            }
                        }
                    }
                    //lambdaSumForEachK.add(lambdaSumK);
                }

                //Update weight to reflect coalescent event}
                double pWeightNew = lambdaSum;
                if (Double.isNaN(pWeightNew)) {
                    System.out.println("WARNING: Some particle weights returned NaN in step 4 at coal event");
                }
                if (pWeightNew <= 0.0) {
                    pWeightNew = Double.MIN_VALUE;
                    //System.out.println("WARNING: Zero  or negative particle weights found");
                }
                //if (pWeightNew > 1.0) {
                    //pWeightNew = 1.0;
                //} 
                //matrix.setQuick(x,time-1,pWeightNew);
                uniqueParticleWeights.add(pWeightNew);

            }
            
            for (int x = 0; x < jParticles; x++) {
                pIndex = finalParticles.get(x);
                uniqueIndex = uniqueParticles.indexOf(pIndex);
                double newWeight = matrix.getQuick(x,time-1) * uniqueParticleWeights.getQuick(uniqueIndex);
                matrix.setQuick(x,time-1,newWeight);
            }
            
        }
    }
    
    public DoubleArrayList computeExpectedWeights(ArrayList<ArrayList<Integer>> ancestryMatrix)
    {
        DoubleArrayList eWeights = new DoubleArrayList();
        int particleIndex;
        for (int p = 0; p < jParticles; p++) {
            double logWeight = 0.0;
            for (int t = 0; t < times; t++) {
                particleIndex = ancestryMatrix.get(t+1).get(p);
                logWeight += Math.log(matrix.getQuick(particleIndex, t));
            }
            eWeights.add(logWeight);
        }
        return eWeights;
        
    }
    
    public DoubleArrayList computeCorrectedWeights()
    {
        DoubleArrayList cWeights = new DoubleArrayList();
        for (int p = 0; p < jParticles; p++) {
            double logWeight = 0.0;
            for (int t = 0; t < times; t++) {
                //int particleIndex = ancestryMatrix.get(t+1).get(p);
                logWeight += Math.log(matrix.getQuick(p, t));
            }
            cWeights.add(logWeight);
        }
        return cWeights;
        
    }
    
    public double computeMargLikelihood()
    {
	double totalLogLike = 0;
        double sumOfWeights;
	double avgWeight;
	for (int n = 0; n < times; ++n) {
            sumOfWeights = 0;
            for (int j = 0; j < jParticles; ++j) {
                sumOfWeights += matrix.getQuick(j,n);
            }
            avgWeight = sumOfWeights / jParticles;
            totalLogLike += Math.log(avgWeight);
	}	
	return totalLogLike;
        
    } //END method
    
    public double computeWeightedMargLikelihood(DoubleArrayList weights)
    {
	double totalLogLike = 0;
        double weightedSum;
	for (int n = 0; n < times; ++n) {
            weightedSum = 0;
            for (int j = 0; j < jParticles; ++j) {
                weightedSum += weights.get(j) * matrix.getQuick(j,n);
            }
            totalLogLike += Math.log(weightedSum);
	}	
	return totalLogLike;
        
    } //END method
    
}
