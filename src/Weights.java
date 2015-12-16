import cern.colt.list.*;
import cern.colt.matrix.*;
import cern.jet.stat.Descriptive;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Update particle weights and compute marginal likelihoods
 * @author David
 */

public class Weights {
    
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
            cumulativeWeights.add(0.0); //log weights
        }
    }
    
    public void updateCumulativeWeights(int time)
    {
        //Cumulative weights tracks log weights to avoid numerical underflow
        double newWeight;
        for (int j = 0; j < jParticles; j++) {
            newWeight = cumulativeWeights.getQuick(j) + Math.log(matrix.getQuick(j,time));
            cumulativeWeights.setQuick(j, newWeight);
        }
    }
    
    
    public void resetCumulativeWeights()
    {
        cumulativeWeights.fillFromToWith(0, (jParticles-1), 0.0); //log weights
    }
    
    //This method is not used with the particle Gibbs sampler
    public void updateWeights(DoubleMatrix3D xMatrix, DoubleArrayList theta, int time, ZVectors dataZ, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes, int zLocStart, int xLocStart, CoalModel coal)
    {
        
        DoubleArrayList pVec = new DoubleArrayList();
        for (int j = 0; j < jParticles; j++) {
            pVec.add(1.0);
        }
        //pVec.fillFromToWith(0, jParticles-1, 1.0);
        
        for (int n = 0; n < xDtTimes.size(); n++) {
            
            int zIndex = zLocStart + n;
            int xIndex = xLocStart + n + 1;
            int lines = (int) dataZ.omegaLineages.getQuick(zIndex); // should this be zIndex + 1?
            
            if (lines >= 2) {
            
                double kLineages = (double) lines; // should this be zIndex + 1?
                double kBinom = kLineages * (kLineages - 1) / 2;

                //This likelihood calculation now does account for multiple coalescence events occuring at one time point!!
                int totalEvents = dataZ.omegaEvents.get(zIndex).size();
                int coalEvents = 0;
                for (int event = 0; event < totalEvents; event++) {
                    if (dataZ.omegaEvents.get(zIndex).get(event)==1) {
                        coalEvents++;
                    }
                }      

                double dt = xDtTimes.getQuick(n);
                double timeNow = xSubTimes.get(n);
                DoubleMatrix2D xCurr = xMatrix.viewColumn(xIndex);
                for (int j = 0; j < jParticles; ++j) {

                    //Compute pairwise coal rate \lambda
                    double lambda = coal.computeRate(theta, xCurr.viewColumn(j), timeNow); 

                    double nextP = 0;
                    if (coalEvents >= 1) {
                        nextP = Math.pow(lambda, coalEvents) * Math.exp(-kBinom * lambda * dt);
                    } else {
                        nextP = Math.exp(-kBinom * lambda * dt);
                    }
                    double newP = pVec.getQuick(j) * nextP; 
                    pVec.setQuick(j,newP);
                }
                
            } else {
                //If lines < 2, weights don't change
            }
            
	}
		
	for (int j = 0; j < jParticles; ++j) {
            double pWeight = pVec.getQuick(j);
            if (pWeight <= 0.0 || Double.isNaN(pWeight)) {
                pWeight = Double.MIN_VALUE;
                //System.out.println("WARNING: zero weights found");
            }
            matrix.setQuick(j,time,pWeight);
        }
    }//END method
    
        //This method is not used with the particle Gibbs sampler
    public void updateWeightsSingleStep(DoubleMatrix3D xMatrix, DoubleArrayList theta, int time, ZVectors dataZ, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes, int zLoc, int xLoc, CoalModel coal)
    {
            
        int zIndex = zLoc + 1;
        int xIndex = xLoc + 1;
        int lines = (int) dataZ.omegaLineages.getQuick(zIndex); // should this be zIndex + 1?

        if (lines >= 2) {

            double kLineages = (double) lines; // should this be zIndex + 1?
            double kBinom = kLineages * (kLineages - 1) / 2;

            //This likelihood calculation now does account for multiple coalescence events occuring at one time point!!
            int totalEvents = dataZ.omegaEvents.get(zIndex).size();
            int coalEvents = 0;
            for (int event = 0; event < totalEvents; event++) {
                if (dataZ.omegaEvents.get(zIndex).get(event)==1) {
                    coalEvents++;
                }
            }      

            double dt = xDtTimes.getQuick(0);
            double timeNow = xSubTimes.get(1);
            DoubleMatrix2D xCurr = xMatrix.viewColumn(xIndex);
            double lambda; double pWeight = 0.0;
            for (int j = 0; j < jParticles; ++j) {
                
                //double infecteds = xCurr.getQuick(0, j);
                //if (kLineages > infecteds) {
                    //System.out.println("Too many lineages!!");
                    //pWeight = 0.0;
                //} else {

                //Compute pairwise coal rate \lambda
                lambda = coal.computeRate(theta, xCurr.viewColumn(j), timeNow);

                //Compute likelihood for tree
                if (coalEvents >= 1) {
                    pWeight = Math.pow(lambda, coalEvents) * Math.exp(-kBinom * lambda * dt);
                } else {
                    pWeight = Math.exp(-kBinom * lambda * dt);
                }
                    
                //}

                if (pWeight <= 0.0 || Double.isNaN(pWeight)) {
                    pWeight = Double.MIN_VALUE;
                    //System.out.println("WARNING: zero weights found");
                }
                
                /**
                 * Compute likelihood for tree2
                 */
                
                /**
                 * 
                 * double pFinal = pWeight1 * pWeight2;
                 * 
                 */
                
                matrix.setQuick(j,time,pWeight);

            }

        } else {
            //If lines < 2, weights don't change
        }
		
    }//END method
    
    public void updateWeightsTimeSeries(DoubleMatrix3D xMatrix, int time, int zLoc, ZVectors dataZ, EpiModel epi, TimeSeries tSeries, int zLocStart)
    {
        
        int zIndex = zLoc + 1;
        double currTime = dataZ.absoluteTimes.get(zIndex);
        int timeSeriesIndex = tSeries.times.indexOf(currTime);
        double obsvVal = tSeries.observations.get(timeSeriesIndex);
        
        //Get previous time Series obsveration time
        int prevTimeIndex;
        if (timeSeriesIndex > 0) {
            double prevObsvTime = tSeries.times.get(timeSeriesIndex-1);
            prevTimeIndex = dataZ.absoluteTimes.indexOf(prevObsvTime) - zLocStart; 
        } else {
            prevTimeIndex = 0;
        }

        
        DoubleMatrix2D xCurr = xMatrix.viewColumn(time+1);
        DoubleMatrix2D xPrev = xMatrix.viewColumn(prevTimeIndex); //for computing incidence
        
        double obsvProb; double pWeight;
        for (int j = 0; j < jParticles; ++j) {
        
            obsvProb = epi.getObservationProb(xCurr.viewColumn(j), xPrev.viewColumn(j), obsvVal);
            pWeight = matrix.getQuick(j,time) * obsvProb;
            //System.out.println(pWeight);
            if (pWeight <= 0.0) {
                pWeight = Double.MIN_VALUE;
            }
            matrix.setQuick(j,time,pWeight);
               
        }
        
    }

    public double computeEffectiveSampleSize()
    {
        //Compute ESS as the inverse of the sum of the squared normalized wieghts
        double maxCumWeights = Descriptive.max(cumulativeWeights);
        double sumOfWeights = 0.0;
        for (int j = 0; j < jParticles; j++) {
            sumOfWeights += Math.exp(cumulativeWeights.getQuick(j) - maxCumWeights);    
        }
        double sumOfSquaredNormWeights = 0.0;
        for (int j = 0; j < jParticles; j++) {
            sumOfSquaredNormWeights += Math.pow((cumulativeWeights.getQuick(j) / sumOfWeights), 2.0);    
        }
        double eSS = 1.0 / sumOfSquaredNormWeights;
        return eSS;
        
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
    
    public double computeMargLikelihoodWithAncestry(ArrayList<ArrayList<Integer>> ancestryMatrix)
    {

        int particleIndex; double sumLogWeights;
        double totalLogLike = 0.0;
        for (int t = 0; t < times; t++) {
            sumLogWeights = 0.0;
            for (int p = 0; p < jParticles; p++) {
                particleIndex = ancestryMatrix.get(t+1).get(p);
                sumLogWeights += Math.log(matrix.getQuick(particleIndex, t));
            }
            totalLogLike += sumLogWeights / jParticles;
        }
        return totalLogLike;
        
    }
    
    //These methods are only used in the old StructParticleFilter
    public void updateWeightsBackNoCoalFromProbs(DoubleMatrix2D xCurr, LineageStateProbs stateProbs, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime)
    {
        
        //THIS HAS TO BE GENERALIZED TO WORK WITH MORE THAN TWO STATES (AT LEAST THE BLOCKING METHOD)

        //For each particle
        int states = coal.Y.size();
        int numLineages = stateProbs.activeLineages.size();
        double lambdaSum; int i; int j; double popCoalRate; double pIIsInK; double pIIsInL; double pJIsInK; double pJIsInL; double linPairCoalRate;
        double Ak; double Al; double rateAkAk; double rateAkAl; double rateAlAl; double lambdaSumAg;
        DoubleArrayList linesAk = new DoubleArrayList();
        
        if (numLineages > 1) {

            for (int x = 0; x < jParticles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j

                double pWeight = 1.0;
                
                coal.updateF(xCurr.viewColumn(x), theta, currTime);
                coal.updateG(xCurr.viewColumn(x), theta);
                coal.updateY(xCurr.viewColumn(x));
                
                linesAk = stateProbs.getLineagesByState(x);
                
                //Check for incompatible population sizes
                boolean sizeFlag = false;
                boolean lineFlag = false;
                for (int st = 0; st < states; st++) {
                    if (linesAk.get(st) > coal.Y.getQuick(st) || Double.isNaN(linesAk.get(st))) {
                        sizeFlag = true;
                        //double error = linesAk.get(st) - coal.Y.getQuick(st);
                        //System.out.println("Error at time " + currTime + " is: " + error);
                        //linesAk.set(st, coal.Y.getQuick(st));
                        
                    }
                    if (linesAk.get(st) < 1.0) {
                        lineFlag = false;
                    }
                }
                
                if (sizeFlag) {
                    //System.out.println("Yk smaller than Ak: " + currTime);
                    pWeight = 0.0;
                } else {
                
                    lambdaSum = 0.0; //total rate of coalescence across all pairs of lineages
                    //lambdaSumAg = 0.0;
                    
                    //if (linesAk.getQuick(0) > 1.0) { //override for 1 pop model
                    //if (linesAk.getQuick(0) > 1.0 & linesAk.getQuick(1) > 1.0) {
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
                        
                        
                        //Special case for two states
//                        Ak = linesAk.getQuick(0);
//                        Al = linesAk.getQuick(1);
//
//                        if (coal.Y.getQuick(0) > 0) {
//                            rateAkAk = ((Ak*(Ak-1))/2) * (2 * coal.F.getQuick(0,0) / (coal.Y.getQuick(0)*coal.Y.getQuick(0)));
//                        } else {
//                            rateAkAk = 0.0;
//                        }
//                        if (coal.Y.getQuick(0) > 0 & coal.Y.getQuick(1) > 0) {
//                            rateAkAl = (Ak * Al) * ((coal.F.getQuick(1,0) + coal.F.getQuick(0,1)) / (coal.Y.getQuick(1)*coal.Y.getQuick(0)));
//                        } else {
//                            rateAkAl = 0.0;
//                        }
//                        if (coal.Y.getQuick(1) > 0) {
//                            rateAlAl = ((Al*(Al-1))/2) * (2 * coal.F.getQuick(1,1) / (coal.Y.getQuick(1)*coal.Y.getQuick(1)));
//                        } else {
//                            rateAlAl = 0.0;
//                        }
//                        lambdaSum = rateAkAk + rateAkAl + rateAlAl;
//                        if (Double.isNaN(lambdaSum)) {
//                            System.out.println();
//                        }
                        //System.out.println();
                        
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
                                                pIIsInK = stateProbs.matrix.getQuick(x, i, k); //(particle, lineage, state)
                                                pIIsInL = stateProbs.matrix.getQuick(x, i, l);
                                                pJIsInK = stateProbs.matrix.getQuick(x, j, k);
                                                pJIsInL = stateProbs.matrix.getQuick(x, j, l);
                                                linPairCoalRate = popCoalRate * (pIIsInK*pJIsInL + pIIsInL*pJIsInK);
                                                if (Double.isNaN(linPairCoalRate)) { //Checking for NaN's here might not be needed
                                                    linPairCoalRate = 0.0;
                                                }
                                                lambdaSum += linPairCoalRate;
                                                //if (Double.isNaN(lambdaSum)) {
                                                    //System.out.println();
                                                //}
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
                if (pWeight > 1.0) {
                    //pWeight = 1.0;
                    //System.out.println("WARNING: Particle weights greater than one found");
                } 
                matrix.setQuick(x,time-1,pWeight);
                
            }
        } else { //less than two lineages
            for (int x = 0; x < jParticles; x++) {
                matrix.setQuick(x,time-1,1.0); //assign weight = 1.0 if less than one lineage present
            }
        }

    }
    
    public void updateWeightsBackCoalFromProbs(ZVectors dataZ, int zLoc, TreeNode[] tree, DoubleMatrix2D xCurr, LineageStateProbs stateProbs, DoubleArrayList theta, int time, StructCoalModel coal, double currTime) 
    {
    
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        //double eventTime = dataZ.absoluteTimes.get(zLoc);
        int states = coal.Y.size();
        
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
        
        int coalEvents = coalNodesAtEvent.size();
        for (int event = 0; event < coalEvents; event++) {
                
            int coalNode = coalNodesAtEvent.get(event);
            int daughterLineage1Index = tree[coalNode].childNodes[0] - 1;
            int daughterLineage2Index = tree[coalNode].childNodes[1] - 1;

            for (int x = 0; x < jParticles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j

                coal.updateF(xCurr.viewColumn(x), theta, currTime);
                coal.updateG(xCurr.viewColumn(x), theta);
                coal.updateY(xCurr.viewColumn(x));

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
                                double pIIsInK = stateProbs.matrix.getQuick(x, daughterLineage1Index, k); //(particle, lineage, state)
                                double pIIsInL = stateProbs.matrix.getQuick(x, daughterLineage1Index, l);
                                double pJIsInK = stateProbs.matrix.getQuick(x, daughterLineage2Index, k);
                                double pJIsInL = stateProbs.matrix.getQuick(x, daughterLineage2Index, l);
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
                double pWeightNew = matrix.get(x,time-1) * lambdaSum;
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
                matrix.setQuick(x,time-1,pWeightNew);

                //Update LineageStateProbs
                stateProbs.updateProbsAfterCoal(x, coalNode, daughterLineage1Index, daughterLineage2Index, lambdaSum, lambdaSumForEachK);
            }

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
    }
    
}
