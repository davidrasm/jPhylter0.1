import cern.colt.list.*;
import cern.colt.matrix.*;
import cern.jet.random.*;
import cern.jet.stat.Descriptive;
import java.util.ArrayList;
/**
 * Run particle filter (with or without resampling)
 * @author David
 */

public class ParticleImportanceSampler
{
    int jParticles; //number of particles
    double margLikelihood;
    StateTrajectory x;
    ParticleAncestry A;
    XTrajectory xTrajSample;
    //LineageStateProbs stateProbs;
    DoubleFactory2D factory2D = DoubleFactory2D.dense;
    DoubleFactory3D factory3D = DoubleFactory3D.dense;
    LineageStateProbsPIS determStateProbs;
    
    public double runPIS(int jNum, EpiModel epi, ZVectors dataZ, TreeNode[] tree, StructCoalModel coal, double startTime, double endTime)
    {
	
        //Run particle filter forward using expected state probs and then resample using importance weights
        
        jParticles = jNum;
        int zLocStart = dataZ.absoluteTimes.indexOf(startTime);
        int zLocEnd = dataZ.absoluteTimes.indexOf(endTime);
        xTrajSample = new XTrajectory(); //holds sampled trajectory
        xTrajSample.zLocStart = zLocStart;
        xTrajSample.zLocEnd = zLocEnd;
        DoubleArrayList xTimes = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zLocStart, zLocEnd);
	int totalTimes = xTimes.size();
        
        DoubleArrayList filterTimes = new DoubleArrayList();
        filterTimes.add(startTime);
        filterTimes.add(endTime);
        int filterTotalTimes = filterTimes.size();
        
        //If starting at specific time
        //double simStartTime = epi.params.get(10);
        //double currTime = startTime;
        //int counter = 0;
        //while (currTime <= simStartTime) {
            //counter++;
            //currTime = dataZ.absoluteTimes.get(counter);
        //}
        //int simStartIndex = counter-1;
        //epi.startTimeIndex = simStartIndex;
        epi.startTimeIndex = 0;
        
	//Set up x matrix
        x = new StateTrajectory();
        x.getMatrix(jParticles, totalTimes, epi.getInitialConditions(), factory3D);
        //OR
        //x.getMatrixWithOffset(jParticles, totalTimes, epi.getInitialConditions(), factory3D, simStartIndex);
        
	//Set up w matrix
        WeightsPIS w = new WeightsPIS();
        w.getMatrix(jParticles, totalTimes, factory2D);
        w.setCumulativeWeights();
        
        //Set up A matrix for particle ancestry - don't need this here since were not resampling
        A = new ParticleAncestry();
        A.getArray(jParticles, totalTimes);
        
        StateTrajectory determX = new StateTrajectory();
        determX.getMatrix(1, totalTimes, epi.getInitialConditions(), factory3D);
        //OR
        //determX.getMatrixWithOffset(jParticles, totalTimes, epi.getInitialConditions(), factory3D, simStartIndex);
        
        /**
         * Step 1: Get deterministic trajectory E(x_1:T)
         */
        DoubleArrayList xSubTimes; DoubleArrayList xDtTimes;
        double timeStartInterval; double timeEndInterval; 	
        int xLocStart; int xLocEnd;
	for (int time = 0; time < (filterTotalTimes - 1); ++time) {
            timeStartInterval = filterTimes.get(time); //absolute start time
            timeEndInterval = filterTimes.get(time+1); //absolute end time
            xLocStart = xTimes.indexOf(timeStartInterval); //index of where to start in xTimes
            xLocEnd = xTimes.indexOf(timeEndInterval); //index of where to end in xTimes
            xSubTimes = (DoubleArrayList) xTimes.partFromTo(xLocStart, xLocEnd);
            xDtTimes = new DoubleArrayList();
            for (int i = 1; i < xSubTimes.size(); ++i) {
                double xDiff = xSubTimes.get(i) - xSubTimes.get(i-1);
		xDtTimes.add(xDiff);
            }
            determX.propogateParticlesDeterministic(epi, xLocStart, xDtTimes, xSubTimes); //Update trajectory
	}
        
        
        
        /**
         * Step 2: Get expected lineage state probs E(s) conditional on E(x_1:T)
         */
        int zLocStartInterval; int zLocEndInterval; int backCounter = 0;
        determStateProbs.activeLineages.removeAll(determStateProbs.activeLineages);
        for (int time = (totalTimes-1); time > 0; --time) {
            
            zLocEndInterval = zLocEnd - backCounter; //zLoc at end of interval (time = t+1)
            zLocStartInterval = zLocEndInterval - 1;
            timeEndInterval = dataZ.absoluteTimes.get(zLocEndInterval); 
            timeStartInterval = dataZ.absoluteTimes.get(zLocStartInterval);
            double dtTime = dataZ.absoluteTimes.get(zLocEndInterval) - dataZ.absoluteTimes.get(zLocStartInterval);
            
            //Get new samples (if samples)
            if (dataZ.omegaEvents.get(zLocEndInterval).contains(2)) { //if this interval begins with a sampling event
                //determStateProbs.addSamples(dataZ, zLocEndInterval, tree, time);
                //OR
                determStateProbs.addSamplesWithPriors(dataZ, zLocEndInterval, tree, time);
            }
            
            //Update lineageStateProbs
            if (determStateProbs.activeLineages.size() > 0) {
                determStateProbs.updateProbs(determX.matrix.viewColumn(time), epi.params, time, coal, dtTime, timeEndInterval);
            }
            
            //Update probs after coalescent events
            if (dataZ.omegaEvents.get(zLocStartInterval).contains(1)) {
                determStateProbs.updateProbsCoalEvent(dataZ, zLocStartInterval, tree, determX.matrix.viewColumn(time-1), epi.params, time, coal, timeStartInterval);
            }
            backCounter++;
        }
        
        
        
        /**
         * Step 3: Run particle to get samples from p(x_1:T | G, E(s_1:T))
         * Have to be careful about how count events at t=0 t=tEnd so they are included in both the forwards and backwards pass
         */
        xSubTimes = new DoubleArrayList();
        xDtTimes = new DoubleArrayList();
        ArrayList<Integer> resamplingTimes = new ArrayList<Integer>();
        
        //resamplingTimes.add(156); resamplingTimes.add(6*52); //for simulated HIV data
        
        //For two-pop model
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 6209.0));
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 6574.0));
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 6939.0));
        
        //For Erik's HIV Trees
        resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 723181)); //Jan-01-1980 
        resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 725008)); //Jan-01-1985
        resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 728660)); //Jan-01-1995
        resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 732313)); //Jan-01-2005
        resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 733408)); //Jan-01-2008
        
        ArrayList<Integer> resampleTimeIndexes = new ArrayList<Integer>();
        resampleTimeIndexes.add(0);
        int resamplingEvents = 0;
	for (int time = 0; time < (totalTimes - 1); ++time) {
            
            timeStartInterval = xTimes.get(time); //absolute start time
            timeEndInterval = xTimes.get(time+1); //absolute end time
            zLocEndInterval = dataZ.absoluteTimes.indexOf(timeEndInterval);
            xSubTimes.add(timeStartInterval); xSubTimes.add(timeEndInterval);
            xDtTimes.add(timeEndInterval - timeStartInterval); //dt has negative values because time is reversed
            
            //Update particle states CHANGE THIS BACK!!!
            x.propogateParticles(epi, time, xDtTimes, xSubTimes); //uses tau leap as is
            
            //Update particle weights
            w.updateWeightsNoCoal(x.matrix.viewColumn(time+1), determStateProbs, epi.params, (time+1), coal, xDtTimes.get(0), timeEndInterval);
            
            //Update weights after coalescent events
            if (dataZ.omegaEvents.get(zLocEndInterval).contains(1)) {
                w.updateWeightsCoal(dataZ, zLocEndInterval, tree, x.matrix.viewColumn(time+1), determStateProbs, epi.params, (time+1), coal, timeEndInterval);
            }
            
            //ESS should be computed from the particles cumulative weights, not just the incremental weights
            w.updateCumulativeWeights(time);
            //double ESS = w.computeEffectiveSampleSize();
            //System.out.println("ESS = " + ESS + " at time: " + timeEndInterval);
            //double meanWeight = Descriptive.mean(w.cumulativeWeights);
            //System.out.println("mean = " + meanWeight);
            //double resamplingESSThreshold = jParticles / 2.0;
            //if (ESS < resamplingESSThreshold) { //Resample here if below resampling ESS threshold
            //if (dataZ.omegaEvents.get(zLocEndInterval).contains(1)) { //resample at coalescent events
            //double meanWeight = Descriptive.mean(w.cumulativeWeights);
            //if (meanWeight < 1e-280) {
            if (resamplingTimes.contains(time)) {
                this.resampleForward(w.cumulativeWeights, time);
                resampleTimeIndexes.add(time+1);
                w.resetCumulativeWeights();
                resamplingEvents++;
            }
            
            xSubTimes.removeAll(xSubTimes); xDtTimes.removeAll(xDtTimes);
	}
        //System.out.println("Number of resampling events: " + resamplingEvents);

        /**
         * Resample particles as part of step 3
         */
        ArrayList<ArrayList<Integer>> ancestryMatrix; // = new ArrayList<ArrayList<Integer>>();
        DoubleArrayList expectedWeights = new DoubleArrayList();
	if (resamplingEvents > 0) {
            if (! resampleTimeIndexes.contains(totalTimes-1)) { //Resample again if didn't resample at last time point
                this.resampleForward(w.cumulativeWeights, (totalTimes-2)); //-1 for indexing and -1 because timeIndex is updated by 1 in resampling method
                resampleTimeIndexes.add(totalTimes-1);
            }
            ancestryMatrix = A.getParticleAncestryMatrix(x.matrix, resampleTimeIndexes); //check if indexing is correct for getSampleForward
            expectedWeights = w.computeExpectedWeights(ancestryMatrix); ////Compute weights for sampled particle trajectories (this is the expected log likelihood)
        } else {
            ancestryMatrix = A.getNullAncestryMatrix(x.matrix); //Make full ancestryMatrix will all particles unique
        }
        ArrayList<ArrayList<Integer>> uniqueParticles = A.uniqueParticles;
        ArrayList<Integer> finalUniqueParticles = uniqueParticles.get(totalTimes-1);
        //System.out.println("Final unique particles: " + finalUniqueParticles);
        ArrayList<Integer> finalParticles = ancestryMatrix.get(totalTimes-1);

        /**
         * Step 4: Re-weight particles based on particle-specific lineage state probs
         */
        w.matrix.assign(1.0); // reset weights
        LineageStateProbsPIS stateProbs = new LineageStateProbsPIS();
        //stateProbs.lineageMap = determStateProbs.lineageMap; stateProbs.absoluteTimes = determStateProbs.absoluteTimes; stateProbs.lineages = determStateProbs.lineages; stateProbs.states = determStateProbs.states;
        stateProbs.states = determStateProbs.states;
        stateProbs.rootNode = determStateProbs.rootNode; stateProbs.particles = jParticles;
        stateProbs.setUpMiniParticleProbArray(finalUniqueParticles);
        //stateProbs.setUpParticleProbArray(uniqueParticles);
        backCounter = 0;
        for (int time = (totalTimes-1); time > 0; --time) {
            
            zLocEndInterval = zLocEnd - backCounter; //zLoc at end of interval (time = t+1)
            zLocStartInterval = zLocEndInterval - 1;
            timeEndInterval = dataZ.absoluteTimes.get(zLocEndInterval); 
            timeStartInterval = dataZ.absoluteTimes.get(zLocStartInterval);
            
            //Update weights after coalescent events
            //if (dataZ.omegaEvents.get(zLocEndInterval).contains(1)) {
                //w.updateParticleWeightsCoal(dataZ, zLocEndInterval, tree, x.matrix.viewColumn(time), stateProbs, epi.params, time, coal, timeEndInterval, ancestryMatrix.get(time), finalParticles, finalUniqueParticles);
                //stateProbs.updateParticleProbsCoalEvent(dataZ, zLocEndInterval, tree, x.matrix.viewColumn(time), epi.params, time, coal, timeEndInterval, ancestryMatrix.get(time), uniqueParticles.get(time));
            //}
            
            //Get new samples (if samples)
            if (dataZ.omegaEvents.get(zLocEndInterval).contains(2)) { //if this interval begins with a sampling event
                //stateProbs.addParticleMiniSamples(dataZ, zLocEndInterval, tree, time);
                //OR
                stateProbs.addParticleMiniSamplesWithPriors(dataZ, zLocEndInterval, tree, time);
                //stateProbs.addParticleSamples(dataZ, zLocEndInterval, tree, time);
            }
            
            //Update weights
            double dtTime = dataZ.absoluteTimes.get(zLocEndInterval) - dataZ.absoluteTimes.get(zLocStartInterval);
            w.updateParticleWeightsNoCoal(x.matrix.viewColumn(time), stateProbs, epi.params, time, coal, dtTime, timeEndInterval, ancestryMatrix.get(time), finalParticles, finalUniqueParticles);
            
            //Update lineage state probs
            stateProbs.updateMiniParticleProbs(x.matrix.viewColumn(time), epi.params, time, coal, dtTime, timeEndInterval, ancestryMatrix.get(time), finalUniqueParticles);
            
            //Update state probs after coalescent events
            if (dataZ.omegaEvents.get(zLocStartInterval).contains(1)) {
                //w.updateWeightsBackCoalFromProbs(dataZ, zLocStartInterval, tree, x.matrix.viewColumn(time-1), stateProbs, epi.params, time, coal, timeStartInterval);
                //w.updateParticleWeightsCoal(dataZ, zLocStartInterval, tree, x.matrix.viewColumn(time-1), stateProbs, epi.params, (time), coal, timeStartInterval, ancestryMatrix.get(time), uniqueParticles.get(time));
                w.updateParticleWeightsCoal(dataZ, zLocStartInterval, tree, x.matrix.viewColumn(time-1), stateProbs, epi.params, time, coal, timeStartInterval, ancestryMatrix.get(time-1), finalParticles, finalUniqueParticles);
                stateProbs.updateMiniParticleProbsCoalEvent(dataZ, zLocStartInterval, tree, x.matrix.viewColumn(time-1), epi.params, (time-1), coal, timeStartInterval, ancestryMatrix.get(time-1), finalUniqueParticles);
            }
            
            //Update back counter
            backCounter++;
        }
        
        //Get (log) importance sampling weights
        DoubleArrayList importanceWeights = new DoubleArrayList(); //unnormalized log weights
        if (resamplingEvents > 0) {
            DoubleArrayList correctedWeights = w.computeCorrectedWeights(); //log likelihood
            for (int p = 0; p < jParticles; p++) {
                importanceWeights.add(correctedWeights.get(p) - expectedWeights.get(p));
            }
            //System.out.println("Difference between expected and corrected weights: ");
            //System.out.println(importanceWeights);
            double covExpCorr = Descriptive.covariance(correctedWeights, expectedWeights);
            //System.out.println("Covariance between expected and corrected weights: ");
            //System.out.println(covExpCorr);
            
        } else {
            importanceWeights = w.computeCorrectedWeights();
            // OR importanceWeights = w.computeExpectedWeights(ancestryMatrix);
        }
        
        //Get weights on linear scale before resampling
        DoubleArrayList expWeights = new DoubleArrayList();
        double maxLogWeights = Descriptive.max(importanceWeights);
        for (int p = 0; p < jParticles; p++) {
            expWeights.add(Math.exp(importanceWeights.get(p) - maxLogWeights)); //recenter before exponentiating to prevent overflow/underflow
        }
        //Normalize
        DoubleArrayList normWeights = new DoubleArrayList();
        double sumExpWeights = Descriptive.sum(expWeights);
        for (int p = 0; p < jParticles; p++) {
            normWeights.add(expWeights.get(p) / sumExpWeights);
        }
          
        margLikelihood = w.computeWeightedMargLikelihood(normWeights);
        if (margLikelihood == Double.POSITIVE_INFINITY) {
            margLikelihood = Double.NEGATIVE_INFINITY;
            System.out.println("Marginal likelihood evaluated as INFINITY");
        }
        
        //For debugging
        //System.out.println("Foward marg likelihood = " + forwardMargLikelihood);
        //System.out.println("Backward marg likelihood = " + margLikelihood);
        //DoubleMatrix2D backwardMatrix = w.matrix.copy();
        //DoubleArrayList diffs = new DoubleArrayList();
        //for (int wTime = 0; wTime < (totalTimes-1); wTime++) {
            //double diff = forwardMatrix.get(0,wTime) - backwardMatrix.get(0,wTime);
            //diffs.add(diff);
            //System.out.println("Time: " + wTime + " error: " + diff);
        //}
        //double sumDiffs = Descriptive.sum(diffs);
        //System.out.println();
        
        //Sample xTraj based on particle weights
        xTrajSample.full = this.sampleParticleTrajectory(normWeights, ancestryMatrix);
        //int randomParticle = (int) Math.round(Math.random() * (jParticles-1));
        //DoubleMatrix2D firstParticle = x.matrix.viewRow(0);
        //xTrajSample.full = x.matrix.viewRow(randomParticle); //check dimensions of returned matrix
        
        //System.out.println("Coal events seen by likelihood: " + w.coalEventsCounted);
        //System.out.println("Move events seen by likelihood: " + w.moveEventsCounted);
        
        return margLikelihood;
    }//End Method
    
    private DoubleMatrix2D sampleParticleTrajectory(DoubleArrayList normWeights, ArrayList<ArrayList<Integer>> ancestryMatrix) 
    {
        
        DoubleArrayList normCumWeights = new DoubleArrayList();
        double cumNormWeight = 0.0;
        for (int j = 0; j < jParticles; ++j) {
            //normWeight = expWeights.getQuick(j) / sumExpWeights;
            cumNormWeight += normWeights.get(j);
            normCumWeights.add(cumNormWeight);
	}
        double uniRand = Math.random();
        double cumVal = 0.0; 
        int counter = 0;
        while (cumVal <= uniRand) {
            cumVal = normCumWeights.getQuick(counter);
            counter++;
        }
        int sample = counter - 1;
        
        int vars = x.matrix.slices();
        int times = x.matrix.columns();
        int pIndex; double state;
        DoubleMatrix2D xSample = factory2D.make(vars, times);
        for (int t = 0; t < times; t++) {
            pIndex = ancestryMatrix.get(t).get(sample);
            for (int v = 0; v < vars; v++) {
                state = x.matrix.getQuick(v, pIndex, t);
                xSample.setQuick(v, t, state);
            }
        }
        return xSample;
        
    }
    
    private void resampleForward(DoubleArrayList cumulativeWeights, int currIndex) 
    {
        //THIS IS THE CONDITIONAL VERSION OF THE ALGORITHM SO RESAMPLING IS ALSO DONE CONDITIONALLY
        //double weightSum = 0.0;
        Integer[] kIndexes = new Integer[jParticles];
        
        //Get weights on linear scale before resampling
        DoubleArrayList expWeights = new DoubleArrayList();
        double maxLogWeights = Descriptive.max(cumulativeWeights);
        for (int p = 0; p < jParticles; p++) {
            expWeights.add(Math.exp(cumulativeWeights.get(p) - maxLogWeights)); //recenter before exponentiating to prevent overflow/underflow
        }
        
        DoubleArrayList normCumWeights = new DoubleArrayList();
        double sumExpWeights = Descriptive.sum(expWeights);
        double normWeight = 0.0;
        double cumNormWeight = 0.0;
        for (int j = 0; j < jParticles; ++j) {
            normWeight = expWeights.getQuick(j) / sumExpWeights;
            cumNormWeight += normWeight;
            normCumWeights.add(cumNormWeight);
	}
        
        //Copy over old xCurr stateVariables (they will be overwrittern upon resampling)
        DoubleMatrix2D xCurrView = x.matrix.viewColumn(currIndex+1);
        DoubleMatrix2D newXCurr = xCurrView.copy();
        
        double uniRand; double cumVal; int counter; int newK;
        for (int n = 0; n < (jParticles); ++n) { //was n < (jParticle - 1) 
            uniRand = Math.random();
            cumVal = 0.0; 
            counter = 0;
            while (cumVal <= uniRand)
            {
                cumVal = normCumWeights.getQuick(counter);
		counter++;
            }
            newK = counter - 1;
            kIndexes[n] = newK; //-1 because counter has already been incremented
            
            //Update state variables
            for (int var = 0; var < newXCurr.rows(); ++var) {
                x.matrix.setQuick(var,n,(currIndex+1), newXCurr.getQuick(var,newK)); //this indexing might be off depending on how viewColumn works on a 3D matrix
            }
	}
        //kIndexes[jParticles-1] = jParticles-1; //CHANGED THIS FOR CONDITIONAL VERSION 
	A.indexes.add(kIndexes);    
    }
    
    public void setUpFilter(ZVectors dataZ, TreeNode[] tree, StructCoalModel coal, double startTime, double endTime)
    {
        
        int zLocStart = dataZ.absoluteTimes.indexOf(startTime);
        int zLocEnd = dataZ.absoluteTimes.indexOf(endTime);
        xTrajSample = new XTrajectory(); //holds sampled trajectory
        xTrajSample.zLocStart = zLocStart;
        xTrajSample.zLocEnd = zLocEnd;
        DoubleArrayList xTimes = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zLocStart, zLocEnd);
        
        determStateProbs = new LineageStateProbsPIS();
        determStateProbs.getLineageMap(xTimes, dataZ, zLocEnd, tree, coal);
        determStateProbs.setUpProbArray();
    }
    
}
