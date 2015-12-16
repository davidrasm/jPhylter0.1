import cern.colt.list.*;
import cern.colt.matrix.*;
import cern.jet.random.*;
import java.util.ArrayList;
/**
 * Run particle filter (with or without resampling)
 * @author David
 */

public class StructParticleFilter
{
    int jParticles; //number of particles
    double margLikelihood;
    StateTrajectory x;
    ParticleAncestry A;
    XTrajectory xTrajSample;
    LineageEntropies entropies;
    LineageStateArray stateArray;
    LineageStateProbs stateProbs;
    DoubleFactory2D factory2D = DoubleFactory2D.dense;
    DoubleFactory3D factory3D = DoubleFactory3D.dense;
    
    /**
    public void runFilterBackSimulated(int jNum, EpiModel epi, ZVectors dataZ, TreeNode[] tree, StructCoalModel coal, double startTime, double endTime, AbstractDistribution stndNorm, DoubleFactory2D factory2D, DoubleFactory3D factory3D)
    {
	
        //This filter backward simulates the state variable trajectories while simulating the lineage states
        //Allows particle resampling
        
        jParticles = jNum;
        
        int zLocStart = dataZ.absoluteTimes.indexOf(startTime);
        int zLocEnd = dataZ.absoluteTimes.indexOf(endTime);
        xTrajSample = new XTrajectory(); //holds sampled trajectory obtained from particle filter
        xTrajSample.zLocStart = zLocStart;
        xTrajSample.zLocEnd = zLocEnd;
        DoubleArrayList xTimes = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zLocStart, zLocEnd);
	int totalTimes = xTimes.size();
        
	//Set up x matrix
        StateTrajectory x = new StateTrajectory();
        x.getMatrixBack(jParticles, totalTimes, epi.getInitialConditions(), factory3D);
        
	//Set up w matrix
        Weights w = new Weights();
        w.getMatrix(jParticles, totalTimes, factory2D);

        //Set up A matrix for particle ancestry - don't need this here since were not resampling
        //A = new ParticleAncestry();
        //A.getArray(jParticles, totalTimes);
        
        stateArray = new LineageStateArray();
        stateArray.getArray(jParticles, tree.length, coal.Y.size());
        int backCounter = 0; int xLocStart;
        DoubleArrayList xDtTimes = new DoubleArrayList();
        DoubleArrayList xSubTimes = new DoubleArrayList();
        int lastResamplingTimeIndex = totalTimes - 2;
        for (int time = (totalTimes-1); time > 0; --time) {
            
            int zLocEndInterval = zLocEnd - backCounter; //zLoc at end of interval (time = t+1)
            int zLocStartInterval = zLocEndInterval - 1;
            double timeEndInterval = dataZ.absoluteTimes.get(zLocEndInterval); 
            double timeStartInterval = dataZ.absoluteTimes.get(zLocStartInterval);
            
            //Propogate particles back one step
            xSubTimes.add(timeEndInterval); xSubTimes.add(timeStartInterval);
            xDtTimes.add(timeStartInterval - timeEndInterval); //dt has negative values because time is reversed
            x.propogateParticles(epi, time, xDtTimes, xSubTimes);
            xSubTimes.removeAll(xSubTimes); xDtTimes.removeAll(xDtTimes);
            
            //Get new samples (if samples)
            if (dataZ.omegaEvents.get(zLocEndInterval).contains(2)) { //if this interval begins with a sampling event
                stateArray.addSamples(dataZ, zLocEndInterval, tree);
            }
            
            //Update weights
            double dtTime = dataZ.absoluteTimes.get(zLocEndInterval) - dataZ.absoluteTimes.get(zLocStartInterval); //dt has positive values
            //w.updateWeightsBackNoCoalFromStates(x.matrix.viewColumn(time), stateArray, epi.params, time, coal, dtTime, timeEndInterval);
            
            //Update lineageStateProbs
            //stateProbs.updateProbs(x.matrix.viewColumn(time), epi.params, time, coal, dtTime, timeEndInterval);
            stateArray.updateStates(x.matrix.viewColumn(time), epi.params, time, coal, dtTime, timeEndInterval);
            
            //Update weights after coalescent events
            if (dataZ.omegaEvents.get(zLocStartInterval).contains(1)) {
                //w.updateWeightsBackCoalFromStates(dataZ, zLocStartInterval, tree, x.matrix.viewColumn(time-1), stateArray, epi.params, time, coal, timeStartInterval);
                //Resample here if resampling
                this.resampleBackSimulated(x.matrix.viewColumn(time-1), w.matrix, (time-1), lastResamplingTimeIndex); //x and w should also be instance variables
                lastResamplingTimeIndex = time - 1;
            }
            
            //Update back counter
            backCounter++;
        }
        
	margLikelihood = w.computeMargLikelihood();
        if (margLikelihood == Double.POSITIVE_INFINITY) {
            margLikelihood = Double.NEGATIVE_INFINITY;
            System.out.println("Marginal likelihood evaluated as INFINITY");
        }
        
        //If resampling
        //xTrajSample.full = A.getSample(x.matrix, zLocs);
        //Otherwise
        int randomParticle = (int) Math.round(Math.random() * (jParticles-1));
        //DoubleMatrix2D firstParticle = x.matrix.viewRow(0);
        xTrajSample.full = x.matrix.viewRow(randomParticle); //check dimensions of returned matrix
        
        //System.out.println("Coal events seen by likelihood: " + w.coalEventsCounted);
        //System.out.println("Move events seen by likelihood: " + w.moveEventsCounted);
    }//End Method
    */
    
    /**
    public void runFilterBackIntegrated(int jNum, EpiModel epi, ZVectors dataZ, TreeNode[] tree, StructCoalModel coal, double startTime, double endTime, AbstractDistribution stndNorm, DoubleFactory2D factory2D, DoubleFactory3D factory3D)
    {
	
        //This filter backward simulates the state variable trajectories while integrating over the lineage states
        
        jParticles = jNum;
        
        int zLocStart = dataZ.absoluteTimes.indexOf(startTime);
        int zLocEnd = dataZ.absoluteTimes.indexOf(endTime);
        xTrajSample = new XTrajectory(); //holds sampled trajectory obtained from particle filter
        xTrajSample.zLocStart = zLocStart;
        xTrajSample.zLocEnd = zLocEnd;
        DoubleArrayList xTimes = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zLocStart, zLocEnd);
	int totalTimes = xTimes.size();
        
	//Set up x matrix for state trajectories
        x = new StateTrajectory();
        x.getMatrixBack(jParticles, totalTimes, epi.getInitialConditions(), factory3D);
        
	//Set up w matrix for particle weights
        Weights w = new Weights();
        w.getMatrix(jParticles, totalTimes, factory2D);
        
        //Set up A matrix for particle ancestry
        A = new ParticleAncestry();
        A.getArray(jParticles, totalTimes);
        
        stateProbs = new LineageStateProbs();
        stateProbs.getMatrix(jParticles, tree.length, coal.Y.size(), factory3D);
        int backCounter = 0; int xLocStart;
        DoubleArrayList xDtTimes = new DoubleArrayList();
        DoubleArrayList xSubTimes = new DoubleArrayList();
        int lastResamplingTimeIndex = totalTimes - 2; //-2 because weights matrix has one less entry than x matrix
        ArrayList<Integer> resampleTimeIndexes = new ArrayList<Integer>();
        resampleTimeIndexes.add(totalTimes - 1);
        for (int time = (totalTimes-1); time > 0; --time) {
            
            int zLocEndInterval = zLocEnd - backCounter; //zLoc at end of interval (time = t+1)
            int zLocStartInterval = zLocEndInterval - 1;
            double timeEndInterval = dataZ.absoluteTimes.get(zLocEndInterval); 
            double timeStartInterval = dataZ.absoluteTimes.get(zLocStartInterval);
            
            //Propogate particles back one step
            xSubTimes.add(timeEndInterval); xSubTimes.add(timeStartInterval);
            xDtTimes.add(timeStartInterval - timeEndInterval); //dt has negative values because time is reversed
            x.propogateParticles(epi, time, xDtTimes, xSubTimes);
            xSubTimes.removeAll(xSubTimes); xDtTimes.removeAll(xDtTimes);
            
            //Get new samples (if samples)
            if (dataZ.omegaEvents.get(zLocEndInterval).contains(2)) { //if this interval begins with a sampling event
                stateProbs.addSamples(dataZ, zLocEndInterval, tree);
            }
            
            //Update weights
            double dtTime = dataZ.absoluteTimes.get(zLocEndInterval) - dataZ.absoluteTimes.get(zLocStartInterval); //dt has positive values
            w.updateWeightsBackNoCoalFromProbs(x.matrix.viewColumn(time), stateProbs, epi.params, time, coal, dtTime, timeEndInterval);
            
            //Update lineageStateProbs
            stateProbs.updateProbsFast(x.matrix.viewColumn(time), epi.params, time, coal, dtTime, timeEndInterval);
            
            //Update weights after coalescent events
            if (dataZ.omegaEvents.get(zLocStartInterval).contains(1)) {
                w.updateWeightsBackCoalFromProbs(dataZ, zLocStartInterval, tree, x.matrix.viewColumn(time-1), stateProbs, epi.params, time, coal, timeStartInterval);
                //Resample here if resampling
                this.resampleBackIntegrated(w.matrix, (time-1), lastResamplingTimeIndex); //x and w should also be instance variables
                lastResamplingTimeIndex = time - 1;
                resampleTimeIndexes.add(time - 1);
                //System.out.println();
            }
            
            //Update back counter
            backCounter++;
        }
        
	margLikelihood = w.computeMargLikelihood();
        if (margLikelihood == Double.POSITIVE_INFINITY) {
            margLikelihood = Double.NEGATIVE_INFINITY;
            System.out.println("Marginal likelihood evaluated as INFINITY");
        }
        
        //If resampling
        xTrajSample.full = A.getSampleBack(x.matrix, resampleTimeIndexes);
        //Otherwise
        //int randomParticle = (int) Math.round(Math.random() * (jParticles-1));
        //xTrajSample.full = x.matrix.viewRow(randomParticle); //check dimensions of returned matrix
        
        //System.out.println("Coal events seen by likelihood: " + w.coalEventsCounted);
        //System.out.println("Move events seen by likelihood: " + w.moveEventsCounted);
    }//End Method
    */
    
    /*
    public void runFilterForwardBackStochastic(int jNum, EpiModel epi, ZVectors dataZ, TreeNode[] tree, StructCoalModel coal, double startTime, double endTime, AbstractDistribution stndNorm, DoubleFactory2D factory2D, DoubleFactory3D factory3D)
    {
	
        //This filter first forward simulates the particle trajectories and then computes weights on the way back while integrating over lineage states
        jParticles = jNum;
        
        int zLocStart = dataZ.absoluteTimes.indexOf(startTime);
        int zLocEnd = dataZ.absoluteTimes.indexOf(endTime);
        xTrajSample = new XTrajectory(); //holds sampled trajectory obtained from particle filter
        xTrajSample.zLocStart = zLocStart;
        xTrajSample.zLocEnd = zLocEnd;
        DoubleArrayList xTimes = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zLocStart, zLocEnd);
	int totalTimes = xTimes.size();
        
        DoubleArrayList filterTimes = new DoubleArrayList();
        filterTimes.add(startTime);
        filterTimes.add(endTime);
        int filterTotalTimes = filterTimes.size();
        
	//Set up x matrix
        StateTrajectory x = new StateTrajectory();
        x.getMatrix(jParticles, totalTimes, epi.getInitialConditions(), factory3D);
        
	//Set up w matrix
        Weights w = new Weights();
        w.getMatrix(jParticles, totalTimes, factory2D);
        
        //Set up A matrix for particle ancestry - don't need this here since were not resampling
        //A = new ParticleAncestry();
        //A.getArray(jParticles, totalTimes);
        
        //Need these to sample a trajectory from x
        int[] xLocs = new int[totalTimes];
        double timeStart; double timeEnd; 	
        int xLocStart = 0; int xLocEnd = 0;
	
        //All particles simulated together
//        for (int time = 0; time < (filterTotalTimes - 1); ++time) {
//            timeStart = filterTimes.get(time); //absolute start time
//            timeEnd = filterTimes.get(time+1); //absolute end time
//            zLocStart = dataZ.absoluteTimes.indexOf(timeStart);
//            xLocStart = xTimes.indexOf(timeStart); //index of where to start in xTimes
//            xLocs[time] = xLocStart;
//            xLocEnd = xTimes.indexOf(timeEnd); //index of where to end in xTimes
//            DoubleArrayList xSubTimes = (DoubleArrayList) xTimes.partFromTo(xLocStart, xLocEnd);
//            DoubleArrayList xDtTimes = new DoubleArrayList();
//            for (int i = 1; i < xSubTimes.size(); ++i) {
//                double xDiff = xSubTimes.get(i) - xSubTimes.get(i-1);
//		xDtTimes.add(xDiff);
//            }
//
//            //Update particle states
//            x.propogateParticles(epi, xLocStart, xDtTimes, xSubTimes);
//	}
//        xLocs[totalTimes-1] = xLocEnd;
        
        
        //Simulate each particle individually
        timeStart = filterTimes.get(0); //absolute start time
        timeEnd = filterTimes.get(1); //absolute end time
        zLocStart = dataZ.absoluteTimes.indexOf(timeStart);
        xLocStart = xTimes.indexOf(timeStart); //index of where to start in xTimes
        xLocs[0] = xLocStart;
        xLocEnd = xTimes.indexOf(timeEnd); //index of where to end in xTimes
        DoubleArrayList xSubTimes = (DoubleArrayList) xTimes.partFromTo(xLocStart, xLocEnd);
        DoubleArrayList xDtTimes = new DoubleArrayList();
        for (int i = 1; i < xSubTimes.size(); ++i) {
            double xDiff = xSubTimes.get(i) - xSubTimes.get(i-1);
            xDtTimes.add(xDiff);
        }
        boolean extinction;
        for (int particle = 0; particle < jParticles; particle++) {
            //Make sure pop doesn't go extinct
            extinction = true;
            while (extinction) {
                x.propogateSingleParticle(epi, xLocStart, xDtTimes, xSubTimes, particle);
                if (epi.alive == true) {
                    extinction = false;
                }     
            }
        }
        xLocs[totalTimes-1] = xLocEnd;
        
        LineageStateProbs stateProbs = new LineageStateProbs();
        stateProbs.getMatrix(jParticles, tree.length, coal.Y.size(), factory3D);
        
        int backCounter = 0;
        for (int time = (totalTimes-1); time > 0; --time) {
            
            int zLocEndInterval = zLocEnd - backCounter; //zLoc at end of interval (time = t+1)
            int zLocStartInterval = zLocEndInterval - 1;
            double timeEndInterval = dataZ.absoluteTimes.get(zLocEndInterval); 
            double timeStartInterval = dataZ.absoluteTimes.get(zLocStartInterval);
            
            //Get new samples (if samples)
            if (dataZ.omegaEvents.get(zLocEndInterval).contains(2)) { //if this interval begins with a sampling event
                stateProbs.addSamples(dataZ, zLocEndInterval, tree);
            }
            
            //Update weights
            double dtTime = dataZ.absoluteTimes.get(zLocEndInterval) - dataZ.absoluteTimes.get(zLocStartInterval);
            w.updateWeightsBackNoCoalFromProbs(x.matrix.viewColumn(time), stateProbs, epi.params, time, coal, dtTime, timeEndInterval);
            
            //Update lineageStateProbs
            stateProbs.updateProbs(x.matrix.viewColumn(time), epi.params, time, coal, dtTime, timeEndInterval);
            
            //Update weights after coalescent events
            if (dataZ.omegaEvents.get(zLocStartInterval).contains(1)) {
                w.updateWeightsBackCoalFromProbs(dataZ, zLocStartInterval, tree, x.matrix.viewColumn(time-1), stateProbs, epi.params, time, coal, timeStartInterval);
            }
            
            //Update back counter
            backCounter++;
        }
        
	margLikelihood = w.computeMargLikelihood();
        if (margLikelihood == Double.POSITIVE_INFINITY) {
            margLikelihood = Double.NEGATIVE_INFINITY;
            System.out.println("Marginal likelihood evaluated as INFINITY");
        }
        
        //If resampling
        //xTrajSample.full = A.getSample(x.matrix, zLocs);
        //Otherwise
        int randomParticle = (int) Math.round(Math.random() * (jParticles-1));
        //DoubleMatrix2D firstParticle = x.matrix.viewRow(0);
        xTrajSample.full = x.matrix.viewRow(randomParticle); //check dimensions of returned matrix
        
        //System.out.println("Coal events seen by likelihood: " + w.coalEventsCounted);
        //System.out.println("Move events seen by likelihood: " + w.moveEventsCounted);
    }//End Method
    * */
    
    public void runFilterForwardBack(int jNum, EpiModel epi, ZVectors dataZ, TreeNode[] tree, StructCoalModel coal, double startTime, double endTime, AbstractDistribution stndNorm, DoubleFactory2D factory2D, DoubleFactory3D factory3D)
    {
	
        //This filter first forward simulates the particle trajectories and then computes weights on the way back while integrating over lineage states
        
        jParticles = jNum;
        
        int zLocStart = dataZ.absoluteTimes.indexOf(startTime);
        int zLocEnd = dataZ.absoluteTimes.indexOf(endTime);
        xTrajSample = new XTrajectory(); //holds sampled trajectory obtained from particle filter
        xTrajSample.zLocStart = zLocStart;
        xTrajSample.zLocEnd = zLocEnd;
        DoubleArrayList xTimes = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zLocStart, zLocEnd);
	int totalTimes = xTimes.size();
        
        DoubleArrayList filterTimes = new DoubleArrayList();
        filterTimes.add(startTime);
        filterTimes.add(endTime);
        int filterTotalTimes = filterTimes.size();
        
	//Set up x matrix
        StateTrajectory x = new StateTrajectory();
        x.getMatrix(jParticles, totalTimes, epi.getInitialConditions(), factory3D);
        
	//Set up w matrix
        Weights w = new Weights();
        w.getMatrix(jParticles, totalTimes, factory2D);
        
        //Set up A matrix for particle ancestry - don't need this here since were not resampling
        //A = new ParticleAncestry();
        //A.getArray(jParticles, totalTimes);
        
        //Need these to sample a trajectory from x
        int[] xLocs = new int[totalTimes];
        double timeStart; double timeEnd; 	
        int xLocStart = 0; int xLocEnd = 0;
	for (int time = 0; time < (filterTotalTimes - 1); ++time) {
            timeStart = filterTimes.get(time); //absolute start time
            timeEnd = filterTimes.get(time+1); //absolute end time
            zLocStart = dataZ.absoluteTimes.indexOf(timeStart);
            xLocStart = xTimes.indexOf(timeStart); //index of where to start in xTimes
            xLocs[time] = xLocStart;
            xLocEnd = xTimes.indexOf(timeEnd); //index of where to end in xTimes
            DoubleArrayList xSubTimes = (DoubleArrayList) xTimes.partFromTo(xLocStart, xLocEnd);
            DoubleArrayList xDtTimes = new DoubleArrayList();
            for (int i = 1; i < xSubTimes.size(); ++i) {
                double xDiff = xSubTimes.get(i) - xSubTimes.get(i-1);
		xDtTimes.add(xDiff);
            }

            //Update particle states
            x.propogateParticlesDeterministic(epi, xLocStart, xDtTimes, xSubTimes);
	}
        xLocs[totalTimes-1] = xLocEnd;
        
        LineageStateProbs stateProbs = new LineageStateProbs();
        stateProbs.getMatrix(jParticles, tree.length, coal.Y.size(), factory3D);
        
        int backCounter = 0;
        for (int time = (totalTimes-1); time > 0; --time) {
            
            int zLocEndInterval = zLocEnd - backCounter; //zLoc at end of interval (time = t+1)
            int zLocStartInterval = zLocEndInterval - 1;
            double timeEndInterval = dataZ.absoluteTimes.get(zLocEndInterval); 
            double timeStartInterval = dataZ.absoluteTimes.get(zLocStartInterval);
            
            //Get new samples (if samples)
            if (dataZ.omegaEvents.get(zLocEndInterval).contains(2)) { //if this interval begins with a sampling event
                //stateProbs.addSamples(dataZ, zLocEndInterval, tree);
                stateProbs.addSamples(dataZ, zLocEndInterval, tree);
            }
            
            //Update weights
            double dtTime = dataZ.absoluteTimes.get(zLocEndInterval) - dataZ.absoluteTimes.get(zLocStartInterval);
            w.updateWeightsBackNoCoalFromProbs(x.matrix.viewColumn(time), stateProbs, epi.params, time, coal, dtTime, timeEndInterval);
            
            //Update lineageStateProbs
            stateProbs.updateProbs(x.matrix.viewColumn(time), epi.params, time, coal, dtTime, timeEndInterval);
            
            //Update weights after coalescent events
            if (dataZ.omegaEvents.get(zLocStartInterval).contains(1)) {
                w.updateWeightsBackCoalFromProbs(dataZ, zLocStartInterval, tree, x.matrix.viewColumn(time-1), stateProbs, epi.params, time, coal, timeStartInterval);
            }
            
            //Get prob
            //double currProb = w.matrix.getQuick(0, time-1);
            //System.out.println("Curr prob at time" + timeEndInterval + " = " + currProb);
            
            //Update back counter
            backCounter++;
        }
        
	margLikelihood = w.computeMargLikelihood();
        if (margLikelihood == Double.POSITIVE_INFINITY) {
            margLikelihood = Double.NEGATIVE_INFINITY;
            System.out.println("Marginal likelihood evaluated as INFINITY");
        }
        
        //If resampling
        //xTrajSample.full = A.getSample(x.matrix, zLocs);
        //Otherwise
        //int randomParticle = (int) Math.round(Math.random() * (jParticles-1));
        //DoubleMatrix2D firstParticle = x.matrix.viewRow(0);
        xTrajSample.full = x.matrix.viewRow(0); //check dimensions of returned matrix
        
        //System.out.println("Coal events seen by likelihood: " + w.coalEventsCounted);
        //System.out.println("Move events seen by likelihood: " + w.moveEventsCounted);
    }//End Method
    
    public double runFilterForwardBackMMH(int jNum, EpiModel epi, ZVectors dataZ, TreeNode[] tree, StructCoalModel coal, double startTime, double endTime)
    {
	
        //This filter first forward simulates the particle trajectories and then computes weights on the way back while integrating over lineage states
        
        jParticles = jNum;
        
        int zLocStart = dataZ.absoluteTimes.indexOf(startTime);
        int zLocEnd = dataZ.absoluteTimes.indexOf(endTime);
        xTrajSample = new XTrajectory(); //holds sampled trajectory obtained from particle filter
        xTrajSample.zLocStart = zLocStart;
        xTrajSample.zLocEnd = zLocEnd;
        DoubleArrayList xTimes = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zLocStart, zLocEnd);
	int totalTimes = xTimes.size();
        
        DoubleArrayList filterTimes = new DoubleArrayList();
        filterTimes.add(startTime);
        filterTimes.add(endTime);
        int filterTotalTimes = filterTimes.size();
        
	//Set up x matrix
        x = new StateTrajectory();
        x.getMatrix(jParticles, totalTimes, epi.getInitialConditions(), factory3D);
        
	//Set up w matrix
        Weights w = new Weights();
        w.getMatrix(jParticles, totalTimes, factory2D);
        
        //Set up A matrix for particle ancestry - don't need this here since were not resampling
        //A = new ParticleAncestry();
        //A.getArray(jParticles, totalTimes);
        
        //Need these to sample a trajectory from x
        int[] xLocs = new int[totalTimes];
        double timeStart; double timeEnd; 	
        int xLocStart = 0; int xLocEnd = 0;
	for (int time = 0; time < (filterTotalTimes - 1); ++time) {
            timeStart = filterTimes.get(time); //absolute start time
            timeEnd = filterTimes.get(time+1); //absolute end time
            zLocStart = dataZ.absoluteTimes.indexOf(timeStart);
            xLocStart = xTimes.indexOf(timeStart); //index of where to start in xTimes
            xLocs[time] = xLocStart;
            xLocEnd = xTimes.indexOf(timeEnd); //index of where to end in xTimes
            DoubleArrayList xSubTimes = (DoubleArrayList) xTimes.partFromTo(xLocStart, xLocEnd);
            DoubleArrayList xDtTimes = new DoubleArrayList();
            for (int i = 1; i < xSubTimes.size(); ++i) {
                double xDiff = xSubTimes.get(i) - xSubTimes.get(i-1);
		xDtTimes.add(xDiff);
            }

            //Update particle states
            x.propogateParticles(epi, xLocStart, xDtTimes, xSubTimes);
	}
        xLocs[totalTimes-1] = xLocEnd;
        
        LineageStateProbs stateProbs = new LineageStateProbs();
        stateProbs.getMatrix(jParticles, tree.length, coal.Y.size(), factory3D);
        
        int backCounter = 0;
        for (int time = (totalTimes-1); time > 0; --time) {
            
            int zLocEndInterval = zLocEnd - backCounter; //zLoc at end of interval (time = t+1)
            int zLocStartInterval = zLocEndInterval - 1;
            double timeEndInterval = dataZ.absoluteTimes.get(zLocEndInterval); 
            double timeStartInterval = dataZ.absoluteTimes.get(zLocStartInterval);
            
            //Get new samples (if samples)
            if (dataZ.omegaEvents.get(zLocEndInterval).contains(2)) { //if this interval begins with a sampling event
                //stateProbs.addSamples(dataZ, zLocEndInterval, tree);
                stateProbs.addSamples(dataZ, zLocEndInterval, tree);
            }
            
            //Update weights
            double dtTime = dataZ.absoluteTimes.get(zLocEndInterval) - dataZ.absoluteTimes.get(zLocStartInterval);
            w.updateWeightsBackNoCoalFromProbs(x.matrix.viewColumn(time), stateProbs, epi.params, time, coal, dtTime, timeEndInterval);
            
            //Update lineageStateProbs
            stateProbs.updateProbs(x.matrix.viewColumn(time), epi.params, time, coal, dtTime, timeEndInterval);
            
            //Update weights after coalescent events
            if (dataZ.omegaEvents.get(zLocStartInterval).contains(1)) {
                w.updateWeightsBackCoalFromProbs(dataZ, zLocStartInterval, tree, x.matrix.viewColumn(time-1), stateProbs, epi.params, time, coal, timeStartInterval);
            }
            
            //Get prob
            //double currProb = w.matrix.getQuick(0, time-1);
            //System.out.println("Curr prob at time" + timeEndInterval + " = " + currProb);
            
            //Update back counter
            backCounter++;
        }
        
	margLikelihood = w.computeMargLikelihood();
        if (margLikelihood == Double.POSITIVE_INFINITY) {
            margLikelihood = Double.NEGATIVE_INFINITY;
            System.out.println("Marginal likelihood evaluated as INFINITY");
        }
        
        //If resampling
        //xTrajSample.full = A.getSample(x.matrix, zLocs);
        //Otherwise
        int randomParticle = (int) Math.round(Math.random() * (jParticles-1));
        xTrajSample.full = x.matrix.viewRow(randomParticle); //check dimensions of returned matrix
        
        //System.out.println("Coal events seen by likelihood: " + w.coalEventsCounted);
        //System.out.println("Move events seen by likelihood: " + w.moveEventsCounted);
        
        return margLikelihood;
    }//End Method
    
    public void recomputeLikelihood(int jNum, DoubleArrayList theta, ZVectors dataZ, DoubleArrayList xInit, TreeNode[] tree, StructCoalModel coal, double startTime, double endTime, AbstractDistribution stndNorm, DoubleFactory2D factory2D, DoubleFactory3D factory3D)
    {
	
        //Recompute marginal likelihoods without resimulating from the model
        
        jParticles = jNum;
        
        int zLocStart = dataZ.absoluteTimes.indexOf(startTime);
        int zLocEnd = dataZ.absoluteTimes.indexOf(endTime);
        DoubleArrayList xTimes = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zLocStart, zLocEnd);
	int totalTimes = xTimes.size();
        
        DoubleArrayList filterTimes = new DoubleArrayList();
        filterTimes.add(startTime);
        filterTimes.add(endTime);
        int filterTotalTimes = filterTimes.size();
        
	//Set up x matrix
        StateTrajectory x = new StateTrajectory();
        x.getMatrix(1, totalTimes, xInit, factory3D);
        for (int sliceI = 0; sliceI < xInit.size(); sliceI++) {
            for (int colJ = 0; colJ < totalTimes; colJ++) {
                x.matrix.setQuick(sliceI, 0, colJ, xTrajSample.full.getQuick(sliceI, colJ));
            }
        }
        
	//Set up w matrix
        Weights w = new Weights();
        w.getMatrix(jParticles, totalTimes, factory2D);
        
        LineageStateProbs stateProbs = new LineageStateProbs();
        stateProbs.getMatrix(jParticles, tree.length, coal.Y.size(), factory3D);
        
        int backCounter = 0;
        for (int time = (totalTimes-1); time > 0; --time) {
            
            int zLocEndInterval = zLocEnd - backCounter; //zLoc at end of interval (time = t+1)
            int zLocStartInterval = zLocEndInterval - 1;
            double timeEndInterval = dataZ.absoluteTimes.get(zLocEndInterval); 
            double timeStartInterval = dataZ.absoluteTimes.get(zLocStartInterval);
            
            //Get new samples (if samples)
            if (dataZ.omegaEvents.get(zLocEndInterval).contains(2)) { //if this interval begins with a sampling event
                stateProbs.addSamples(dataZ, zLocEndInterval, tree); //can get rid of lineageMap since just need sampling states
            }
            
            //Update weights
            double dtTime = dataZ.absoluteTimes.get(zLocEndInterval) - dataZ.absoluteTimes.get(zLocStartInterval);
            w.updateWeightsBackNoCoalFromProbs(x.matrix.viewColumn(time), stateProbs, theta, time, coal, dtTime, timeEndInterval);
            
            //Update lineageStateProbs
            stateProbs.updateProbs(x.matrix.viewColumn(time), theta, time, coal, dtTime, timeEndInterval);
            
            //Update weights after coalescent events
            if (dataZ.omegaEvents.get(zLocStartInterval).contains(1)) {
                w.updateWeightsBackCoalFromProbs(dataZ, zLocStartInterval, tree, x.matrix.viewColumn(time-1), stateProbs, theta, time, coal, timeStartInterval);
            }
            
            //Update back counter
            backCounter++;
        }
        
	margLikelihood = w.computeMargLikelihood();
        if (margLikelihood == Double.POSITIVE_INFINITY) {
            margLikelihood = Double.NEGATIVE_INFINITY;
            System.out.println("Marginal likelihood evaluated as INFINITY");
        }
        
    }//End Method
    
    /*
    public void runFilterBackEntropy(int jNum, EpiModel epi, ZVectors dataZ, TreeNode[] tree, StructCoalModel coal, double startTime, double endTime, AbstractDistribution stndNorm, DoubleFactory2D factory2D, DoubleFactory3D factory3D)
    {
	
        //This filter backward simulates the state variable trajectories while integrating over the lineage states
        
        jParticles = jNum;
        
        int zLocStart = dataZ.absoluteTimes.indexOf(startTime);
        int zLocEnd = dataZ.absoluteTimes.indexOf(endTime);
        xTrajSample = new XTrajectory(); //holds sampled trajectory obtained from particle filter
        xTrajSample.zLocStart = zLocStart;
        xTrajSample.zLocEnd = zLocEnd;
        DoubleArrayList xTimes = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zLocStart, zLocEnd);
	int totalTimes = xTimes.size();
        
	//Set up x matrix for state trajectories
        x = new StateTrajectory();
        x.getMatrixBack(jParticles, totalTimes, epi.getInitialConditions(), factory3D);
        
	//Set up w matrix for particle weights
        Weights w = new Weights();
        w.getMatrix(jParticles, totalTimes, factory2D);
        
        //Set up A matrix for particle ancestry
        A = new ParticleAncestry();
        A.getArray(jParticles, totalTimes);
        
        //Set up stateProbs matrix
        stateProbs = new LineageStateProbs();
        stateProbs.getMatrix(jParticles, tree.length, coal.Y.size(), factory3D);
        
        //Set up matrix for lineage entropies
        entropies = new LineageEntropies();
        entropies.getArrays(tree.length, coal.Y.size());
        
        int backCounter = 0; int xLocStart;
        DoubleArrayList xDtTimes = new DoubleArrayList();
        DoubleArrayList xSubTimes = new DoubleArrayList();
        int lastResamplingTimeIndex = totalTimes - 2; //-2 because weights matrix has one less entry than x matrix
        ArrayList<Integer> resampleTimeIndexes = new ArrayList<Integer>();
        resampleTimeIndexes.add(totalTimes - 1);
        for (int time = (totalTimes-1); time > 0; --time) {
            
            int zLocEndInterval = zLocEnd - backCounter; //zLoc at end of interval (time = t+1)
            int zLocStartInterval = zLocEndInterval - 1;
            double timeEndInterval = dataZ.absoluteTimes.get(zLocEndInterval); 
            double timeStartInterval = dataZ.absoluteTimes.get(zLocStartInterval);
            
            //Propogate particles back one step
            xSubTimes.add(timeEndInterval); xSubTimes.add(timeStartInterval);
            xDtTimes.add(timeStartInterval - timeEndInterval); //dt has negative values because time is reversed
            x.propogateParticles(epi, time, xDtTimes, xSubTimes);
            xSubTimes.removeAll(xSubTimes); xDtTimes.removeAll(xDtTimes);
            
            //Get new samples (if samples)
            if (dataZ.omegaEvents.get(zLocEndInterval).contains(2)) { //if this interval begins with a sampling event
                stateProbs.addSamples(dataZ, zLocEndInterval, tree);
                entropies.addSamples(dataZ, zLocEndInterval, tree);
            }
            
            //Update weights
            double dtTime = dataZ.absoluteTimes.get(zLocEndInterval) - dataZ.absoluteTimes.get(zLocStartInterval); //dt has positive values
            w.updateWeightsBackNoCoalFromProbs(x.matrix.viewColumn(time), stateProbs, epi.params, time, coal, dtTime, timeEndInterval);
            
            //Update lineageStateProbs
            stateProbs.updateProbsFast(x.matrix.viewColumn(time), epi.params, time, coal, dtTime, timeEndInterval);
            entropies.updateEntropies(dataZ.absoluteTimes.get(zLocStartInterval), stateProbs);
            
            //Update weights after coalescent events
            if (dataZ.omegaEvents.get(zLocStartInterval).contains(1)) {
                w.updateWeightsBackCoalFromProbs(dataZ, zLocStartInterval, tree, x.matrix.viewColumn(time-1), stateProbs, epi.params, time, coal, timeStartInterval);
                entropies.updateEntropiesCoal(dataZ, zLocStartInterval, tree, stateProbs);
                //Resample here if resampling
                //this.resampleBackIntegrated(w.matrix, (time-1), lastResamplingTimeIndex); //x and w should also be instance variables
                //lastResamplingTimeIndex = time - 1;
                //resampleTimeIndexes.add(time - 1);
                //System.out.println();
            }
            
            //Update back counter
            backCounter++;
        }
        
	margLikelihood = w.computeMargLikelihood();
        if (margLikelihood == Double.POSITIVE_INFINITY) {
            margLikelihood = Double.NEGATIVE_INFINITY;
            System.out.println("Marginal likelihood evaluated as INFINITY");
        }
        
        //If resampling
        //xTrajSample.full = A.getSampleBack(x.matrix, resampleTimeIndexes);
        //Otherwise
        int randomParticle = (int) Math.round(Math.random() * (jParticles-1));
        xTrajSample.full = x.matrix.viewRow(randomParticle); //check dimensions of returned matrix
        
        //System.out.println("Coal events seen by likelihood: " + w.coalEventsCounted);
        //System.out.println("Move events seen by likelihood: " + w.moveEventsCounted);
    }//End Method
    * */
    
    public void runFilterForwardEntropy(int jNum, EpiModel epi, ZVectors dataZ, TreeNode[] tree, StructCoalModel coal, double startTime, double endTime, AbstractDistribution stndNorm, DoubleFactory2D factory2D, DoubleFactory3D factory3D)
    {
	jParticles = jNum;
        
        int zLocStart = dataZ.absoluteTimes.indexOf(startTime);
        int zLocEnd = dataZ.absoluteTimes.indexOf(endTime);
        xTrajSample = new XTrajectory(); //holds sampled trajectory obtained from particle filter
        xTrajSample.zLocStart = zLocStart;
        xTrajSample.zLocEnd = zLocEnd;
        DoubleArrayList xTimes = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zLocStart, zLocEnd);
	int totalTimes = xTimes.size();
	
        DoubleArrayList filterTimes = new DoubleArrayList();
        filterTimes.add(startTime);
        filterTimes.add(endTime);
        int filterTotalTimes = filterTimes.size();
        
	//Set up x matrix
        x = new StateTrajectory();
        x.getMatrix(jParticles, totalTimes, epi.getInitialConditions(), factory3D);
        
	//Set up w matrix
        Weights w = new Weights();
        w.getMatrix(jParticles, totalTimes, factory2D);
        
        //Need these to sample a trajectory from x
        int[] xLocs = new int[totalTimes];
        double timeStart; double timeEnd; 	
        int xLocStart = 0; int xLocEnd = 0;
	for (int time = 0; time < (filterTotalTimes - 1); ++time) {
            timeStart = filterTimes.get(time); //absolute start time
            timeEnd = filterTimes.get(time+1); //absolute end time
            zLocStart = dataZ.absoluteTimes.indexOf(timeStart);
            xLocStart = xTimes.indexOf(timeStart); //index of where to start in xTimes
            xLocs[time] = xLocStart;
            xLocEnd = xTimes.indexOf(timeEnd); //index of where to end in xTimes
            DoubleArrayList xSubTimes = (DoubleArrayList) xTimes.partFromTo(xLocStart, xLocEnd);
            DoubleArrayList xDtTimes = new DoubleArrayList();
            for (int i = 1; i < xSubTimes.size(); ++i) {
                double xDiff = xSubTimes.get(i) - xSubTimes.get(i-1);
		xDtTimes.add(xDiff);
            }

            //Update particle states
            x.propogateParticlesDeterministic(epi, time, xDtTimes, xSubTimes);
            //x.updateStates(theta, xLocStart, xDtTimes, xSubTimes, stndNorm);
	}
        xLocs[totalTimes-1] = xLocEnd;
        
        stateProbs = new LineageStateProbs();
        stateProbs.getMatrix(jParticles, tree.length, coal.Y.size(), factory3D);
        
        entropies = new LineageEntropies();
        entropies.getArrays(tree.length, coal.Y.size());
        
        int backCounter = 0;
        for (int time = (totalTimes-1); time > 0; --time) {
            
            int zLocEndInterval = zLocEnd - backCounter; //zLoc at end of interval (time = t+1)
            int zLocStartInterval = zLocEndInterval - 1;
            double timeEndInterval = dataZ.absoluteTimes.get(zLocEndInterval); 
            double timeStartInterval = dataZ.absoluteTimes.get(zLocStartInterval);
            
            //Get new samples (if samples)
            if (dataZ.omegaEvents.get(zLocEndInterval).contains(2)) { //if this interval begins with a sampling event
                stateProbs.addSamples(dataZ, zLocEndInterval, tree); //can get rid of lineageMap since just need sampling states
                entropies.addSamples(dataZ, zLocEndInterval, tree);
            }
            
            //Update weights
            double dtTime = dataZ.absoluteTimes.get(zLocEndInterval) - dataZ.absoluteTimes.get(zLocStartInterval);
            //w.updateWeightsBackNoCoalFromProbs(x.matrix.viewColumn(time), stateProbs, theta, time, coal, dtTime, timeEndInterval);
            w.updateWeightsBackNoCoalFromProbs(x.matrix.viewColumn(time), stateProbs, epi.params, time, coal, dtTime, timeEndInterval);
            
            //Update lineageStateProbs
            //stateProbs.updateProbs(x.matrix.viewColumn(time), theta, time, coal, dtTime, timeEndInterval);
            stateProbs.updateProbs(x.matrix.viewColumn(time), epi.params, time, coal, dtTime, timeEndInterval);
            entropies.updateEntropies(dataZ.absoluteTimes.get(zLocStartInterval), stateProbs);
            
            //Update weights after coalescent events
            if (dataZ.omegaEvents.get(zLocStartInterval).contains(1)) {
                //System.out.println();
                //w.updateWeightsBackCoalFromProbs(dataZ, zLocStartInterval, tree, x.matrix.viewColumn(time-1), stateProbs, theta, time, coal, timeStartInterval);
                w.updateWeightsBackCoalFromProbs(dataZ, zLocStartInterval, tree, x.matrix.viewColumn(time-1), stateProbs, epi.params, time, coal, timeStartInterval);
                entropies.updateEntropiesCoal(dataZ, zLocStartInterval, tree, stateProbs);
            }
            
            //Update back counter
            backCounter++;
        }
		
	margLikelihood = w.computeMargLikelihood();
        
        //If resampling
        //xTrajSample.full = A.getSample(x.matrix, zLocs);
        //Otherwise
        int randomParticle = (int) Math.round(Math.random() * (jParticles-1));
        //DoubleMatrix2D firstParticle = x.matrix.viewRow(0);
        xTrajSample.full = x.matrix.viewRow(randomParticle); //check dimensions of returned matrix
        
        //System.out.println("Coal events seen by likelihood: " + w.coalEventsCounted);
        //System.out.println("Move events seen by likelihood: " + w.moveEventsCounted);
    }//End Method

    
    private void resampleBackIntegrated(DoubleMatrix2D wCurr, int currIndex, int lastIndex) {
        
        //Take product of weights over resampling interval to compute current weights
        double iWeight;
        double weightSum = 0.0;
        double newStateProb;
        DoubleArrayList intervalWeights = new DoubleArrayList();
        Integer[] kIndexes = new Integer[jParticles];
        
        for (int j = 0; j < jParticles; ++j) {
            iWeight = 1.0;
            for (int t = lastIndex; t >= currIndex; t--) {
                iWeight *= wCurr.getQuick(j, t);    
            }
            intervalWeights.add(iWeight);
            weightSum += iWeight;
	}
        
        DoubleArrayList normCumWeights = new DoubleArrayList();
        double normWeight = 0.0;
        double cumNormWeight = 0.0;
        for (int j = 0; j < jParticles; ++j) {
            normWeight = intervalWeights.getQuick(j) / weightSum;
            cumNormWeight += normWeight;
            normCumWeights.add(cumNormWeight);
	}
        
        //Copy over old xCurr stateVariables (they will be overwrittern upon resampling)
        DoubleMatrix2D newXCurr = x.matrix.viewColumn(currIndex);
        
        //Copy over activeLineages in stateProbs (they will be overwritten upon resampling)
        DoubleMatrix3D stateProbsCopy = stateProbs.matrix.copy();
        
        for (int n = 0; n < jParticles; ++n) {
            double uniRand = Math.random();
            double cumVal = 0.0; 
            int counter = 0;
            while (cumVal <= uniRand)
            {
                cumVal = normCumWeights.getQuick(counter);
		counter++;
            }
            int newK = counter - 1;
            kIndexes[n] = newK; //-1 because counter has already been incremented
            
            //Update state variables
            for (int var = 0; var < newXCurr.rows(); ++var) {
                x.matrix.setQuick(var,n,currIndex, newXCurr.getQuick(var,newK)); //this indexing might be off depending on how viewColumn works on a 3D matrix
            }
            
            //Update lineage probs
            for (int i = 0; i < stateProbs.activeLineages.size(); i++) {
                int line = stateProbs.activeLineages.get(i);
                for (int stateK = 0; stateK < stateProbs.states; stateK++) {
                    newStateProb = stateProbsCopy.getQuick(newK, line, stateK);
                    stateProbs.matrix.setQuick(n,line,stateK,newStateProb);
                }
            }
	}
	A.indexes.add(kIndexes);    
    }
    
    private void resampleBackSimulated(DoubleMatrix2D xCurr, DoubleMatrix2D wCurr, int currIndex, int lastIndex) {
        
        //Take product of weights over resampling interval to compute current weights
        double iWeight;
        double weightSum = 0.0;
        DoubleArrayList intervalWeights = new DoubleArrayList();
        Integer[] kIndexes = new Integer[jParticles];
        
        for (int j = 0; j < jParticles; ++j) {
            iWeight = 1.0;
            for (int t = lastIndex; t >= currIndex; t--) {
                iWeight *= wCurr.getQuick(j, t);    
            }
            intervalWeights.add(iWeight);
            weightSum += iWeight;
	}
        
        DoubleArrayList normCumWeights = new DoubleArrayList();
        double normWeight = 0.0;
        double cumNormWeight = 0.0;
        for (int j = 0; j < jParticles; ++j) {
            normWeight = intervalWeights.getQuick(j) / weightSum;
            cumNormWeight += normWeight;
            normCumWeights.add(cumNormWeight);
	}
        
        //Copy over old xCurr stateVariables (they will be overwrittern upon resampling)
        DoubleMatrix2D newXCurr = xCurr.copy();
        
        //Copy over activeLineages in stateArray (they will be overwritten upon resampling)
        ArrayList<ArrayList<Integer>> stateArrayCopy = new ArrayList<ArrayList<Integer>>();
        int lines = stateArray.array.size(); 
        for (int i = 0; i < lines; i++) {
            //stateArrayCopy.add(new ArrayList<Integer>());
            //stateArrayCopy.get(i).addAll(stateArray.array.get(i));
            stateArrayCopy.add(stateArray.array.get(i));
        }
        
        for (int n = 0; n < jParticles; ++n) {
            double uniRand = Math.random();
            double cumVal = 0.0; 
            int counter = 0;
            while (cumVal <= uniRand)
            {
                cumVal = normCumWeights.getQuick(counter);
		counter++;
            }
            int newK = counter - 1; //-1 because counter has already been incremented
            kIndexes[n] = newK; 
            
            //Update state variables
            for (int var = 0; var < xCurr.rows(); ++var) {
                newXCurr.setQuick(var,n, xCurr.getQuick(var,newK)); //this indexing might be off depending on how viewColumn works on a 3D matrix
            }
            
            //Update lineage states
            for (int i = 0; i < stateArray.activeLineages.size(); i++) {
                int line = stateArray.activeLineages.get(i);
                stateArray.array.get(line).set(n,stateArrayCopy.get(line).get(newK));
            }
            xCurr = newXCurr;
	}
	//A.indexes.add(kIndexes);    
    }
    
    public void resample(DoubleMatrix2D xCurr, DoubleMatrix1D wCurr)
    {
        //Currently only implements multinomial resampling with replacement
        
	double wSum = wCurr.zSum();
	double currSum = 0;	
	DoubleMatrix1D wNorm = wCurr.copy(); //factory1D.make(jParticles);
	DoubleMatrix1D wCum =  wCurr.copy(); //factory1D.make(jParticles);
        Integer[] kIndexes = new Integer[jParticles];
		
	for (int j = 0; j < jParticles; ++j) {
            wNorm.set(j,(wCurr.getQuick(j)/wSum)); //normalized particle weights
            currSum += wNorm.getQuick(j);
            wCum.setQuick(j,currSum);
	}	
	
        DoubleMatrix2D newXCurr = xCurr.copy();
	for (int n = 0; n < jParticles; ++n)
	{
            double uniRand = Math.random();
            double cumVal = 0.0; 
            int counter = 0;
            while (cumVal <= uniRand)
            {
                cumVal = wCum.getQuick(counter);
		counter++;
            }
            kIndexes[n] = counter - 1; //-1 because counter has already been incremented
            for (int var = 0; var < xCurr.rows(); ++var) {
                newXCurr.setQuick(var,n, xCurr.getQuick(var,(counter-1))); //this indexing might be off depending on how viewColumn works on a 3D matrix
            }
            xCurr = newXCurr;
	}
	A.indexes.add(kIndexes);
    }//END method
}