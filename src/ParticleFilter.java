import cern.colt.list.*;
import cern.colt.matrix.*;
import cern.jet.random.*;
import cern.jet.stat.Descriptive;
import java.util.ArrayList;

/**
 * Run particle filter (with or without resampling)
 * @author David
 */

public class ParticleFilter
{
    int jParticles; //number of particles
    double margLikelihood;
    ParticleAncestry A;
    StateTrajectory x;
    XTrajectory xTrajSample;
    DoubleFactory2D factory2D = DoubleFactory2D.dense;
    DoubleFactory3D factory3D = DoubleFactory3D.dense;
    
    public void runFilter(int jNum, EpiModel epi, ZVectors dataZ, CoalModel coal, double startTime, double endTime, TimeSeries timeSeries)
    {
	jParticles = jNum;
        
        //If startTime endTime equal negative/positive infinity have to update start and end times
        int zLocStart = dataZ.absoluteTimes.indexOf(startTime);
        int zLocEnd = dataZ.absoluteTimes.indexOf(endTime);
        xTrajSample = new XTrajectory(); //holds sampled trajectory obtained from particle filter
        xTrajSample.zLocStart = zLocStart;
        xTrajSample.zLocEnd = zLocEnd;
        DoubleArrayList xTimes = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zLocStart, zLocEnd);
	int totalTimes = xTimes.size();
        
        //If starting at specific time:
        double simStartTime = epi.params.get(2);
        double currTime = startTime;
        int counter = 0;
        while (currTime <= simStartTime) {
            counter++;
            currTime = dataZ.absoluteTimes.get(counter);
        }
        int simStartIndex = counter-1;
        epi.startTimeIndex = simStartIndex;
        //epi.startTimeIndex = 0;
        
        //Set resampling times
        //double firstTime = 734015.0;
        ArrayList<Integer> resamplingTimes = new ArrayList<Integer>();
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 0.0 + 0.25));
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 0.0 + 0.5));
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 0.0 + 1.0));
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 0.0 + 2.0));
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 0.0 + 4.0));
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 0.0 + 6.0));
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 0.0 + 8.0));
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, 0.0 + 10.0));
        
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, firstTime + 28.0));
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, firstTime + 48.0));
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, firstTime + 76.0));
        //resamplingTimes.add(ArrayUtils.getIndexOfClosestDouble(dataZ.absoluteTimes, firstTime + 104.0));
        
	//Set up x matrix
        x = new StateTrajectory();
        //x.getMatrix(jParticles, totalTimes, epi.getInitialConditions(), factory3D);
        //OR
        x.getMatrixWithOffset(jParticles, totalTimes, epi.getInitialConditions(), factory3D, simStartIndex);
        
	//Set up w matrix
        Weights w = new Weights();
        w.getMatrix(jParticles, totalTimes, factory2D);
        w.setCumulativeWeights();
        
        //Set up A matrix for particle ancestry
        A = new ParticleAncestry();
        A.getArray(jParticles, totalTimes);
        
        //Run filter
        DoubleArrayList xSubTimes = new DoubleArrayList(); DoubleArrayList xDtTimes = new DoubleArrayList();
        double timeStart; double timeEnd; 	
        ArrayList<Integer> resampleTimeIndexes = new ArrayList<Integer>();//Need these to sample a trajectory from x 
        resampleTimeIndexes.add(0);
	for (int time = 0; time < (totalTimes - 1); ++time) {
            
            timeStart = xTimes.get(time); //absolute start time
            timeEnd = xTimes.get(time+1); //absolute end time
            zLocStart = dataZ.absoluteTimes.indexOf(timeStart);
            
            xSubTimes.add(timeStart); xSubTimes.add(timeEnd);
            xDtTimes.add(timeEnd - timeStart);

            //Update particle states
            x.matrix = epi.updateStates(jParticles, x.matrix, time, xDtTimes, xSubTimes); //simulate from process model
            //x.updateStates(theta, xLocStart, xDtTimes, xSubTimes, stndNorm);
            
            w.updateWeightsSingleStep(x.matrix, epi.params, time, dataZ, xDtTimes, xSubTimes, zLocStart, time, coal); //pass ZVectors for 2nd tree 
            if (dataZ.omegaEvents.get(zLocStart+1).contains(3)) {
                //interval ends in observation event in timeSeries
                w.updateWeightsTimeSeries(x.matrix, time, zLocStart, dataZ, epi, timeSeries, xTrajSample.zLocStart);
            }
            
            
            w.updateCumulativeWeights(time); //cumulative weights are log weights to avoid underflow
            
            //double avgCumulativeWeight = Descriptive.mean(w.cumulativeWeights);
            //System.out.println(time + ": " + avgCumulativeWeight);
            //double ESS = w.computeEffectiveSampleSize();
            //System.out.println("ESS = " + ESS + " at time: " + timeEnd);
            
            //Resample here if resampling
            if (resamplingTimes.contains(time)) {
                
                this.resample(w.cumulativeWeights, time);
                resampleTimeIndexes.add(time+1);
                w.resetCumulativeWeights();
            }
            xSubTimes.removeAll(xSubTimes); xDtTimes.removeAll(xDtTimes);
            
	}
        //xLocs[filterTotalTimes-1] = xLocEnd; necessary
		
	//Compute marg likelihood estimate (if not resampling)
        //margLikelihood = w.computeMargLikelihood();
        
        //Compute marg likelihood estimate when resampling
        if (! resampleTimeIndexes.contains(totalTimes-1)) { //Resample again if didn't resample at last time point
            this.resample(w.cumulativeWeights, (totalTimes-2)); //-2 b/c -1 for indexing and -1 because timeIndex is timeIndex+1 in resampling method
            resampleTimeIndexes.add(totalTimes-1);
        }
        ArrayList<ArrayList<Integer>> ancestryMatrix = A.getParticleAncestryMatrix(x.matrix, resampleTimeIndexes);
        margLikelihood = w.computeMargLikelihoodWithAncestry(ancestryMatrix);
        if (margLikelihood == Double.POSITIVE_INFINITY) {
            margLikelihood = Double.NEGATIVE_INFINITY;
        }
        
        //Sample particle state trajectory when resampling
        xTrajSample.full = A.getSample(x.matrix, resampleTimeIndexes);
        
        //Otherwise
        //int randomParticle = (int) Math.round(Math.random() * (jParticles-1));
        //xTrajSample.full = x.matrix.viewRow(randomParticle); //check dimensions of returned matrix
        //System.out.println();
        
    }//End Method
	
    public void resample(DoubleArrayList cumulativeWeights, int time)
    {
        
        //Get weights on linear scale before resampling
        DoubleArrayList expWeights = new DoubleArrayList();
        double maxLogWeights = Descriptive.max(cumulativeWeights);
        for (int p = 0; p < jParticles; p++) {
            expWeights.add(Math.exp(cumulativeWeights.get(p) - maxLogWeights)); //recenter before exponentiating to prevent overflow/underflow
        }
        
        DoubleArrayList normCumWeights = new DoubleArrayList();
        double sumExpWeights = Descriptive.sum(expWeights);
        double normWeight;
        double cumNormWeight = 0.0;
        for (int j = 0; j < jParticles; ++j) {
            normWeight = expWeights.getQuick(j) / sumExpWeights;
            cumNormWeight += normWeight;
            normCumWeights.add(cumNormWeight);
	}
        
        Integer[] kIndexes = new Integer[jParticles];

        
        //Copy over old xCurr stateVariables (they will be overwrittern upon resampling)
        DoubleMatrix2D xCurrView = x.matrix.viewColumn(time+1);
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
                x.matrix.setQuick(var,n,(time+1), newXCurr.getQuick(var,newK)); //this indexing might be off depending on how viewColumn works on a 3D matrix
            }
	} 
	A.indexes.add(kIndexes);
	
    }//END method
			
} //END class