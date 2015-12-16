import cern.jet.random.*;
import cern.colt.list.*;
import cern.colt.matrix.*;

/**
 * Particle state trajectories
 * Modify updateStates according to process model
 * @author David
 */
public class StateTrajectory {
    
    DoubleMatrix3D matrix;
    int jParticles;
    int times;
    int slices;
    
    public void getMatrix(int jNum, int ts, DoubleArrayList xInit, DoubleFactory3D factory3D)
    {
        slices = xInit.size();
        jParticles = jNum;
        times = ts;
	matrix = factory3D.make(slices, jParticles, times);
	matrix.assign(Double.NaN);
	for (int p = 0; p < slices; ++p) {
            for (int m = 0; m < jParticles; ++m) {
                matrix.setQuick(p, m, 0, xInit.get(p));
            }
	}
    }
    
    public void getMatrixWithOffset(int jNum, int ts, DoubleArrayList xInit, DoubleFactory3D factory3D, int offsetIndex)
    {
        slices = xInit.size();
        jParticles = jNum;
        times = ts;
	matrix = factory3D.make(slices, jParticles, times);
	matrix.assign(0.0);
	for (int p = 0; p < slices; ++p) {
            for (int m = 0; m < jParticles; ++m) {
                matrix.setQuick(p, m, offsetIndex, xInit.get(p));
            }
	}
    }
    
    public void getMatrixWithRandomInits(int jNum, int ts, EpiModel epi, DoubleFactory3D factory3D)
    {

//        DoubleArrayList determInits = epi.getInitialConditions();
//        slices = determInits.size();
//        jParticles = jNum;
//        times = ts;
//	matrix = factory3D.make(slices, jParticles, times);
//	matrix.assign(Double.NaN);
//        for (int m = 0; m < jParticles; ++m) {
//            DoubleArrayList particleInits = epi.getRandomInitialConditions();
//            for (int p = 0; p < slices; ++p) {
//                matrix.setQuick(p, m, 0, particleInits.get(p));
//            }
//        }
        
    }
    
    /**
    public void getConditionalMatrix(XTrajectory prevXTraj) {
        
        for(int t = 0; t < times; t++) {
            for(int p = 0; p < slices; p++) {
                matrix.setQuick(p,(jParticles-1),t, prevXTraj.full.getQuick(p,t));
            }
        }
        
    } */
    
    /**
    public void getMatrixBack(int jNum, int ts, DoubleArrayList xInit, DoubleFactory3D factory3D)
    {
        
        //Same as getMatrix but sets initial conditions at the final endTime
        
        slices = xInit.size();
        jParticles = jNum;
        times = ts;
	matrix = factory3D.make(slices, jParticles, times);
	matrix.assign(Double.NaN);
	for (int p = 0; p < slices; ++p) {
            for (int m = 0; m < jParticles; ++m) {
                matrix.setQuick(p, m, (times - 1), xInit.get(p));
            }
	}
    } */
    
    public void propogateParticles(EpiModel epi, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes)
    {   
        //Pass xMatrix to EpiModel to update states
        matrix = epi.updateStatesTauLeap(jParticles, matrix, xLocStart, xDtTimes, xSubTimes);
    }
    
    public void propogateParticlesDeterministic(EpiModel epi, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes)
    {   
        //Pass xMatrix to EpiModel to update states
        matrix = epi.updateStatesDeterministic(jParticles, matrix, xLocStart, xDtTimes, xSubTimes);
    }
    
    /**
    public void propogateParticlesConditional(EpiModel epi, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes)
    {   
        //Pass xMatrix to EpiModel to update states
        //matrix = epi.updateStates((jParticles-1), matrix, xLocStart, xDtTimes, xSubTimes, stndNorm);
        matrix = epi.updateStatesTauLeap((jParticles-1), matrix, xLocStart, xDtTimes, xSubTimes);
    } */
    
    
    /**
    public void propogateSingleParticle(EpiModel epi, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes, int particle)
    {   
        //Pass xMatrix to EpiModel to update states
        //matrix = epi.updateStatesTauLeapSingleParticle(particle, matrix, xLocStart, xDtTimes, xSubTimes);
    } */
    
    
    //These last methods were just for debugging
    
    public void viewParticleStates(int n)
    {
        DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
	DoubleMatrix2D nMatrix = factory2D.make(jParticles, slices);
	nMatrix = matrix.viewColumn(n);
	System.out.println("Particle states at time =" + n);
	System.out.println(nMatrix);
    } //END method
	
	
    public void viewParticleTrajectory(int n)
    {
	DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
	DoubleMatrix2D nMatrix = factory2D.make(slices, times);
	nMatrix = matrix.viewRow(n);
	System.out.println("Trajectory for particle number" + n);
	System.out.println(nMatrix);
    } //END method
    
}
