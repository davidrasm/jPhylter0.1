import java.util.ArrayList;
import cern.colt.matrix.*;
import cern.colt.list.*;
import java.util.Collections;
import java.util.HashSet;

/**
 * Track ancestry of particles for use in sampling a particle trajectory when resampling
 * @author David
 */
public class ParticleAncestry {
    
    ArrayList<Integer[]> indexes = new ArrayList<Integer[]>(); 
    int jParticles;
    int times;
    ArrayList<ArrayList<Integer>> uniqueParticles;
    
    public void getArray(int jNum, int ts)
    {
        jParticles = jNum;
        times = ts-1;
    }
    
    
    //THIS METHOD IS NOT USED IN THE PARTICLE GIBBS SAMPLER
    public DoubleMatrix2D getSample(DoubleMatrix3D x, ArrayList<Integer> resampleTimeIndexes)
    {
        
        int finalTimeIndex = x.viewSlice(0).viewRow(0).size();
        if (!resampleTimeIndexes.contains(finalTimeIndex-1)) {
            resampleTimeIndexes.add(finalTimeIndex-1);
        }
        times = resampleTimeIndexes.size();
        int slices = x.slices();
        
        //Randomly sample particle by final weights
        int randomParticle = (int) Math.round(Math.random() * (jParticles - 1));
        
        //Trace particle ancestry backwards in time
        int[] path = new int[times-1];
        path[times - 2] = randomParticle;
        ArrayList<Integer> pathArray = new ArrayList<Integer>();
        pathArray.add(randomParticle);
        for (int t = times - 2; t > 0; --t) {
            int parentIndex = path[t];
            path[t-1] = indexes.get(t-1)[parentIndex];
            pathArray.add(indexes.get(t-1)[parentIndex]);
        }
        //System.out.println("Start particle: " + path[0]);
        //System.out.println("Particle path: " + pathArray);
        
        //Get particle state trajectory
        DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
        DoubleMatrix2D xSample = factory2D.make(slices, x.columns());
        for (int n = 0; n < (times-1); n++) {
            int xLocStart = resampleTimeIndexes.get(n);
            int xLocEnd = resampleTimeIndexes.get(n+1);
            for (int t = xLocStart; t < xLocEnd; t++) {
                for (int i = 0; i < slices; ++i) {
                    xSample.set(i, t, x.getQuick(i, path[n], t));
                }   
            }
        }
        for (int i = 0; i < slices; ++i) {
            xSample.set(i, (finalTimeIndex-1), x.getQuick(i, path[times-2], (finalTimeIndex-1)));
        } 
        
        return xSample;
        
    }
    
    //THIS METHOD IS NOT USED IN THE PARTICLE GIBBS SAMPLER
    public DoubleMatrix2D getSampleBack(DoubleMatrix3D x, ArrayList<Integer> resampleTimeIndexes)
    {
        
        if (resampleTimeIndexes.contains(0)) {
            //Then don't add it
        } else {
            resampleTimeIndexes.add(0);
        }
        Collections.reverse(resampleTimeIndexes);
        Collections.reverse(indexes);
        times = resampleTimeIndexes.size();
        int slices = x.slices();
        
        
        int randomParticle = (int) Math.round(Math.random() * (jParticles - 1));
        int[] path = new int[times - 1];
        path[0] = randomParticle;
        for (int t = 0; t < (times - 2); t++) {
            int childIndex = path[t];
            path[t+1] = indexes.get(t)[childIndex];    
        }
        
        DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
        DoubleMatrix2D xSample = factory2D.make(slices, x.columns());
        for (int n = 0; n < (times - 1); ++n) {
            int xLocStart = resampleTimeIndexes.get(n);
            int xLocEnd = resampleTimeIndexes.get(n+1);
            for (int t = xLocStart; t <= xLocEnd; ++t) {
                for (int i = 0; i < slices; ++i) {
                    xSample.set(i, t, x.getQuick(i, path[n], t));
                }   
            }  
        }
        //System.out.println();
        return xSample;
    }
    
    
    public ArrayList<ArrayList<Integer>> getParticleAncestryMatrix(DoubleMatrix3D x, ArrayList<Integer> resampleTimeIndexes)
    {
        int finalTimeIndex = x.viewSlice(0).viewRow(0).size();
        //if (resampleTimeIndexes.contains(finalTimeIndex-1)) {
            
        //} else {
            //resampleTimeIndexes.add(finalTimeIndex-1);
        //}
        times = resampleTimeIndexes.size();
        
        int xLocStart; int xLocEnd;
        ArrayList<ArrayList<Integer>> ancestryMatrix = new ArrayList<ArrayList<Integer>>();
        for (int t = 0; t < finalTimeIndex; t++) {
            ancestryMatrix.add(new ArrayList<Integer>());
        }
        int parentIndex; int endParticle;
        for (int p = 0; p < jParticles; p++) {
            endParticle = indexes.get(times-2)[p];
            int[] path = new int[times-1];
            path[times-2] = endParticle;
            for (int t = times - 2; t > 0; --t) {
                parentIndex = path[t];
                path[t-1] = indexes.get(t-1)[parentIndex];
            }
            for (int n = 0; n < (times-1); n++) {
                xLocStart = resampleTimeIndexes.get(n);
                xLocEnd = resampleTimeIndexes.get(n+1);
                if (n == times-2) {
                    xLocEnd = resampleTimeIndexes.get(n+1) + 1;
                }
                for (int t = xLocStart; t < xLocEnd; t++) {
                    ancestryMatrix.get(t).add(path[n]); 
                }
            }
            
        }
        
        //Determine unique particles in ancestryMatrix---these won't necessarily be in order
        uniqueParticles = new ArrayList<ArrayList<Integer>>();
        for (int t = 0; t < finalTimeIndex; t++) {
            HashSet<Integer> uniqueSet = new HashSet<Integer>(ancestryMatrix.get(t));
            uniqueParticles.add(new ArrayList<Integer>());
            uniqueParticles.get(t).addAll(uniqueSet);  
        }
        
        return ancestryMatrix;
        
    }
    
    public ArrayList<ArrayList<Integer>> getNullAncestryMatrix(DoubleMatrix3D x) 
    {
        ArrayList<ArrayList<Integer>> ancestryMatrix = new ArrayList<ArrayList<Integer>>();
        
        for (int t = 0; t < (times+1); t++) {
            ancestryMatrix.add(new ArrayList<Integer>());
            for (int p = 0; p < jParticles; p++) {
                ancestryMatrix.get(t).add(p);
            }
        }
        
        uniqueParticles = new ArrayList<ArrayList<Integer>>();
        for (int t = 0; t < (times+1); t++) {
            HashSet<Integer> uniqueSet = new HashSet<Integer>(ancestryMatrix.get(t));
            uniqueParticles.add(new ArrayList<Integer>());
            uniqueParticles.get(t).addAll(uniqueSet);  
        }
        
        return ancestryMatrix;
        
    }
    
}
