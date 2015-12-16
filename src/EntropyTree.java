import cern.jet.random.*;
import cern.colt.list.*;
import cern.colt.matrix.*;
import java.io.*; //IOException;

/**
 * Get Entropy tree for params in EpiModel
 * @author David
 */

public class EntropyTree
{
    public void getTree(int particles, TreeNode[] treeStart, EpiModel epi, double startTime, double endTime, double dt, DoubleArrayList thetaEstStart, DoubleArrayList thetaFixed) throws IOException
    {
        //Set random number generator and matrix factories
        cern.jet.random.engine.RandomEngine engine = new cern.jet.random.engine.MersenneTwister(); 
        AbstractDistribution stndNorm = new Normal(0.0, 1.0, engine);

        //Initialize colt matrix factories
        DoubleFactory2D factory2D;
        factory2D = DoubleFactory2D.dense;
        DoubleFactory3D factory3D;
        factory3D = DoubleFactory3D.dense;

        //Get coal model
        StructCoalModel coal = new StructCoalModel2Pop(); //modify coal model to fit process model
        coal.make();

        //Get initial  ZVector
        ZVectors dataZNow = new ZVectors(); //Vector data structure for tree and observation data
        dataZNow.setZVectors(treeStart, startTime, endTime, dt); //do this first without the moveTimes
        startTime = dataZNow.startTime;
        endTime = dataZNow.endTime;

        StructParticleFilter filter = new StructParticleFilter();
        particles = 1; //only need one trajectory
        //filter.runFilterBackIntegrated(particles, epi, dataZNow, treeStart, coal, startTime, endTime, stndNorm, factory2D, factory3D);
        //filter.runFilterForwardBack(particles, epi, dataZNow, treeStart, coal, startTime, endTime, stndNorm, factory2D, factory3D);
        //filter.runFilterForwardEntropy(particles, epi, dataZNow, treeStart, coal, startTime, endTime, stndNorm, factory2D, factory3D);
        //LineageEntropies entropies = filter.entropies;

        //Write out entropy tree
        NewickWriter treeWriter = new NewickWriter();
        FileWriter fw = new FileWriter("PISTestTree_2PopStochRhoH_entropyTree.tre");
        PrintWriter pw = new PrintWriter(fw);
        filter.runFilterForwardEntropy(particles, epi, dataZNow, treeStart, coal, startTime, endTime, stndNorm, factory2D, factory3D);
        LineageEntropies entropiesStart = filter.entropies;
        treeWriter.oldEntropyTreeToString(treeStart, entropiesStart);
        pw.print(treeWriter.treeString);
        pw.println();
        pw.close();
        fw.close();

    }//END method

}
