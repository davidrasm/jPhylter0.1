import cern.jet.random.*;
import cern.colt.list.*;
import cern.colt.matrix.*;
import java.io.*; //IOException;

/**
 * Run MCMC sampling trees, epi params and evo params
 * @author David
 */

public class LikelihoodSurface
{
        public void getSurface(int particles, TreeNode[] treeStart, EpiModel epi, double startTime, double endTime, double dt, DoubleArrayList thetaEstStart, DoubleArrayList thetaFixed) throws IOException
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
                StructCoalModel coal = new StructCoalModelSECA(); //modify coal model to fit process model
                coal.make();
                
                //Get initial  ZVector
                ZVectors dataZNow = new ZVectors(); //Vector data structure for tree and observation data
                dataZNow.setZVectors(treeStart, startTime, endTime, dt); //do this first without the moveTimes
                startTime = dataZNow.startTime;
                endTime = dataZNow.endTime;
                
                //If settin lineage states
                //LineageStateMap lineageMapNow = new LineageStateMap();
                //lineageMapNow.getMoveTimesFromTree(treeStart, dataZNow); //not sure how this handles the root node
                //dataZNow.setZVectors(treeStart, lineageMapNow, startTime, endTime, dt); //do this again with moveTimes
                //lineageMapNow.getStatesFromTree(treeStart, dataZNow);

                StructParticleFilter filter = new StructParticleFilter();
                
                //Set up for param1 (betaW)
                DoubleArrayList param1List = new DoubleArrayList();
                double param1LowerBound = 0.08;
                double param1UpperBound = 0.15;
                double param1Increment = 0.001;
                double currValueParam1 = param1LowerBound;
                param1List.add(currValueParam1);
                while (currValueParam1 < param1UpperBound) {
                    currValueParam1 += param1Increment;
                    param1List.add(currValueParam1);
                }
                
                //Set up for param2
                DoubleArrayList param2List = new DoubleArrayList();
                double param2LowerBound = 0.0055; //0.0;
                double param2UpperBound = 0.0055;
                double param2Increment = 0.001;
                double currValueParam2 = param2LowerBound;
                param2List.add(currValueParam2);
                while (currValueParam2 < param2UpperBound) {
                    currValueParam2 += param2Increment;
                    param2List.add(currValueParam2);
                }
                //param2List.set(0, 0.001);
                
                //Set up params and likelihood matrices
                int param1Vals = param1List.size();
                int param2Vals = param2List.size();
                DoubleMatrix2D param1Matrix = factory2D.make(param2Vals, param1Vals);
                DoubleMatrix2D param2Matrix = factory2D.make(param2Vals, param1Vals);
                DoubleMatrix2D likelihoodMatrix = factory2D.make(param2Vals, param1Vals);
                
                int counter = 0;
                int totalCounts = param1Vals * param2Vals;
                for (int param1Index = 0; param1Index < param1Vals; param1Index++) {
                    for (int param2Index = 0; param2Index < param2Vals; param2Index++) {
                        
                        counter++;
                        System.out.println("Current likelihood evaluation: " + counter + " of " + totalCounts);
                        
                        double param1Now = param1List.get(param1Index);
                        double param2Now = param2List.get(param2Index);
                        DoubleArrayList newParamList = new DoubleArrayList();
                        newParamList.add(param1Now); newParamList.add(param2Now);
                        epi.updateEstParams(newParamList);
                        
                        //filter.runFilterBackIntegrated(particles, epi, dataZNow, treeStart, coal, startTime, endTime, stndNorm, factory2D, factory3D);
                        filter.runFilterForwardBack(particles, epi, dataZNow, treeStart, coal, startTime, endTime, stndNorm, factory2D, factory3D);
                        
                        //XTrajectory xTrajNow = filter.xTrajSample;
                        //String xTrajLabel = "IfTrajStart.jpg";
                        //MCMCPlotter.getParticleTrace(1, xTrajNow.full, xTrajLabel);
                        //xTrajLabel = "IgTrajStart.jpg";
                        //MCMCPlotter.getParticleTrace(3, xTrajNow.full, xTrajLabel);
                        //LineageEntropies entropiesStart = filter.entropies;
                        
                        double pGThetaNow = filter.margLikelihood;
                        System.out.println("Current log likelihood: " + pGThetaNow);
                        param1Matrix.set(param2Index, param1Index, param1Now);
                        param2Matrix.set(param2Index, param1Index, param2Now);
                        likelihoodMatrix.set(param2Index, param1Index, pGThetaNow);
                        
                    }
                }
                //Write out thetaSamples and psiSamples to CSV
                WriteMatrixToCSV csvWriter = new WriteMatrixToCSV();
                String outputFile1 = "LikeSurfParam1MatrixHIVSimTree14"; 
                csvWriter.writeCSV(param1Matrix, outputFile1);

                String outputFile2 = "LikeSurfParam2MatrixHIVSimTree14"; 
                csvWriter.writeCSV(param2Matrix, outputFile2);
                
                String outputFile3 = "LikeSurfLikeMatrixHIVSimTree14"; 
                csvWriter.writeCSV(likelihoodMatrix, outputFile3);
                
                //Write out entropy tree
//                NewickWriter treeWriter = new NewickWriter();
//                FileWriter fw = new FileWriter("EntropyTree2PatchNoNoise02.tre");
//                PrintWriter pw = new PrintWriter(fw);
//                filter.runFilterForwardEntropy(particles, epi, dataZNow, treeStart, coal, startTime, endTime, stndNorm, factory2D, factory3D);
//                LineageEntropies entropiesStart = filter.entropies;
//                treeWriter.oldEntropyTreeToString(treeStart, entropiesStart);
//                pw.print(treeWriter.treeString);
//                pw.println();
//                pw.close();
//                fw.close();
                
	}//END method

}
