import cern.jet.random.*;
import cern.colt.list.*;
import cern.colt.matrix.*;
import java.util.ArrayList;
import java.io.*; //IOException;

/**
 * Run MCMC using deterministic models
 * @author David
 */

public class ClassicMCMC
{

    //public void runMCMC(int iterations, int particles, TreeNode[] treeStart, SeqAlign alignment, GTREvoModel model, double startTime, double endTime, double dt, DoubleArrayList thetaEstStart, DoubleArrayList thetaFixed, ProposalDist thetaProposal, DoubleArrayList psiEstNow, DoubleArrayList psiFixed, ProposalDist psiProposal, DoubleArrayList piNow, ProposalStep pStep, EpiModel epi) throws IOException
    public void runMCMC(int iterations, int particles, TreeNode[] treeStart, SeqAlign alignment, double startTime, double endTime, double dt, DoubleArrayList thetaEstStart, DoubleArrayList thetaFixed, ProposalDist thetaProposal, DoubleArrayList psiEstNow, DoubleArrayList psiFixed, ProposalDist psiProposal, DoubleArrayList piNow, EpiModel epi) throws IOException
    {
        
        boolean verbose = true;

        //Tree sample output file
        String fileNameStem = "DENV1_subMixed1130"; // + Integer.toString(run+1);
        String fileNamePostFix = "structVectored_082513";
        NewickWriter treeWriter = new NewickWriter();
        String treeFileName = fileNameStem + "_trees_" + fileNamePostFix + ".tre";
        FileWriter treeFileWriter = new FileWriter(treeFileName);
        PrintWriter treePrintWriter = new PrintWriter(treeFileWriter);

        //Theta sample output file
        SampleWriter sampleWriter = new SampleWriter();
        String paramFileName = fileNameStem + "_params_" + fileNamePostFix;
        FileWriter thetaFileWriter = new FileWriter(paramFileName);
        PrintWriter thetaPrintWriter = new PrintWriter(thetaFileWriter);

        String xTrajFileName = fileNameStem + "_xTrajs_" + fileNamePostFix;
        FileWriter xTrajFileWriter = new FileWriter(xTrajFileName);
        PrintWriter xTrajPrintWriter = new PrintWriter(xTrajFileWriter);

        //Likelihood output file
        String likeFileName = fileNameStem + "_likes_" + fileNamePostFix;
        FileWriter likeFileWriter = new FileWriter(likeFileName);
        PrintWriter likePrintWriter = new PrintWriter(likeFileWriter);

        //Initialize matrix for parameter samples
        int paramThin = 500; //when to save param samples

        DoubleArrayList thetaEstNow = thetaEstStart.copy(); //DONT NEED THIS
        DoubleArrayList thetaAllNow = thetaEstNow.copy(); 
        thetaAllNow.addAllOf(thetaFixed);
        sampleWriter.samplesToString(thetaEstNow);
        thetaPrintWriter.print(sampleWriter.sampleString); thetaPrintWriter.flush(); thetaPrintWriter.println();

        //Initialize MCMC timer and counters
        double startClockTime = System.currentTimeMillis();
        int thetaProposalCount = 0;
        int thetaAcceptCount = 0;

        //Get coal model
        StructCoalModel coal = new StructCoalModelDengueVectored(); //modify coal model to fit process model
        coal.make();

        //Set initial ZVectors
        ZVectors dataZNow = new ZVectors(); //Vector data structure for tree and observation data
        dataZNow.setZVectors(treeStart, startTime, endTime, dt);
        startTime = dataZNow.startTime;
        endTime = dataZNow.endTime;
        
        //Set random number generator and matrix factories
        cern.jet.random.engine.RandomEngine engine = new cern.jet.random.engine.MersenneTwister(); 
        AbstractDistribution stndNorm = new Normal(0.0, 1.0, engine);

        //Initialize colt matrix factories
        DoubleFactory2D factory2D;
        factory2D = DoubleFactory2D.dense;
        DoubleFactory3D factory3D;
        factory3D = DoubleFactory3D.dense;
        
        StructParticleFilter filter = new StructParticleFilter();
        filter.runFilterForwardBack(particles, epi, dataZNow, treeStart, coal, startTime, endTime, stndNorm, factory2D, factory3D);
        XTrajectory xTrajNow = filter.xTrajSample; //this should not be updated after this point!!
        double pGNow = filter.margLikelihood;

        //Initialize xTrajSamples for xTraj samples
        XTrajSamples xTrajs = new XTrajSamples();
        double sampleFreq = 30.0; //do this in weeks to align with sim times
        int sampleVar1 = 22; //HCMC
        int sampleVar2 = 23; //non-HCMC
        
        //For prevalence
//       xTrajs.getSampleTimes(startTime, endTime, sampleFreq);
//        DoubleArrayList xTrajSampleNow = xTrajs.getSample(sampleVar1, dataZNow, xTrajNow);
//        sampleWriter.samplesToString(xTrajs.sampleTimes);
//        xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
//        sampleWriter.samplesToString(xTrajSampleNow);
//        xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
//        xTrajSampleNow = xTrajs.getSample(sampleVar2, dataZNow, xTrajNow);
//        sampleWriter.samplesToString(xTrajSampleNow);
//        xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();

        //Print times
        xTrajs.getSampleTimes(startTime, endTime, sampleFreq);
        sampleWriter.samplesToString(xTrajs.sampleTimes);
        xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();

        DoubleArrayList xTrajSampleNow = xTrajs.getSampleIncidence(sampleVar1, dataZNow, xTrajNow);
        sampleWriter.samplesToString(xTrajSampleNow);
        xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();

        xTrajSampleNow = xTrajs.getSampleIncidence(sampleVar2, dataZNow, xTrajNow);
        sampleWriter.samplesToString(xTrajSampleNow);
        xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
        
        
        //Can view state trajectories for initial particle filter run
        //String xTrajLabel = "xTrajStart.jpg";
        //MCMCPlotter.getParticleTrace(1, xTrajNow.full, xTrajLabel);

        TreeNode[] treeNow = treeStart;

        //Write out starting likelihoods
        DoubleArrayList likeSamples = new DoubleArrayList();
        //double totalLikeNow = pGThetaNow + pDGPsiNow;
        likeSamples.add(pGNow); //likeSamples.add(pXNow); likeSamples.add(pLgNow); //likeSamples.add(totalLikeNow);
        sampleWriter.samplesToString(likeSamples);
        likePrintWriter.print(0 + ", "); likePrintWriter.print(sampleWriter.sampleString); likePrintWriter.flush(); likePrintWriter.println();

        //Set up autoOptimization routines
        AutoOptimizer optimizer = new AutoOptimizer();
        optimizer.autoOptimize = true;
        if (optimizer.autoOptimize) {
            ArrayList<Integer> runPoints = new ArrayList<Integer>();
            runPoints.add(200); runPoints.add(500); runPoints.add(1000); runPoints.add(2000); runPoints.add(10000); runPoints.add(50000); runPoints.add(100000);
            optimizer.setRunPoints(runPoints, thetaEstNow.size());
            optimizer.getSampleMatrix();
        }

        //Run MCMC
        int paramCounter = 0;
        for (int n = 1; n <= iterations; ++n) {

                paramCounter++;
                if (verbose) {
                    System.out.println("MCMC iteration = " + n);
                }

                //Pick proposal type
                //int proposalType = pStep.getPropStep();
                int proposalType = 1;

                //Standard MH update for parameters in \theta
                if (proposalType == 1) {
                    
                    //Propose and accept/reject new \theta* with MH step
                    thetaProposalCount++;
                    DoubleArrayList thetaEstNew = thetaProposal.nextProposal(thetaEstNow, epi);
                    //DoubleArrayList thetaEstNew = thetaEstNow.copy(); //if not proposing new params
                    DoubleArrayList thetaAllNew = thetaEstNew.copy();
                    thetaAllNew.addAllOf(thetaFixed);
                    epi.updateEstParams(thetaEstNew); //This isn't actually necessary since the particles get update automatically
                    
                    filter.runFilterForwardBack(particles, epi, dataZNow, treeStart, coal, startTime, endTime, stndNorm, factory2D, factory3D);
                    double pGNew = filter.margLikelihood;
                    
                    if (verbose) {
                        System.out.println("thetaNow: " + thetaEstNow);
                        System.out.println("thetaNew " + thetaEstNew);
                        System.out.println("pGNow: " + pGNow);
                        System.out.println("pGNew: " + pGNew);
                    }
                    
                    double probAcceptTheta = Math.exp(pGNew - pGNow);
                    if (probAcceptTheta >= Math.random()) {
                        thetaAllNow = thetaAllNew.copy(); //check this
                        thetaEstNow = thetaEstNew.copy();
                        pGNow = pGNew;
                        xTrajNow = filter.xTrajSample;
                        if (verbose) {
                            System.out.println("Accepted theta proposal");
                        }
                        thetaAcceptCount++;

                    } else {
                        epi.updateEstParams(thetaEstNow);
                    }

                }

                if (paramCounter == paramThin) {

                    sampleWriter.samplesToString(thetaEstNow);
                    thetaPrintWriter.print(sampleWriter.sampleString); thetaPrintWriter.flush(); thetaPrintWriter.println();

                    //totalLikeNow = pGThetaNow + pDGPsiNow;
                    likeSamples = new DoubleArrayList();
                    likeSamples.add(pGNow); //likeSamples.add(pXNow); likeSamples.add(pLgNow); //likeSamples.set(2,totalLikeNow);
                    sampleWriter.samplesToString(likeSamples);
                    likePrintWriter.print(n + ", "); likePrintWriter.print(sampleWriter.sampleString); likePrintWriter.flush(); likePrintWriter.println();

                    //Get xTraj samples
//                    xTrajSampleNow = xTrajs.getSample(sampleVar1, dataZNow, xTrajNow);
//                    sampleWriter.samplesToString(xTrajSampleNow);
//                    xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
//                    xTrajSampleNow = xTrajs.getSample(sampleVar2, dataZNow, xTrajNow);
//                    sampleWriter.samplesToString(xTrajSampleNow);
//                    xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
                    
                    //For incidence 
                    xTrajSampleNow = xTrajs.getSampleIncidence(sampleVar1, dataZNow, xTrajNow);
                    sampleWriter.samplesToString(xTrajSampleNow);
                    xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
                    xTrajSampleNow = xTrajs.getSampleIncidence(sampleVar2, dataZNow, xTrajNow);
                    sampleWriter.samplesToString(xTrajSampleNow);
                    xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
                    
                    paramCounter = 0;
                }

                if (optimizer.autoOptimize) {   
                    optimizer.addSample(thetaEstNow);
                    if (optimizer.runPoints.contains(n)) {
                        DoubleMatrix2D newCovMatrix = optimizer.autoOptimize(); 
                        thetaProposal.covMatrix = newCovMatrix;
                    }
                }
                
                if (verbose) {
                    System.out.println();
                }

        }
        //MCMC statistics
        double elapsedClockTime = System.currentTimeMillis() - startClockTime;
        double elapsedTimeSecs = elapsedClockTime / 1000;
        System.out.println("Total run time was:" + elapsedTimeSecs + "seconds");
        double iters = (double) iterations;
        double avgTimePerIteration = elapsedTimeSecs/ iters;
        System.out.println("Average time per iteration: " + avgTimePerIteration + " secs");

        //System.out.println("Tree acceptance rate:");
        //double accepts = (double) treeAcceptCount;
        //double its = (double) treeProposalCount;
        //double acceptRate = accepts / its;
        //System.out.println(acceptRate);

        System.out.println("Theta acceptance rate:");
        double accepts = (double) thetaAcceptCount;
        double its = (double) thetaProposalCount;
        double acceptRate = accepts / its;
        System.out.println(acceptRate);
        
        WriteMatrixToCSV covWriter = new WriteMatrixToCSV();
        String covFileName = fileNameStem + "_proposalMatrix_" + fileNamePostFix;
        covWriter.writeCSV(thetaProposal.covMatrix, covFileName);

        treePrintWriter.close();
        treeFileWriter.close();
        thetaPrintWriter.close();
        thetaFileWriter.close();
        likePrintWriter.close();
        likeFileWriter.close();
        xTrajPrintWriter.close();
        xTrajFileWriter.close();

    }//END method

}
