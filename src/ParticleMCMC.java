import cern.colt.list.*;
import cern.colt.matrix.*;
import java.util.ArrayList;
import java.io.*; //IOException;
//import cern.jet.stat.Descriptive;
//import cern.jet.random.*;

/**
 * Runs pMCMC algorithms
 * @author David
 */

public class ParticleMCMC
{
    
        //For unstructured models
        public void runPMMH(int iterations, int particles, int outputFreq, boolean verbose, TreeNode[] tree, double startTime, double endTime, double dt, ProposalDist thetaProposal, EpiModel epi, CoalModel coal, TimeSeries timeSeries, String fileNameStem, String fileNamePostFix) throws IOException
        {

                //Set up for writing MCMC samples to files:
            
                //Parameter sample output file
                SampleWriter sampleWriter = new SampleWriter();
                String paramFileName = fileNameStem + "_params_" + fileNamePostFix;
                FileWriter thetaFileWriter = new FileWriter(paramFileName);
                PrintWriter thetaPrintWriter = new PrintWriter(thetaFileWriter);

                //Pop state trajectory output file
                String xTrajFileName = fileNameStem + "_xTrajs_" + fileNamePostFix;
                FileWriter xTrajFileWriter = new FileWriter(xTrajFileName);
                PrintWriter xTrajPrintWriter = new PrintWriter(xTrajFileWriter);

                //Likelihood output file
                String likeFileName = fileNameStem + "_likes_" + fileNamePostFix;
                FileWriter likeFileWriter = new FileWriter(likeFileName);
                PrintWriter likePrintWriter = new PrintWriter(likeFileWriter);
                
                //Write initial params to output file
                DoubleArrayList thetaEstNow = epi.getEstParams();
                DoubleArrayList thetaAllNow = thetaEstNow.copy();
                DoubleArrayList thetaFixed = epi.getFixedParams();
                thetaAllNow.addAllOf(thetaFixed);
                sampleWriter.samplesToString(thetaEstNow);
                thetaPrintWriter.print(sampleWriter.sampleString); thetaPrintWriter.flush(); thetaPrintWriter.println();

                //Set initial ZVectors - this is the main data structure used to compute the coalescent likelihood
                ZVectors dataZNow = new ZVectors(); //Vector data structure for tree and observation data
                if (timeSeries.active) {
                    dataZNow.setZVectorsWithTimeSeries(tree, startTime, endTime, dt, timeSeries);
                } else {
                    dataZNow.setZVectors(tree, startTime, endTime, dt);
                }
                startTime = dataZNow.startTime;
                endTime = dataZNow.endTime;

                //Run particle filter to get initial xTraj and likelihood
                ParticleFilter filter = new ParticleFilter();
                filter.runFilter(particles, epi, dataZNow, coal, startTime, endTime, timeSeries); 
                double pGNow = filter.margLikelihood;//returns log marginal likelihood estimate pGNow
                double pThetaNow = epi.getPriorProb(); //get log prior on theta pThetaNow
                XTrajectory xTrajNow = filter.xTrajSample;

                //Initialize xTrajSamples for xTraj samples
                XTrajSamples xTrajs = new XTrajSamples();
                double sampleFreq = 5.0; //how often to sample pop state variables
                int sampleVar1 = 0; //set variable to sample
                //int sampleVar2 = 2; //set variable to sample
                xTrajs.getSampleTimes(startTime, endTime, sampleFreq);
                sampleWriter.samplesToString(xTrajs.sampleTimes);
                xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
                //DoubleArrayList xTrajSampleNow = xTrajs.getSample(sampleVar1, dataZNow, xTrajNow);
                //sampleWriter.samplesToString(xTrajSampleNow);
                //xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
                //xTrajSampleNow = xTrajs.getSample(sampleVar2, dataZNow, xTrajNow);
                //sampleWriter.samplesToString(xTrajSampleNow);
                //xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();

                //Write out starting likelihoods and priors
                DoubleArrayList likeSamples = new DoubleArrayList();
                likeSamples.add(pGNow); //likeSamples.add(pThetaNow);
                sampleWriter.samplesToString(likeSamples);
                likePrintWriter.print(0 + ", "); likePrintWriter.print(sampleWriter.sampleString); likePrintWriter.flush(); likePrintWriter.println();

                //Set up autoOptimization routines for parameter proposal density
                AutoOptimizer optimizer = new AutoOptimizer();
                optimizer.autoOptimize = true;
                if (optimizer.autoOptimize) {
                    ArrayList<Integer> runPoints = new ArrayList<Integer>();
                    runPoints.add(2000); runPoints.add(5000); runPoints.add(10000); runPoints.add(25000); runPoints.add(50000);
                    optimizer.setRunPoints(runPoints, thetaEstNow.size());
                    optimizer.getSampleMatrix();
                }
                
                //Initialize MCMC timer and counters
                double startClockTime = System.currentTimeMillis();
                int thetaAcceptCount = 0;
                int thinCounter = 0; //counts iterations to next sample

                //Run MCMC
                for (int n = 1; n <= iterations; ++n) {

                        thinCounter++;
                        if (verbose) {
                            System.out.println("MCMC iteration = " + n);
                        }

                        XTrajectory xTrajBacked = filter.xTrajSample.copy(); //should probably deep copy this!!
                        
                        //Propose new theta and run filter to evaluate marginal likelihood
                        DoubleArrayList thetaEstNew = thetaProposal.nextProposal(thetaEstNow, epi);
                        DoubleArrayList thetaAllNew = thetaEstNew.copy();
                        thetaAllNew.addAllOf(thetaFixed);
                        epi.updateEstParams(thetaEstNew);

                        //Run particle filter to get new marginal likelihood estimate
                        filter.runFilter(particles, epi, dataZNow, coal, startTime, endTime, timeSeries);
                        double pGNew = filter.margLikelihood;
                        double pThetaNew = epi.getPriorProb();

                        if (verbose) {
                            System.out.println("thetaEstNow: " + thetaEstNow);
                            System.out.println("thetaEstNew " + thetaEstNew);
                            System.out.println("posteriorProbNow: " + (pGNow + pThetaNow));
                            System.out.println("posteriorProbNew " + (pGNew + pThetaNew));
                        }

                        double probAcceptTheta = Math.exp((pGNew + pThetaNew) - (pGNow + pThetaNow)); //exp b/c likelihood and priors are in log values
                        if (probAcceptTheta >= Math.random()) {
                            thetaAllNow = thetaAllNew.copy(); //check this
                            thetaEstNow = thetaEstNew.copy();
                            pGNow = pGNew;
                            pThetaNow = pThetaNew;
                            if (pGNow > -363.0) {
                                pGNow = -363.0;
                            }
                            if (verbose) {
                                System.out.println("Accepted theta proposal");
                            }
                            thetaAcceptCount++;

                        } else {
                            epi.updateEstParams(thetaEstNow);
                            filter.xTrajSample = xTrajBacked; //reset xTraj to old sample
                        }
                        xTrajNow = filter.xTrajSample;
                            
                        //Print out MCMC samples
                        if (thinCounter == outputFreq) {
                            
                            System.out.println("MCMC iteration = " + n);
                            
                            //Write parameter sample to file
                            sampleWriter.samplesToString(thetaEstNow);
                            thetaPrintWriter.print(sampleWriter.sampleString); thetaPrintWriter.flush(); thetaPrintWriter.println();

                            //Write likelihood to file
                            likeSamples = new DoubleArrayList();
                            likeSamples.add(pGNow); //likeSamples.add(pThetaNow);;
                            sampleWriter.samplesToString(likeSamples);
                            likePrintWriter.print(n + ", "); likePrintWriter.print(sampleWriter.sampleString); likePrintWriter.flush(); likePrintWriter.println();
                            
                            //xTrajSampleNow = xTrajs.getSample(sampleVar1, dataZNow, xTrajNow);
                            //sampleWriter.samplesToString(xTrajSampleNow);
                            //xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
                            //xTrajSampleNow = xTrajs.getSample(sampleVar2, dataZNow, xTrajNow);
                            //sampleWriter.samplesToString(xTrajSampleNow);
                            //xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();

                            thinCounter = 0;
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
                
                //double likeVariance = Descriptive.sampleVariance(margLikes, Descriptive.mean(margLikes));
                //System.out.println("Likelihood variance = " + likeVariance);

                System.out.println("Theta acceptance rate:");
                double accepts = (double) thetaAcceptCount;
                double acceptRate = accepts / iters;
                System.out.println(acceptRate);
                
                //Write final proposal covariance matrix to file
                WriteMatrixToCSV covWriter = new WriteMatrixToCSV();
                String covFileName = fileNameStem + "_proposalMatrix_" + fileNamePostFix;
                covWriter.writeCSV(thetaProposal.covMatrix, covFileName);
                
                thetaPrintWriter.close();
                thetaFileWriter.close();
                likePrintWriter.close();
                likeFileWriter.close();
                xTrajPrintWriter.close();
                xTrajFileWriter.close();


        }//END method

    
        //For structured models
        public void runStructPMMH(int iterations, int particles, int outputFreq, boolean verbose, TreeNode[] tree, double startTime, double endTime, double dt, ProposalDist thetaProposal, EpiModel epi, StructCoalModel coal, String fileNameStem, String fileNamePostFix) throws IOException
        {

                //Set up for writing MCMC samples to files:
            
                //Parameter sample output file
                SampleWriter sampleWriter = new SampleWriter();
                String paramFileName = fileNameStem + "_params_" + fileNamePostFix;
                FileWriter thetaFileWriter = new FileWriter(paramFileName);
                PrintWriter thetaPrintWriter = new PrintWriter(thetaFileWriter);

                //Pop state trajectory output file
                String xTrajFileName = fileNameStem + "_xTrajs_" + fileNamePostFix;
                FileWriter xTrajFileWriter = new FileWriter(xTrajFileName);
                PrintWriter xTrajPrintWriter = new PrintWriter(xTrajFileWriter);

                //Likelihood output file
                String likeFileName = fileNameStem + "_likes_" + fileNamePostFix;
                FileWriter likeFileWriter = new FileWriter(likeFileName);
                PrintWriter likePrintWriter = new PrintWriter(likeFileWriter);
                
                //Write initial params to output file
                DoubleArrayList thetaEstNow = epi.getEstParams();
                DoubleArrayList thetaAllNow = thetaEstNow.copy();
                DoubleArrayList thetaFixed = epi.getFixedParams();
                thetaAllNow.addAllOf(thetaFixed);
                sampleWriter.samplesToString(thetaEstNow);
                thetaPrintWriter.print(sampleWriter.sampleString); thetaPrintWriter.flush(); thetaPrintWriter.println();

                //Set initial ZVectors - this is the main data structure used to compute the coalescent likelihood
                ZVectors dataZNow = new ZVectors(); //Vector data structure for tree and observation data
                dataZNow.setZVectors(tree, startTime, endTime, dt);
                startTime = dataZNow.startTime;
                endTime = dataZNow.endTime;

                //Run particle filter to get initial xTraj and likelihood
                ParticleImportanceSampler filter = new ParticleImportanceSampler();
                filter.setUpFilter(dataZNow, tree, coal, startTime, endTime);
                double pGNow = filter.runPIS(particles, epi, dataZNow, tree, coal, startTime, endTime); //returns log marginal likelihood estimate pGNow
                double pThetaNow = epi.getPriorProb(); //get log prior on theta pThetaNow
                XTrajectory xTrajNow = filter.xTrajSample;

                //Initialize xTrajSamples for xTraj samples
                XTrajSamples xTrajs = new XTrajSamples();
                double sampleFreq = 28.0; //how often to sample pop state variables
                int sampleVar1 = 1; //set variable to sample
                int sampleVar2 = 2; //set variable to sample
                xTrajs.getSampleTimes(startTime, endTime, sampleFreq);
                sampleWriter.samplesToString(xTrajs.sampleTimes);
                xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
                DoubleArrayList xTrajSampleNow = xTrajs.getSample(sampleVar1, dataZNow, xTrajNow);
                sampleWriter.samplesToString(xTrajSampleNow);
                xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
                xTrajSampleNow = xTrajs.getSample(sampleVar2, dataZNow, xTrajNow);
                sampleWriter.samplesToString(xTrajSampleNow);
                xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();

                //Write out starting likelihoods and priors
                DoubleArrayList likeSamples = new DoubleArrayList();
                likeSamples.add(pGNow); //likeSamples.add(pThetaNow);
                sampleWriter.samplesToString(likeSamples);
                likePrintWriter.print(0 + ", "); likePrintWriter.print(sampleWriter.sampleString); likePrintWriter.flush(); likePrintWriter.println();

                //Set up autoOptimization routines for parameter proposal density
                AutoOptimizer optimizer = new AutoOptimizer();
                optimizer.autoOptimize = true;
                if (optimizer.autoOptimize) {
                    ArrayList<Integer> runPoints = new ArrayList<Integer>();
                    runPoints.add(200); runPoints.add(500); runPoints.add(1000); runPoints.add(2000); runPoints.add(5000);
                    optimizer.setRunPoints(runPoints, thetaEstNow.size());
                    optimizer.getSampleMatrix();
                }
                
                //Initialize MCMC timer and counters
                double startClockTime = System.currentTimeMillis();
                int thetaAcceptCount = 0;
                int thinCounter = 0; //counts iterations to next sample

                //Run MCMC
                for (int n = 1; n <= iterations; ++n) {

                        thinCounter++;
                        if (verbose) {
                            System.out.println("MCMC iteration = " + n);
                        }

                        XTrajectory xTrajBacked = filter.xTrajSample.copy(); //should probably deep copy this!!
                        
                        //Propose new theta and run filter to evaluate marginal likelihood
                        DoubleArrayList thetaEstNew = thetaProposal.nextProposal(thetaEstNow, epi);
                        DoubleArrayList thetaAllNew = thetaEstNew.copy();
                        thetaAllNew.addAllOf(thetaFixed);
                        epi.updateEstParams(thetaEstNew);

                        //Run particle filter to get new marginal likelihood estimate    
                        double pGNew = filter.runPIS(particles, epi, dataZNow, tree, coal, startTime, endTime);
                        double pThetaNew = epi.getPriorProb();

                        if (verbose) {
                            System.out.println("thetaEstNow: " + thetaEstNow);
                            System.out.println("thetaEstNew " + thetaEstNew);
                            System.out.println("posteriorProbNow: " + (pGNow + pThetaNow));
                            System.out.println("posteriorProbNew " + (pGNew + pThetaNew));
                        }

                        double probAcceptTheta = Math.exp((pGNew + pThetaNew) - (pGNow + pThetaNow)); //exp b/c likelihood and priors are in log values
                        if (probAcceptTheta >= Math.random()) {
                            thetaAllNow = thetaAllNew.copy(); //check this
                            thetaEstNow = thetaEstNew.copy();
                            pGNow = pGNew;
                            pThetaNow = pThetaNew;
                            if (verbose) {
                                System.out.println("Accepted theta proposal");
                            }
                            thetaAcceptCount++;

                        } else {
                            epi.updateEstParams(thetaEstNow);
                            filter.xTrajSample = xTrajBacked; //reset xTraj to old sample
                        }
                        xTrajNow = filter.xTrajSample;
                            

                        if (thinCounter == outputFreq) {
                            
                            System.out.println("MCMC iteration = " + n);
                            
                            //Write parameter sample to file
                            sampleWriter.samplesToString(thetaEstNow);
                            thetaPrintWriter.print(sampleWriter.sampleString); thetaPrintWriter.flush(); thetaPrintWriter.println();

                            //Write likelihood to file
                            likeSamples = new DoubleArrayList();
                            likeSamples.add(pGNow); //likeSamples.add(pThetaNow);;
                            sampleWriter.samplesToString(likeSamples);
                            likePrintWriter.print(n + ", "); likePrintWriter.print(sampleWriter.sampleString); likePrintWriter.flush(); likePrintWriter.println();
                            
                            xTrajSampleNow = xTrajs.getSample(sampleVar1, dataZNow, xTrajNow);
                            sampleWriter.samplesToString(xTrajSampleNow);
                            xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
                            xTrajSampleNow = xTrajs.getSample(sampleVar2, dataZNow, xTrajNow);
                            sampleWriter.samplesToString(xTrajSampleNow);
                            xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();

                            thinCounter = 0;
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
                
                //double likeVariance = Descriptive.sampleVariance(margLikes, Descriptive.mean(margLikes));
                //System.out.println("Likelihood variance = " + likeVariance);

                System.out.println("Theta acceptance rate:");
                double accepts = (double) thetaAcceptCount;
                double acceptRate = accepts / iters;
                System.out.println(acceptRate);
                
                //Write final proposal covariance matrix to file
                WriteMatrixToCSV covWriter = new WriteMatrixToCSV();
                String covFileName = fileNameStem + "_proposalMatrix_" + fileNamePostFix;
                covWriter.writeCSV(thetaProposal.covMatrix, covFileName);
                
                thetaPrintWriter.close();
                thetaFileWriter.close();
                likePrintWriter.close();
                likeFileWriter.close();
                xTrajPrintWriter.close();
                xTrajFileWriter.close();


        }//END method

}