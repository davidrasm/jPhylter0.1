import cern.colt.list.*;
import java.io.IOException;
import cern.jet.stat.Descriptive;
import java.util.ArrayList;

/**
 * runPHYLTer is the main method used to set all initial values and load data
 * @author David
 */
public class runPHYLter {
    
        
    	public static void main(String args[]) throws IOException
        {
            
        //Load genealogy - tip labels nead to be formated like name_sampleTime
        String treeFileName = "./example/HIV_DetroitMSM_tree1.tre";
        TreeNode[] tree = NewickReader.getNewickTree(treeFileName); //parse tree from Newick input
        //TreeNode[] tree = NewickReader.getSimMapTree(treeFileName); //parse tree from SimMap format
        
        //For HIV tree - need to shift node heights
        TreeUtils.shiftNodeHeights(tree, 723496.0); //so time is in days since 01-01-0000
        CSVReader reader = new CSVReader();
        ArrayList<ArrayList<String>> tipStates = reader.readMatrixFromCSV("./example/sidSampleStates4-12-forDavid.csv");
        TreeUtils.setTipPriors(tree, tipStates);
        
        // Check root time
        DoubleArrayList coalTimesStart = TreeUtils.getInternalNodeTimes(tree);
        double rootTime = Descriptive.min(coalTimesStart);
        System.out.println("Root time: " + rootTime);
        
        //Load time series if using one
        //String timeSeriesFileName = "NY_fluSDI_2009_2010.csv";
        TimeSeries timeSeries = new TimeSeries();
        timeSeries.active = false;
        //ArrayList<ArrayList<String>> epiData = reader.readMatrixFromCSV(timeSeriesFileName);
        //timeSeries.setTimeSeries(epiData);
        
        //Set output file names
        String outputFileStem = "HIVTest_tree1";
        String outputFilePostFix = "_run1";
	
	//Set start and end times and dt integration step
	double startTime = 721355; // start time for particle filter
        //double = Double.NEGATIVE_INFINITY; //If NEGATIVE_INFINITY, startTime is root time 
	double endTime = Double.POSITIVE_INFINITY; //If POSITIVE_INFINITY, endTime is terminal sampling time
	double dt = 7; //Integration dt time step
        
        //Specify epidemiological model
        EpiModel epi = new EpiModelSECA();
        epi.setInitParams();
        epi.setRandomGenerator(); //seed random number generator
        epi.rootTime = rootTime;
        
	
	//Specify proposal density for epi parameters
        ProposalDist thetaProposal = new ProposalDist(); //this is always multivariate normal!
        thetaProposal.setRandomGenerator(); //seed random number generator
        //double [][] thetaCov = {{0.0000000001,0,0,0,0},{0,0.0000000005,0,0,0},{0,0,0.0000000005,0,0},{0,0,0,0.00000000001,0},{0,0,0,0,10.0}}; //SECA model for HIV with I_init
        double [][] thetaCov = {{0.0000000001,0,0},{0,0.0000000005,0},{0,0,0.0000000005}}; //SECA model for HIV without alpha
        thetaProposal.setCov(thetaCov, 3); // set covariance matrix for proposals
        
        //If specifying proposal cov matrix from file
        //String proposalCovFile = "BDTest_tree1_proposalMatrix_originEst_run1";
        //thetaProposal.covMatrix = WriteMatrixToCSV.readMatrixFromCSV(proposalCovFile);
	
        //Set pMCMC control params
	int particles = 10; //number of particles for particle filtering
	int iterations = 20000; //number of MCMC iterations
        int outputFreq = 100; //frequency at which MCMC samples are written to file
        boolean verbose = true; //print MCMC output?
        
        //Run pMCMC (choose between structured or unstructured version)
        ParticleMCMC mcmc = new ParticleMCMC();
        boolean structured = true; //is coalescent model structured?
        if (structured) {
            
            /**
             * Run modified PMMH algorithm for structured models (Rasmussen et al., PLoS Comp Bio, 2014)
             */
            
            //Specify coalescent model (needs to match EpiModel)
            StructCoalModel coal = new StructCoalModelSECA(); //modify coal model to fit process model
            coal.make();
            
            mcmc.runStructPMMH(iterations, particles, outputFreq, verbose, tree, startTime, endTime, dt, thetaProposal, epi, coal, outputFileStem, outputFilePostFix);
            
        } else {
            
            /**
             * Run original PMMH algorithm (Rasmussen et al., PLoS Comp Bio, 2011)
             */
            
            //Specify coalescent model (needs to match EpiModel)
            CoalModel coal = new CoalModelExpGrowth(); //an unstructured model
            
            mcmc.runPMMH(iterations, particles, outputFreq, verbose, tree, startTime, endTime, dt, thetaProposal, epi, coal, timeSeries, outputFileStem, outputFilePostFix);
        }
        
        //For fitting deterministic models
        //ClassicMCMC mcmc = new ClassicMCMC();
	//mcmc.runMCMC(iterations, particles, startTree, alignment, model, startTime, endTime, dt, thetaEst, thetaFixed, thetaProposal, psiParamsEst, psiParamsFixed, psiProposal, piParams, pStep, epi);
        
        //For computing likelihood surfaces:
        //startTime = 0.0;
        //LikelihoodSurface likeSurf = new LikelihoodSurface();
        //likeSurf.getSurface(particles, startTree, epi, startTime, endTime, dt, thetaEst, thetaFixed);
        //EntropyTree eTree = new EntropyTree();
        //eTree.getTree(particles, startTree, epi, startTime, endTime, dt, thetaEst, thetaFixed);
        
          
    }
    
        
    //Batch version
    /**
    public static void main(String args[]) throws IOException
    {
            
        //Load genealogy - tip labels nead to be formated like name_sampleTime
        //System.out.println("Root time: " + rootTime);
        
        //Load time series
        //String timeSeriesFileName = "NY_fluSDI_2009_2010.csv";
        TimeSeries timeSeries = new TimeSeries();
        timeSeries.active = false;
        CSVReader reader = new CSVReader();
        //ArrayList<ArrayList<String>> epiData = reader.readMatrixFromCSV(timeSeriesFileName);
        //timeSeries.setTimeSeries(epiData);
        
        //Assign priors on tip states (only for structured coalescent models)
        //CSVReader reader = new CSVReader();
        //ArrayList<ArrayList<String>> tipStates = reader.readMatrixFromCSV("sidSampleStates4-12-forDavid.csv");
        //TreeUtils.setTipPriors(tree, tipStates);
	
	//Set start and end times and dt integration step
	double startTime = 0.0; //734015.0;
        //double = Double.NEGATIVE_INFINITY; //If NEGATIVE_INFINITY, startTime is root time 
	double endTime = Double.POSITIVE_INFINITY; //734134.0; //If POSITIVE_INFINITY, endTime is terminal sampling time
	double dt = 0.05; //Integration dt time step
        
        //Specify epidemiological model
        //EpiModel epi = new EpiModelSuperFlu();
        //EpiModel epi = new EpiModelExpGrowth();
        //epi.setInitParams();
        //epi.setRandomGenerator(); //seed random number generator
        //epi.rootTime = rootTime;
        
        //Specify coalescent model (needs to match EpiModel)
        boolean structured = false; //is coalescent model structured?
        //if (structured) {
            //StructCoalModel structCoal = new StructCoalModel2PopDengue(); //for vector-borne models
            //structCoal.make();
        //} else {
            CoalModel coal = new CoalModelExpGrowth();
        //}
	
	//Specify proposal density for epi parameters
        //ProposalDist thetaProposal = new ProposalDist(); //this is always multivariate normal!
        //thetaProposal.setRandomGenerator(); //seed random number generator
        //double [][] thetaCov = {{0.0000001}};
        //double [][] thetaCov = {{0.0000001,0.0},{0.0,0.000000001}};
        //double [][] thetaCov = {{0.0000001,0.0,0.0},{0.0,0.0000001,0.0},{0.0,0.0,0.05}};
        //double [][] thetaCov = {{0.0000001,0.0,0.0,0.0},{0.0,0.00000001,0.0,0.0},{0.0,0.0,0.00000001,0.0},{0.0,0.0,0.0,0.0000000001}};
        //double [][] thetaCov = {{0.0000001,0.0,0.0,0.0,0.0},{0.0,0.00000001,0.0,0.0,0.0},{0.0,0.0,0.0000000001,0.0,0.0},{0.0,0.0,0.0,0.0000000001,0.0},{0.0,0.0,0.0,0.0,0.0000000000001}};
        //double [][] thetaCov = {{0.0000001,0.0,0.0,0.0,0.0,0.0},{0.0,0.00000001,0.0,0.0,0.0,0.0},{0.0,0.0,0.0000000001,0.0,0.0,0.0},{0.0,0.0,0.0,0.0000000001,0.0,0.0},{0.0,0.0,0.0,0.0,0.0000000000001,0.0},{0.0,0.0,0.0,0.0,0.0,0.0000005}};
        //thetaProposal.setCov(thetaCov, 3);
        
        //If specifying proposal cov matrix from file
        //String proposalCovFile = "BDTest_tree1_proposalMatrix_originEst_run1";
        //thetaProposal.covMatrix = WriteMatrixToCSV.readMatrixFromCSV(proposalCovFile);
	
        //Set pMCMC control params
	int particles = 100; //number of particles for particle filtering
	int iterations = 300000; //number of MCMC iterations
        int outputFreq = 100; //frequency at which MCMC samples are written to file
        boolean verbose = false; //print MCMC output?
        
        int nSims = 100;
        ArrayList<Integer> skips = new ArrayList<Integer>();
        skips.add(1); skips.add(2);
        
        ArrayList<Integer> toDo = new ArrayList<Integer>();
        toDo.add(73);
        
        
        //Run pMCMC
        ParticleMCMC mcmc = new ParticleMCMC();
        if (structured) {
            //Run modified PMMH algorithm for structured models
            //mcmc.runStructPMMH(iterations, particles, outputFreq, verbose, tree, startTime, endTime, dt, thetaProposal, epi, coal, outputFileStem, outputFilePostFix);
        } else {
            //Run original PMMH algorithm (Rasmussen et al., PLoS Comp Bio, 2011)
            for (int sim = 1; sim < nSims; sim++) {
            
                int treeNum = sim + 1;
                
                //if(!skips.contains(treeNum)) {
                if(toDo.contains(treeNum)) {
                    
                    System.out.println("Starting tree number: " + treeNum);
                    
                    String treeFileName = "BDTest_tree" + Integer.toString(treeNum) + "Dated.tre";
                    TreeNode[] tree = NewickReader.getNewickTree(treeFileName); //parse tree from Newick input
        
                    TreeUtils.shiftNodeHeights(tree, 100.0); //shift so not at time zero
                    DoubleArrayList coalTimesStart = TreeUtils.getInternalNodeTimes(tree);
                    double rootTime = Descriptive.min(coalTimesStart);
                    
                    EpiModel epi = new EpiModelExpGrowth();
                    epi.setInitParams();
                    epi.setRandomGenerator(); //seed random number generator
                    epi.rootTime = rootTime;
                    
                    ProposalDist thetaProposal = new ProposalDist(); //this is always multivariate normal!
                    thetaProposal.setRandomGenerator(); //seed random number generator
                    double [][] thetaCov = {{0.005,0.0,0.0},{0.0,0.005,0.0},{0.0,0.0,0.05}};
                    thetaProposal.setCov(thetaCov, 3);
                    //String proposalCovFile = "BDTest_tree1_proposalMatrix_originEst_run1";
                    //thetaProposal.covMatrix = WriteMatrixToCSV.readMatrixFromCSV(proposalCovFile);
                    
                    //Set output file names
                    String outputFileStem = "BDTest_tree" + Integer.toString(treeNum);
                    String outputFilePostFix = "originEst_run1";
            
                    mcmc.runPMMH(iterations, particles, outputFreq, verbose, tree, startTime, endTime, dt, thetaProposal, epi, coal, timeSeries, outputFileStem, outputFilePostFix);
                    
                }
            }
        }
        
        //For fitting deterministic models
        //ClassicMCMC mcmc = new ClassicMCMC();
	//mcmc.runMCMC(iterations, particles, startTree, alignment, model, startTime, endTime, dt, thetaEst, thetaFixed, thetaProposal, psiParamsEst, psiParamsFixed, psiProposal, piParams, pStep, epi);
        
        //For computing likelihood surfaces:
        //startTime = 0.0;
        //LikelihoodSurface likeSurf = new LikelihoodSurface();
        //likeSurf.getSurface(particles, startTree, epi, startTime, endTime, dt, thetaEst, thetaFixed);
        //EntropyTree eTree = new EntropyTree();
        //eTree.getTree(particles, startTree, epi, startTime, endTime, dt, thetaEst, thetaFixed);
        
          
    }
    */
}
