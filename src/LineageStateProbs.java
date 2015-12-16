import cern.colt.list.*;
import cern.colt.matrix.*;
import java.util.ArrayList;

/**
 *
 * @author David
 */
public class LineageStateProbs {
    
    int particles;
    int states;
    int lineages;
    DoubleMatrix3D matrix;
    ArrayList<Integer> activeLineages = new ArrayList<Integer>();
    ArrayList<Integer> lineagesRemoved = new ArrayList<Integer>();
    ArrayList<Integer> lineagesAdded = new ArrayList<Integer>();
    
    public void getMatrix(int jParticles, int lines, int sts, DoubleFactory3D factory3D) {
        
        states = sts;
        particles = jParticles;
        lineages = lines;
        matrix = factory3D.make(particles, lineages, sts);
        
    }
    
    public void addCoalLineages(ZVectors dataZ, int zLoc, TreeNode[] tree, DoubleMatrix2D xCurr, DoubleArrayList theta, StructCoalModel coal)
    {
    
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        double currTime = dataZ.absoluteTimes.get(zLoc);
        for (int event = 0; event < totalEvents; event++) {
            
            if (dataZ.omegaEvents.get(zLoc).get(event) == 1) { //if coalescent event
                
                int coalLineage = dataZ.nodePointers.get(zLoc).get(event);
                int daughterLineage1 = tree[coalLineage].childNodes[0] - 1;
                int daughterLineage2 = tree[coalLineage].childNodes[1] - 1;
                //Add daughter lineages to activeLineages and remove parent
                activeLineages.add(daughterLineage1);
                activeLineages.add(daughterLineage2);
                int listIndexCoalLineage = activeLineages.indexOf(coalLineage);
                if (listIndexCoalLineage == -1) {
                    //DoubleArrayList coalTimes = TreeUtils.getInternalNodeTimes(tree);
                    System.out.println("Cannot remove lineage");
                } else {
                    activeLineages.remove(listIndexCoalLineage);
                }
                
                for (int x = 0; x < particles; x++) {
                
                    coal.updateF(xCurr.viewColumn(x), theta, currTime);
                    coal.updateG(xCurr.viewColumn(x), theta);
                    coal.updateY(xCurr.viewColumn(x));
                    double totalRateSum = 0.0;
                    for (int stateK = 0; stateK < states; stateK++) {
                        for (int stateL = 0; stateL < states; stateL++) {
                            double thisRate = (coal.F.get(stateL, stateK) / coal.Y.get(stateL)) * matrix.get(x, coalLineage, stateL);
                            totalRateSum += thisRate; 
                        }
                    }

                    for (int stateK = 0; stateK < states; stateK++) {
                        
                        double nonTransProbInStateK = 0.0;
                        for (int stateL = 0; stateL < states; stateL++) {
                            double thisRate = ((coal.F.get(stateL, stateK) / coal.Y.get(stateL)) * matrix.get(x, coalLineage, stateL));
                            nonTransProbInStateK += thisRate;
                        }
                        
                        double transProbInStateK = 0.0;
                        for (int stateL = 0; stateL < states; stateL++) {
                            double thisRate = ((coal.F.get(stateK, stateL) / coal.Y.get(stateK)) * matrix.get(x, coalLineage, stateK));
                            transProbInStateK += thisRate;
                        }
                        
                        double probStateK = (nonTransProbInStateK + transProbInStateK) / (2*totalRateSum);
                        
                        //Assign probabilities to daughter lineages
                        matrix.set(x, daughterLineage1, stateK, probStateK);
                        matrix.set(x, daughterLineage2, stateK, probStateK);
                    }
                    //System.out.println();
                } //END FOR particles    
            }
        }    
    }
        
    public void removeSamples(ZVectors dataZ, int zLoc) 
    {
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        double eventTime = dataZ.absoluteTimes.get(zLoc);
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                
                int sampledLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                int sampleListIndex = activeLineages.indexOf(sampledLineage);
                activeLineages.remove(sampleListIndex);
                
            }
        }      
    }
    
    public void addSamples(ZVectors dataZ, int zLoc, TreeNode[] tree) 
    {
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        double eventTime = dataZ.absoluteTimes.get(zLoc);
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                
                int newLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                activeLineages.add(newLineage);
                lineagesAdded.add(newLineage);
                int newState = tree[newLineage].nodeLabel - 1; //nodeLabels are indexed from one
                for (int j = 0; j < particles; j++) {
                    if (newState == -1) {
                        System.out.println("Sampled lineage not in defined state");
                    }
                    matrix.set(j, newLineage, newState, 1.0); //set new state prob to 1.0
                }
            }
        }
        
    }
    
    public void addSamplesWithPriors(ZVectors dataZ, int zLoc, TreeNode[] tree) 
    {
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        double statePrior = 0.0;
        //double eventTime = dataZ.absoluteTimes.get(zLoc);
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                
                int newLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                activeLineages.add(newLineage);
                lineagesAdded.add(newLineage);
                for (int j = 0; j < particles; j++) {
                    for (int s = 0; s < states; s++) {
                        statePrior = tree[newLineage].statePriors.get(s);
                        matrix.setQuick(j, newLineage, s, statePrior); //set new state prob to 1.0
                    }
                }
            }
        }
        
    }
    
    public void updateProbs(DoubleMatrix2D xCurr, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime) 
    {
        
        int numLineages = activeLineages.size();
        
        //For each particle
        for (int x = 0; x < particles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j
        
            coal.updateF(xCurr.viewColumn(x), theta, currTime);
            coal.updateG(xCurr.viewColumn(x), theta);
            coal.updateY(xCurr.viewColumn(x));
            
            //Precompute A(l) terms by summing probs in matrix
            DoubleArrayList lineagesInAByState = this.getLineagesByState(x);
            
            double omegaIn = 0.0; double omegaOut = 0.0;
            double omegaGain; double omegaLoss;
            
            //double lambdaSum = 0.0; //total rate of coalescence across all pairs of lineages
            for (int lin = 0; lin < numLineages; lin++) { //for all lineages i
                
                int i = activeLineages.get(lin);
                DoubleArrayList newLineageProbs = new DoubleArrayList();
                
                //Doing extra work here: what stateL loses to stateK is what stateK gains from stateL
                //Just have to compute gains for each state in one matrix and can compute loses from same matrix
                for (int k = 0; k < states; k++) {
                    
                    omegaGain = 0.0;
                    omegaLoss = 0.0;
                    
                    //Solve the master equation as would an ODE
                    for (int l = 0; l < states; l++) {    
                        if (l != k) {
                            omegaIn = ((coal.G.getQuick(k, l)/coal.Y.getQuick(l)) + (coal.F.getQuick(k, l)/coal.Y.getQuick(l))*((coal.Y.getQuick(k) - lineagesInAByState.getQuick(k))/coal.Y.getQuick(k))) * dtTime * matrix.get(x,i,l);
                            if (Double.isNaN(omegaIn)) {
                                //omegaGain = 0.0;
                            } else {
                                omegaGain += omegaIn;
                            }
                            omegaOut = ((coal.G.getQuick(l, k)/coal.Y.getQuick(k)) + (coal.F.getQuick(l, k)/coal.Y.getQuick(k))*((coal.Y.getQuick(l) - lineagesInAByState.getQuick(l))/coal.Y.getQuick(l))) * dtTime * matrix.get(x,i,k);
                            if (Double.isNaN(omegaOut)) {
                                //omegaLoss = 0.0;
                            } else {
                                omegaLoss += omegaOut;
                            }
                        }
                    }
                    
                    double newProbLinIStateK = matrix.get(x, i, k) + omegaGain - omegaLoss;
                    if (Double.isNaN(newProbLinIStateK)) {
                        //System.out.println();
                    }
                    newLineageProbs.add(newProbLinIStateK);
                }
            
                for (int k = 0; k < states; k++) {
                    matrix.set(x, i, k, newLineageProbs.get(k));
                }   
            }
        }
        //System.out.println();
    }
    
    public void updateProbsMarginal(DoubleMatrix1D xCurr, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime) 
    {
        
        int numLineages = activeLineages.size();
        
        //For each particle
        //for (int x = 0; x < particles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j
        
            coal.updateF(xCurr, theta, currTime);
            coal.updateG(xCurr, theta);
            coal.updateY(xCurr);
            
            //Precompute A(l) terms by summing probs in matrix
            DoubleArrayList lineagesInAByState = this.getLineagesByState(0);
            
            double omegaIn; double omegaOut;
            
            //double lambdaSum = 0.0; //total rate of coalescence across all pairs of lineages
            for (int lin = 0; lin < numLineages; lin++) { //for all lineages i
                
                int i = activeLineages.get(lin);
                DoubleArrayList newLineageProbs = new DoubleArrayList();
                
                //Doing extra work here: what stateL loses to stateK is what stateK gains from stateL
                //Just have to compute gains for each state in one matrix and can compute loses from same matrix
                for (int k = 0; k < states; k++) {
                    
                    double omegaGain = 0.0;
                    double omegaLoss = 0.0;
                    
                    //Solve the master equation as would an ODE
                    for (int l = 0; l < states; l++) {    
                        if (l != k) {
                            omegaIn = ((coal.G.getQuick(k, l)/coal.Y.getQuick(l)) + (coal.F.getQuick(k, l)/coal.Y.getQuick(l))*((coal.Y.getQuick(k) - lineagesInAByState.getQuick(k))/coal.Y.getQuick(k))) * dtTime * matrix.get(0,i,l);
                            //System.out.println("omegaGain = " + omegaGain);
                            if (Double.isNaN(omegaIn)) {
                                //omegaGain = 0.0;
                            } else {
                                omegaGain += omegaIn;
                            }
                            omegaOut = ((coal.G.getQuick(l, k)/coal.Y.getQuick(k)) + (coal.F.getQuick(l, k)/coal.Y.getQuick(k))*((coal.Y.getQuick(l) - lineagesInAByState.getQuick(l))/coal.Y.getQuick(l))) * dtTime * matrix.get(0,i,k);
                            //System.out.println("omegaLoss = " + omegaGain);
                            if (Double.isNaN(omegaOut)) {
                                //omegaLoss = 0.0;
                            } else {
                                omegaLoss += omegaOut;
                            }
                        }
                    }
                    
                    double newProbLinIStateK = matrix.get(0, i, k) + omegaGain - omegaLoss;
                    //System.out.println();
                    if (Double.isNaN(newProbLinIStateK)) {
                        //System.out.println();
                    }
                    newLineageProbs.add(newProbLinIStateK);
                }
            
                //System.out.println();
                for (int k = 0; k < states; k++) {
                    matrix.set(0, i, k, newLineageProbs.get(k));
                }   
            }
        //}
        
    }
    
    public void updateProbsFast(DoubleMatrix2D xCurr, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime) 
    {
        
        int numLineages = activeLineages.size();
        
        DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
        DoubleMatrix2D omegaGainBirths = factory2D.make(states, states);
        DoubleMatrix2D omegaGainMigrations = factory2D.make(states, states);
        
        //For each particle
        for (int x = 0; x < particles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j
        
            coal.updateF(xCurr.viewColumn(x), theta, currTime);
            coal.updateG(xCurr.viewColumn(x), theta);
            coal.updateY(xCurr.viewColumn(x));
            
            //Precompute A(l) terms by summing probs in matrix
            DoubleArrayList lineagesInAByState = this.getLineagesByState(x);
            
            //double lambdaSum = 0.0; //total rate of coalescence across all pairs of lineages
            for (int lin = 0; lin < numLineages; lin++) { //for all lineages i
                
                int i = activeLineages.get(lin);
                DoubleArrayList newLineageProbs = new DoubleArrayList();
                
                //Update omegaGainBirths
                double omegaGainKFromL = 0.0;
                for (int k = 0; k < states; k++) {
                    if (coal.Y.getQuick(k) <= 0.0) { //no gains are possible if k is <= 0                  
                        for (int l = 0; l < states; l++) {    
                            omegaGainBirths.setQuick(k, l, 0.0);
                        }   
                    } else {  
                        for (int l = 0; l < states; l++) {       
                            if (l == k | coal.Y.getQuick(l) <= 0.0) { //no gains are possible if l is <= 0
                                omegaGainBirths.setQuick(k, l, 0.0);
                            } else {
                                if (lineagesInAByState.getQuick(k) > coal.Y.getQuick(k)) {
                                    omegaGainKFromL = 0.0;
                                } else {
                                    omegaGainKFromL = ((coal.F.getQuick(k, l)/coal.Y.getQuick(l))*((coal.Y.getQuick(k) - lineagesInAByState.getQuick(k))/coal.Y.getQuick(k))) * dtTime * matrix.get(x,i,l);
                                }
                                omegaGainBirths.setQuick(k,l, omegaGainKFromL);
                            }  
                        }
                    }
                }
                
                //Update omegaGainMigration
                for (int k = 0; k < states; k++) {    
                    for (int l = 0; l < states; l++) {    
                        if (l == k | coal.Y.getQuick(l) <= 0.0) { //no gains are possible if l is <= 0
                            omegaGainMigrations.setQuick(k, l, 0.0);
                        } else {
                            omegaGainKFromL = (coal.G.getQuick(k, l)/coal.Y.getQuick(l)) * dtTime * matrix.get(x,i,l);
                            omegaGainMigrations.setQuick(k,l, omegaGainKFromL);
                        }
                    }  
                }
                
                double newProbLinIStateK = 0.0;
                boolean positivityCheckFail = false;
                for (int k = 0; k < states; k++) {
                    double omegaChange = 0.0;
                    for (int l = 0; l < states; l++) {
                        omegaChange += omegaGainBirths.getQuick(k,l) + omegaGainMigrations.getQuick(k,l) - omegaGainBirths.getQuick(l,k) - omegaGainMigrations.getQuick(l,k);
                    }
                    newProbLinIStateK = matrix.getQuick(x, i, k) + omegaChange;
                    if (newProbLinIStateK < 0) {
                        //System.out.println("WARNING: Lineage state probs returned negative after updating");
                        newProbLinIStateK = 0.0;
                        positivityCheckFail = true;
                    }
                    matrix.setQuick(x, i, k, newProbLinIStateK);
                }
                
                if (positivityCheckFail) {
                    this.renormalizeLineageProbs(x,i); //renormalize so lineage probs sum to one
                }

                //System.out.println();
            }
        }
    }
    
    public void updateProbsAfterCoal(int x, int coalNode, int daughterLineage1Index, int daughterLineage2Index, double lambdaSum, DoubleArrayList lambdaSumForEachK) 
    {
        
        for (int k = 0; k < states; k++) {
            double pThisState = lambdaSumForEachK.get(k) / lambdaSum;
            if (Double.isNaN(pThisState)) {
                //System.out.println("WARNING: Lineage state probs returned NaN after coalescent event");
                pThisState = 0.0;
            }
            if (pThisState < 0) {
                //System.out.println("WARNING: Lineage state probs returned negative after coalescence event");
                pThisState = 0.0;
            }
            matrix.set(x, coalNode, k, pThisState);
        }
        //System.out.println(); 
    }
    
    public DoubleArrayList getLineagesByState(int x) 
    {
        DoubleArrayList lineagesByState = new DoubleArrayList();
        int numLineages = activeLineages.size();
        for (int k = 0; k < states; k++) {
            double pSum = 0.0;
            for (int lin = 0; lin < numLineages; lin++) {
                int i = activeLineages.get(lin);
                double pMass = matrix.get(x, i, k);
                pSum += pMass;    
            }
            lineagesByState.add(pSum);
        }
        return lineagesByState;
    }
    
    private void renormalizeLineageProbs(int particle, int lineage) 
    {
        double newProb = 0.0;
        double probSum = 0.0;
        for (int k = 0; k < states; k++) {
            probSum += matrix.getQuick(particle, lineage, k);
        }
        
        for (int k = 0; k < states; k++) {
            newProb = matrix.getQuick(particle, lineage, k) / probSum;
            matrix.setQuick(particle, lineage, k, newProb);
        }
        
    }
    
    
    
    
}