import java.util.ArrayList;

/**
 * Sample from a multinomial density
 * @author David
 */
public class Multinomial{
    
    public int sample(ArrayList<Double> probs) //let probs be a vector holding the probability of each event or type
    {
        //Sample from a multinomial density, elements in probs vector should already sum to one
        
        int sample = 0; //sample event or type from multinomial
        ArrayList<Double> cumulativeProbs = new ArrayList<Double>();
	
        //Find cumulative probs
        double currSum = 0.0;
	for (int j = 0; j < probs.size(); ++j) {
            //wNorm.set(j,(wCurr.getQuick(j)/wSum)); //normalized particle weights
            currSum += cumulativeProbs.get(j);
            cumulativeProbs.add(currSum);
	}	
	
        //Draw a random sample
	for (int n = 0; n < probs.size(); ++n)
	{
            double uniRand = Math.random();
            double cumVal = 0.0; 
            int counter = 0;
            while (cumVal <= uniRand)
            {
                cumVal = cumulativeProbs.get(counter);
		counter++;
            }
            sample = counter - 1; //-1 because already incremented counter
            
	}
	
        return sample;
        
    }//END method
    
}
