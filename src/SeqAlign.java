import java.util.ArrayList;

public class SeqAlign
{
	int sequences[][];
	String labels[];
	int sites;
	int seqs;
        
        public void makeAlignmentNoGaps(String[] seqArray, String[] labelArray)
	{
		int totalSites = 1000; //seqArray[1].length();
                ArrayList<Integer> goodSites = new ArrayList<Integer>(); 
		seqs = seqArray.length;
		
		char siteChar; char siteNuc; int siteInt = 0;
		char aUpper = 'A'; char gUpper = 'G'; char cUpper = 'C'; char tUpper = 'T';
                char aLower = 'a'; char gLower = 'g'; char cLower = 'c'; char tLower = 't';
                char hyphen = '-';
                
                //First check all sites to make sure their legit
                boolean foundBadSite;
                for (int j = 0; j < totalSites; ++j) {
                    foundBadSite = false;
                    for (int i = 0; i < seqs; ++i) {
                        if (seqArray[i].charAt(j) == hyphen) {
                            foundBadSite = true;
                        }
                    }
                    if (foundBadSite == false) {
                        goodSites.add(j);
                    }   
                }
                
                sites = goodSites.size();
                System.out.println("Generating sequence alignment:");
                System.out.println("Found " + sites + " sites out of " + totalSites + " with acceptable site patterns");        
                sequences = new int[seqs][sites];
		labels = new String[seqs];
                
                int currSite;
		for (int i = 0; i < seqs; ++i) {
			labels[i] = labelArray[i];
			for (int j = 0; j < sites; ++j) {
                            currSite = goodSites.get(j);
				siteNuc = seqArray[i].charAt(currSite);
				if (siteNuc == aUpper | siteNuc == aLower) {
					siteInt = 1;
				}
				if (siteNuc == gUpper | siteNuc == gLower) {
					siteInt = 2;
				}
				if (siteNuc == cUpper | siteNuc == cLower) {
					siteInt = 3;
				}
				if (siteNuc == tUpper | siteNuc == tLower) {
					siteInt = 4;
				}
				sequences[i][j] = siteInt;
			}
		}
	}//END method

	public void makeAlignmentAllSites(String[] seqArray, String[] labelArray)
	{
		sites = seqArray[1].length();
		seqs = seqArray.length;
		sequences = new int[seqs][sites];
		labels = new String[seqs];
		
		char siteNuc; int siteInt = 0;
		char aUpper = 'A'; char gUpper = 'G'; char cUpper = 'C'; char tUpper = 'T';
                char aLower = 'a'; char gLower = 'g'; char cLower = 'c'; char tLower = 't';
		for (int i = 0; i < seqs; ++i) {
			labels[i] = labelArray[i];
			for (int j = 0; j < sites; ++j) {
				siteNuc = seqArray[i].charAt(j);
				if (siteNuc == aUpper | siteNuc == aLower) {
					siteInt = 1;
				}
				if (siteNuc == gUpper | siteNuc == gLower) {
					siteInt = 2;
				}
				if (siteNuc == cUpper | siteNuc == cLower) {
					siteInt = 3;
				}
				if (siteNuc == tUpper | siteNuc == tLower) {
					siteInt = 4;
				}
				sequences[i][j] = siteInt;
			}
		}
	}//END method
}//END class