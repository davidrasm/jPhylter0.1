import cern.colt.list.*;
import java.util.ArrayList;

public class TreeUtils
{
	public static void getTreeInfo(TreeNode[] tree)
	{
		int nodeCount = tree.length;
		for (int k = 0; k < nodeCount; ++k) {
			System.out.println("Node number is");
			System.out.println(tree[k].nodeNumber);
			//System.out.println("Node label is:");
			//System.out.println(tree[k].nodeLabel);
			System.out.println("Node name is:");
			System.out.println(tree[k].nodeName);
			System.out.println("Node type is:");
			System.out.println(tree[k].nodeType);
			System.out.println("Node parent is:");
			System.out.println(tree[k].parentNode);
			System.out.println("Distance to parent is:");
			System.out.println(tree[k].distanceToParent);
			System.out.println("Children are:");
			System.out.println(tree[k].childNodes[0]);
			System.out.println(tree[k].childNodes[1]);
			System.out.println("Node height is:");
			System.out.println(tree[k].nodeHeight);
                        System.out.println();
		}
	}//END method
        
        /**
         * 
         * Should rewrite these methods to take advantage of tips then internal node order
         *
         */
	
	public static TreeNode[] parseTipTimes(TreeNode[] tree)
	{
		String delim = "_";
		int nodeCount = tree.length;
		for (int k = 1; k <= nodeCount; ++k) {
			if (tree[k-1].nodeType.equals("tip")) {
				String tipName = tree[k-1].nodeName;
				String[] na = tipName.split(delim);
				double tipTime = Double.parseDouble(na[1]);
				tree[k-1].nodeHeight = tipTime;
			}	
		}
		return tree;
	}//END method
	
	public static TreeNode[] parseTipLabels(TreeNode[] tree)
	{
		String delim = "-";
		int nodeCount = tree.length;
		for (int k = 1; k <= nodeCount; ++k) {
			if (tree[k-1].nodeType.equals("tip")) {
				String tipName = tree[k-1].nodeName;
				String[] na = tipName.split(delim);
				String labelStr = na[0];
				int label = Integer.parseInt(labelStr);
				tree[k-1].nodeLabel = label;
			}	
		}
		return tree;
	}//END method
	
	public static DoubleArrayList getTipTimes(TreeNode[] tree)
	{
		DoubleArrayList tipTimes = new DoubleArrayList();
		int nodeCount = tree.length;
		for (int k = 1; k <= nodeCount; ++k) {
			if (tree[k-1].nodeType.equals("tip")) {
				double tipTime = tree[k-1].nodeHeight;
				tipTimes.add(tipTime);
			}	
		}
		return tipTimes;
	}//END method
	
		
	public static ArrayList<Integer> getTipLocs(TreeNode[] tree)
	{
		ArrayList<Integer> tipLocs = new ArrayList<Integer>();
		int nodeCount = tree.length;
		for (int k = 1; k <= nodeCount; ++k) {
			if (tree[k-1].nodeType.equals("tip")) {
				tipLocs.add(k-1);
			}	
		}
		return tipLocs;
	}//END method
	
	
	
	public static DoubleArrayList getInternalNodeTimes(TreeNode[] tree)
	{
		//String delim = "_";
		DoubleArrayList internalTimes = new DoubleArrayList();
		int nodeCount = tree.length;
		for (int k = 1; k <= nodeCount; ++k) {
			if (tree[k-1].nodeType.equals("internal")) {
				double nodeHeight = tree[k-1].nodeHeight;
				internalTimes.add(nodeHeight);
			}
			if (tree[k-1].nodeType.equals("root")) {
				double nodeHeight = tree[k-1].nodeHeight;
				internalTimes.add(nodeHeight);
			}	
		}
		return internalTimes;
	}//END method
	
	public static ArrayList<Integer> getInternalNodeLocs(TreeNode[] tree)
	{
		//String delim = "_";
		ArrayList<Integer> internalLocs = new ArrayList<Integer>();
		int nodeCount = tree.length;
		for (int k = 1; k <= nodeCount; ++k) {
			if (tree[k-1].nodeType.equals("internal")) {
				internalLocs.add(k-1);
			}
			if (tree[k-1].nodeType.equals("root")) {
				internalLocs.add(k-1);
			}	
		}
		return internalLocs;
	}//END method
	
	
	public static TreeNode[] setNodeHeights(TreeNode[] tree)
	{
		//Written for nodeHeights in timeSincePastEvent (i.e. forwards time)
		//Modify if want nodeHeights in timeSincePresent (i.e. backwards time)
		int nodeCount = tree.length;
		//double tipHeight = 0.0;
		for (int k = 1; k <= nodeCount; ++k) {
			if (tree[k-1].nodeType.equals("tip")) {
				//Node height equals tip time
			} else {
				double tipHeight = 0.0;
				int childTip = 0;
				int l = k;
				//Move down tree until get to a tip
				while(childTip == 0) {
					l = tree[l-1].childNodes[0]; //always take left child node
					if (tree[l-1].nodeType.equals("tip")) {
						childTip = l;
						tipHeight = tree[l-1].nodeHeight;
					}
				}
				//Compute nodeHeight
				double heightNow = tipHeight;
				//Starting at tip, move back up tree to parent
				int currNode = l;
				while(currNode != k) {
					heightNow = heightNow - tree[currNode-1].distanceToParent; //if timeSincePastEvent
					//heightNow = heightNow + tree{currNode-1].distanceToParent; //if timeSincePresent
					currNode = tree[currNode-1].parentNode;
				}
				tree[k-1].nodeHeight = heightNow;
			}
		}
		return tree;
	}//END method
        
        public static TreeNode[] shiftNodeHeights(TreeNode[] tree, double amount)
	{
		int nodeCount = tree.length;
		for (int k = 0; k < nodeCount; k++) {
                    tree[k].nodeHeight += amount;
		}
		return tree;
	}//END method
        
        public static TreeNode[] setTipPriors(TreeNode[] tree, ArrayList<ArrayList<String>> priors)
	{
                int tipCount = priors.size();
                int states = priors.get(0).size() - 1;
                ArrayList<String> tipNames = new ArrayList<String>();
                for (int row = 0; row < tipCount; row++) {
                    tipNames.add(priors.get(row).get(0));
                }
                
                //For each tip in tree, look up priors
                for (int k = 0; k < tree.length; k++) {
                    if (tree[k].nodeType.equals("tip")) {
                        int indexInPriors = tipNames.indexOf(tree[k].nodeName);
                        if (indexInPriors == -1) {
                            System.out.println("Could not find priors for leaf: " + k);
                        } else {
                            for (int s = 0; s < states; s++) {
                                tree[k].statePriors.add(Double.parseDouble(priors.get(indexInPriors).get(s+1)));
                            }
                            
                        }
                    }
                }
                
                //System.out.println();
                return tree;

	}//END method
        
}