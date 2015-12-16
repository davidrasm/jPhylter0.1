public class NewickReader
{
	public static TreeNode[] getNewickTree(String treeFileName)
	{
	
		String treeString = treeFileReader.readTree(treeFileName);
		System.out.println("New Netwick tree:");
		System.out.println(treeString);
	
		int treeStringLength = treeString.length();
	
		//Find number of internal nodes
		int intNodeCount = treeString.replaceAll("[^)]", "").length(); //includes the root node
 
		//Find number of tip nodes
		int tipNodeCount = treeString.replaceAll("[^,]", "").length() + 1;
		int nodeCount = intNodeCount + tipNodeCount;
	
		//Find number of edges
		int edgeCount = intNodeCount + tipNodeCount - 1;
		
		//Define characters
		char L_PAREN = '(';
		char COMMA = ',';
		char R_PAREN = ')';
		char COLON = ':';
		
		//Make tree object
		TreeNode[] tree = new TreeNode[nodeCount]; //change to ArrayList<TreeNode>
		for (int j = 0; j < nodeCount; ++j) {
			tree[j] = new TreeNode();
			tree[j].nodeNumber = j+1;
		}
		
		//Parse tree from newick format
		//Nodes 1:#Tips = tips
		//Nodes #Tips+1:#Nodes = internal nodes
		int currNode = tipNodeCount; //start bellow root
		int nextInternalNode = currNode + 1;
		int nextTipNode = 1;
		int stringIndex = 0;
		for (int i = 0; i < treeStringLength; ++i) {
			char currChar = treeString.charAt(stringIndex);
			
			if (currChar == L_PAREN) {
				//Create daughter node
				tree[nextInternalNode-1].setParent(currNode);
				tree[nextInternalNode-1].setNodeName("internal");
				tree[nextInternalNode-1].setNodeType("internal");
				//Move to daughter
				currNode = nextInternalNode; //4
				nextInternalNode = nextInternalNode + 1; //5
				tree[currNode-1].setChildren(0, nextInternalNode);
				if (currNode == (tipNodeCount + 1)) { // then its the root
					tree[currNode].setNodeType("root");
				}
				char nextChar = treeString.charAt(stringIndex+1);
				if (Character.isLetterOrDigit(nextChar)) {
					stringIndex = stringIndex + 1;
					//Get node name
					int indexLabelEnd = treeString.indexOf(COLON, (stringIndex));
					String tipLabel = treeString.substring(stringIndex, indexLabelEnd);
					stringIndex = indexLabelEnd;
					tree[nextTipNode-1].setParent(currNode);
					tree[nextTipNode-1].setNodeName(tipLabel);
					tree[nextTipNode-1].setNodeType("tip");
					tree[currNode-1].setChildren(0, nextTipNode);
					currNode = nextTipNode; //1
					nextTipNode = nextTipNode + 1; //2
					
					//Get edge length
					int indexNextParenth = treeString.indexOf(R_PAREN, stringIndex);
					int indexNextComma = treeString.indexOf(COMMA, stringIndex);
					int indexEdgeLengthEnd = 0;
					if (indexNextParenth < indexNextComma || indexNextComma < 0) {
						indexEdgeLengthEnd = indexNextParenth;
					} else {
						indexEdgeLengthEnd = indexNextComma;
					}
					String edgeLength = treeString.substring(stringIndex+1, indexEdgeLengthEnd);
					tree[currNode-1].setDistance(edgeLength);
					stringIndex = indexEdgeLengthEnd;
				} else {
					stringIndex = stringIndex + 1;
				}
			}
			
			if (currChar == COMMA) {
				//Go back to parent
				currNode = tree[currNode-1].getParent();
				char nextChar = treeString.charAt(stringIndex+1);
				if (Character.isLetterOrDigit(nextChar)) {
					stringIndex = stringIndex + 1;
					//Get node name
					int indexLabelEnd = treeString.indexOf(COLON, (stringIndex));
					String tipLabel = treeString.substring(stringIndex, indexLabelEnd);
					stringIndex = indexLabelEnd;
					tree[nextTipNode-1].setParent(currNode);
					tree[nextTipNode-1].setNodeName(tipLabel);
					tree[nextTipNode-1].setNodeType("tip");
					tree[currNode-1].setChildren(1, nextTipNode);
					currNode = nextTipNode;
					nextTipNode = nextTipNode + 1;
					
					//Get edge length
					int indexNextParenth = treeString.indexOf(R_PAREN, stringIndex);
					int indexNextComma = treeString.indexOf(COMMA, stringIndex);
					int indexEdgeLengthEnd = 0;
					if (indexNextParenth < indexNextComma || indexNextComma < 0) {
						indexEdgeLengthEnd = indexNextParenth;
					} else {
						indexEdgeLengthEnd = indexNextComma;
					}
					String edgeLength = treeString.substring(stringIndex+1, indexEdgeLengthEnd);
					tree[currNode-1].setDistance(edgeLength);
					stringIndex = indexEdgeLengthEnd;
				} else {
					tree[currNode-1].setChildren(1, nextInternalNode);
					stringIndex = stringIndex + 1;
				}
			}
			
			if (currChar == R_PAREN) {
				currNode = tree[currNode-1].getParent();
				char nextChar = treeString.charAt(stringIndex+1);
				if (nextChar == COLON) {
					//get distance
					stringIndex = stringIndex + 1; //now at COLON
					int indexNextParenth = treeString.indexOf(R_PAREN, stringIndex);
					int indexNextComma = treeString.indexOf(COMMA, stringIndex);
					int indexEdgeLengthEnd = 0;
					if (indexNextParenth < indexNextComma || indexNextComma < 0) {
						indexEdgeLengthEnd = indexNextParenth;
					} else {
						indexEdgeLengthEnd = indexNextComma;
					}
					String edgeLength = treeString.substring(stringIndex+1, indexEdgeLengthEnd);
					tree[currNode-1].setDistance(edgeLength);
					stringIndex = indexEdgeLengthEnd;
				} else {
					stringIndex = stringIndex + 1;
				}	
			}
		}
		TreeUtils.parseTipTimes(tree);
		TreeUtils.setNodeHeights(tree);
		//TreeUtils.parseTipLabels(tree); //if tips have labeled states
		return tree;
	} //END method
        
        public static TreeNode[] getSimMapTree(String treeFileName)
	{
	
		String treeString = treeFileReader.readTree(treeFileName);
		System.out.println("New Netwick tree:");
		System.out.println(treeString);
	
		int treeStringLength = treeString.length();
	
		//Find number of internal nodes
		int intNodeCount = treeString.replaceAll("[^)]", "").length(); //includes the root node
 
		//Find number of tip nodes
		int tipNodeCount = intNodeCount + 1;//treeString.replaceAll("[^,]", "").length() + 1;
		int nodeCount = intNodeCount + tipNodeCount;
	
		//Find number of edges
		int edgeCount = intNodeCount + tipNodeCount - 1;
		
		//Define characters
		char L_PAREN = '(';
		char COMMA = ',';
		char R_PAREN = ')';
		char COLON = ':';
                char L_CURLY = '{';
                char R_CURLY = '}';
		
		//Make tree object
		TreeNode[] tree = new TreeNode[nodeCount]; //change to ArrayList<TreeNode>
		for (int j = 0; j < nodeCount; ++j) {
			tree[j] = new TreeNode();
			tree[j].nodeNumber = j+1;
		}
		
		//Parse tree from newick format
		//Nodes 1:#Tips = tips
		//Nodes #Tips+1:#Nodes = internal nodes
		int currNode = tipNodeCount; //start bellow root
		int nextInternalNode = currNode + 1;
		int nextTipNode = 1;
		int stringIndex = 0;
                
		for (int i = 0; i < treeStringLength; ++i) {
			char currChar = treeString.charAt(stringIndex);
			
			if (currChar == L_PAREN) {
				//Create daughter node
				tree[nextInternalNode-1].setParent(currNode);
				tree[nextInternalNode-1].setNodeName("internal");
				tree[nextInternalNode-1].setNodeType("internal");
				//Move to daughter
				currNode = nextInternalNode; //4
				nextInternalNode = nextInternalNode + 1; //5
				tree[currNode-1].setChildren(0, nextInternalNode);
				if (currNode == (tipNodeCount + 1)) { // then its the root
					tree[currNode].setNodeType("root");
				}
				char nextChar = treeString.charAt(stringIndex+1);
				if (Character.isLetterOrDigit(nextChar)) {
					
                                        stringIndex = stringIndex + 1;
					
                                        //Get node name
					int indexLabelEnd = treeString.indexOf(COLON, (stringIndex));
					String tipLabel = treeString.substring(stringIndex, indexLabelEnd);
					stringIndex = indexLabelEnd;
					tree[nextTipNode-1].setParent(currNode);
					tree[nextTipNode-1].setNodeName(tipLabel);
					tree[nextTipNode-1].setNodeType("tip");
					tree[currNode-1].setChildren(0, nextTipNode);
					currNode = nextTipNode; //1
					nextTipNode = nextTipNode + 1; //2
					
                                        //Parse SimMap branch segment lengths:
                                        int indexNextCurly = treeString.indexOf(R_CURLY, stringIndex);
                                        String branchAnnotation = treeString.substring(stringIndex+2, indexNextCurly);
                                        String delimiter = ":";
                                        String[] branchSegments = branchAnnotation.split(delimiter);
                                        double totalBranchLength = 0;
                                        for (int seg = 0; seg < branchSegments.length; ++seg) {
                                            int indexComma = branchSegments[seg].indexOf(COMMA, 0);
                                            String segmentLength = branchSegments[seg].substring(indexComma+1);
                                            totalBranchLength += Double.parseDouble(segmentLength);
                                        }
                                        //if (totalBranchLength > 3652.0) {
                                            //System.out.println();
                                        //}
                                        tree[currNode-1].distanceToParent = totalBranchLength;
                                        tree[currNode-1].setNodeAnnotation(branchAnnotation);
                                        stringIndex = indexNextCurly + 1;
					
				} else {
					stringIndex = stringIndex + 1;
				}
			}
			
			if (currChar == COMMA) {
				//Go back to parent
				currNode = tree[currNode-1].getParent();
				char nextChar = treeString.charAt(stringIndex+1);
				if (Character.isLetterOrDigit(nextChar)) {
					stringIndex = stringIndex + 1;
					//Get node name
					int indexLabelEnd = treeString.indexOf(COLON, (stringIndex));
					String tipLabel = treeString.substring(stringIndex, indexLabelEnd);
					stringIndex = indexLabelEnd;
					tree[nextTipNode-1].setParent(currNode);
					tree[nextTipNode-1].setNodeName(tipLabel);
					tree[nextTipNode-1].setNodeType("tip");
					tree[currNode-1].setChildren(1, nextTipNode);
					currNode = nextTipNode;
					nextTipNode = nextTipNode + 1;
                                        
                                        //Parse SimMap branch segment lengths:
                                        int indexNextCurly = treeString.indexOf(R_CURLY, stringIndex);
                                        String branchAnnotation = treeString.substring(stringIndex+2, indexNextCurly);
                                        String delimiter = ":";
                                        String[] branchSegments = branchAnnotation.split(delimiter);
                                        double totalBranchLength = 0;
                                        for (int seg = 0; seg < branchSegments.length; ++seg) {
                                            int indexComma = branchSegments[seg].indexOf(COMMA, 0);
                                            String segmentLength = branchSegments[seg].substring(indexComma+1);
                                            totalBranchLength += Double.parseDouble(segmentLength);
                                        }
                                        //if (totalBranchLength > 3652.0) {
                                            //System.out.println();
                                        //}
                                        tree[currNode-1].distanceToParent = totalBranchLength;
                                        tree[currNode-1].setNodeAnnotation(branchAnnotation);
                                        stringIndex = indexNextCurly + 1;
                                        

				} else {
					tree[currNode-1].setChildren(1, nextInternalNode);
					stringIndex = stringIndex + 1;
				}
			}
			
			if (currChar == R_PAREN) {
				currNode = tree[currNode-1].getParent();
				char nextChar = treeString.charAt(stringIndex+1);
				if (nextChar == COLON) {
					
					stringIndex = stringIndex + 1; //now at COLON
                                        //Parse SimMap branch segment lengths:
                                        int indexNextCurly = treeString.indexOf(R_CURLY, stringIndex);
                                        String branchAnnotation = treeString.substring(stringIndex+2, indexNextCurly);
                                        String delimiter = ":";
                                        String[] branchSegments = branchAnnotation.split(delimiter);
                                        double totalBranchLength = 0;
                                        for (int seg = 0; seg < branchSegments.length; ++seg) {
                                            int indexComma = branchSegments[seg].indexOf(COMMA, 0);
                                            String segmentLength = branchSegments[seg].substring(indexComma+1);
                                            totalBranchLength += Double.parseDouble(segmentLength);
                                        }
                                        //if (totalBranchLength > 3652.0) {
                                            //System.out.println();
                                        //}
                                        tree[currNode-1].distanceToParent = totalBranchLength;
                                        tree[currNode-1].setNodeAnnotation(branchAnnotation);
                                        stringIndex = indexNextCurly + 1;

				} else {
					stringIndex = stringIndex + 1;
				}	
			}
		}
                
		TreeUtils.parseTipTimes(tree);
		TreeUtils.setNodeHeights(tree);
		TreeUtils.parseTipLabels(tree); //if tips have labeled states
		return tree;
	} //END method
	
} //END class