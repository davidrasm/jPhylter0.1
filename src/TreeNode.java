import java.util.ArrayList;

public class TreeNode implements Comparable<TreeNode>
{
	int nodeNumber;
	int parentNode;
	String nodeType;
	int[] childNodes = new int[2];
	double distanceToParent = 0.0;
	double nodeHeight;
	String nodeName;
	int nodeLabel;
        String nodeAnnotation;
        ArrayList<Double> statePriors = new ArrayList<Double>();
        
    //@Override
        public int compareTo(TreeNode e)
        {
            return Double.compare(nodeHeight, e.nodeHeight);
        }
	
	public void setNodeNumber(int number)
	{
		nodeNumber = number;
	}
	
	public void setChildren(int position, int child)
	{
		childNodes[position] = child;	
	}
	
	public void setNodeType(String type)
	{
		nodeType = type;
	}
	
	public void setParent(int nodeP)
	{
		parentNode = nodeP;
	}
	
	public void setNodeHeight(double height)
	{
		nodeHeight = height;
	}
	
	public void setNodeName(String name)
	{
		nodeName = name;
	}
	
	public void setDistance(String edgeDistance)
	{
		double doubleDistance = Double.parseDouble(edgeDistance);
                if (doubleDistance < 0.0) {
                    System.out.println("WARNING: Negative branch length found! Likelihood analysis will fail.");
                }
		distanceToParent = doubleDistance;
	}
	
	public int getParent()
	{
		return parentNode;
	}
	
	public void setNodeLabel(int label)
	{
		nodeLabel = label;
	}
        
        public void setNodeAnnotation(String annotation)
        {
            nodeAnnotation = annotation;
        }
	
}