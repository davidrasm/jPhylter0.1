import cern.colt.matrix.*;

/**
 * Set or get current state trajectory obtained from the particle filter
 * @author David
 */
public class XTrajectory {
    
    DoubleMatrix2D full;
    DoubleMatrix2D sub;
    int zLocStart;
    int zLocEnd;
    
    public void setXTrajFull(int xs, int ts) {
        
        //xStates = xs;
        //times = ts;
        //DoubleFactory2D factory2D;
	//factory2D = DoubleFactory2D.dense;
	//full = factory2D.make(xStates, times);
        
        //full.assign(100);
    }
    
    public void getSubTraj(int subStartLoc, int subEndLoc) {
        
        int subTimes = subEndLoc - subStartLoc + 1;
        int xStates = full.rows();
        
        DoubleFactory2D factory2D;
	factory2D = DoubleFactory2D.dense;
        sub = factory2D.make(xStates, subTimes);
        
        int subLoc = 0;
        for (int i = 0; i < subTimes; ++i) {
            for (int j = 0; j < xStates; ++j) {
                if (subStartLoc + i <= zLocStart) {
                    subLoc = zLocStart;
                    //System.out.println("SubLoc is: " + subLoc);
                }
                if (subStartLoc + i > zLocStart & subStartLoc + i <= zLocEnd) {
                    subLoc = (subStartLoc + i) - zLocStart;
                    //System.out.println("SubLoc is: " + subLoc);
                }
                if (subStartLoc + i > zLocEnd) {
                    subLoc = zLocEnd;
                    //System.out.println("SubLoc is: " + subLoc);
                }
                sub.setQuick(j, i, full.getQuick(j, subLoc));
            }
        }
    }
    
    public XTrajectory copy() 
    {
        XTrajectory copy = new XTrajectory();
        copy.full = this.full.copy();
        //copy.sub = this.sub.copy();
        copy.zLocStart = zLocStart;
        copy.zLocEnd = zLocEnd;
        return copy;
        
    }
    
    
}
