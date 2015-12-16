import cern.colt.list.*;

/**
 *
 * @author dar24
 */
public class ArrayUtils {
    
    public static int getIndexOfClosestDouble(DoubleArrayList array, double target) 
    {
        double cDistance = 0;;
        double distance = Math.abs(array.getQuick(0) - target);
        int index = 0;
        for (int c = 1; c < array.size(); c++) {
            cDistance = Math.abs(array.getQuick(c) - target);
            if (cDistance < distance) {
                index = c;
                distance = cDistance;
            }
                    
        }
        return index;
    }
    
    
}
