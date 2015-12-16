import java.util.ArrayList;
/**
 * 
 * @author David
 */
public class Event implements Comparable<Event> {
    
    double absoluteTime;
    int eventID;
    ArrayList<Integer> eventTypes = new ArrayList<Integer>();
    ArrayList<Integer> nodes = new ArrayList<Integer>();
    
    //@Override
    public int compareTo(Event e)
    {
        return Double.compare(absoluteTime, e.absoluteTime);
    }
}
