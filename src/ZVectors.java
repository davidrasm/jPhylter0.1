import cern.colt.list.*;
import java.util.ArrayList;
import java.util.Collections;

/**
 * ZVectors is the main data structure in PHYLter used to store the "events" in the tree (e.g. sampling, coalescence, ect.) and observation data (observation times)
 * Multiple events of different types can occur at the same time point
 * @author David
 */

public class ZVectors
{
        ArrayList<ArrayList<Integer>> omegaEvents;
	DoubleArrayList absoluteTimes;
	DoubleArrayList omegaLineages;
	ArrayList<ArrayList<Integer>> nodePointers;
        DoubleArrayList inactiveLineages;
        double startTime;
        double endTime;
        int startZSubLoc;
        int endZSubLoc;
        int thirdChildZSubLoc;
	
	public void setZVectors(TreeNode[] tree, double startTimeIn, double endTimeIn, double dt)
	{ 
                startTime = startTimeIn; //this can change based on the position of the root
                endTime = endTimeIn; //this probably should never change
            
                omegaEvents = new ArrayList<ArrayList<Integer>>(); //type of event
                absoluteTimes = new DoubleArrayList(); //absolute times of events
                omegaLineages = new DoubleArrayList(); //number of lineages in the genealogy at this time
                nodePointers = new ArrayList<ArrayList<Integer>>(); //pointer to tree node for events in the genealogy
            
                DoubleArrayList sampleTimes = TreeUtils.getTipTimes(tree);
		ArrayList<Integer> sampleLocs = TreeUtils.getTipLocs(tree);
		DoubleArrayList coalTimes = TreeUtils.getInternalNodeTimes(tree);
		ArrayList<Integer> coalLocs = TreeUtils.getInternalNodeLocs(tree);
                
                ArrayList<Event> events = new ArrayList<Event>();
                DoubleArrayList eventTimes = new DoubleArrayList();
                
                //Add sampling events to events array
                for (int sample = 0; sample < sampleTimes.size(); ++sample) {
                    double thisTime = sampleTimes.get(sample);
                    int eventLoc = eventTimes.indexOf(thisTime);
                    if (eventLoc == -1) {
                        //new event time
                        Event thisEvent = new Event();
                        thisEvent.eventID = events.size();
                        thisEvent.absoluteTime = thisTime;
                        thisEvent.eventTypes.add(2);
                        thisEvent.nodes.add(sampleLocs.get(sample)); 
                        events.add(thisEvent);
                        eventTimes.add(thisTime);
                    } else {
                        //event time already exists
                        events.get(eventLoc).eventTypes.add(2);
                        events.get(eventLoc).nodes.add(sampleLocs.get(sample));
                    }      
                }
                
                //Add coalescence events to events array
                for (int coal = 0; coal < coalTimes.size(); ++coal) {
                    double thisTime = coalTimes.get(coal);
                    int eventLoc = eventTimes.indexOf(thisTime);
                    if (eventLoc == -1) {
                        Event thisEvent = new Event();
                        thisEvent.eventID = events.size();
                        thisEvent.absoluteTime = thisTime;
                        thisEvent.eventTypes.add(1);
                        thisEvent.nodes.add(coalLocs.get(coal));
                        events.add(thisEvent);
                        eventTimes.add(thisTime);
                    } else {
                        events.get(eventLoc).eventTypes.add(1);
                        events.get(eventLoc).nodes.add(coalLocs.get(coal)); 
                    }
                }
                
                
                //Get dt observation times;
		coalTimes.sort();
                sampleTimes.sort();
		double treeStartTime = coalTimes.get(0);
		double treeEndTime = sampleTimes.get(sampleTimes.size()-1);
		if (startTime == Double.NEGATIVE_INFINITY) {
			startTime = treeStartTime;
		}
		if (endTime == Double.POSITIVE_INFINITY) {
			endTime = treeEndTime;
		}
		
		DoubleArrayList obsvTimes = new DoubleArrayList(); //also includes dt times
		if (dt == Double.POSITIVE_INFINITY) { //startTime is root time
                    
                    //if (startTime != treeStartTime) {
                        obsvTimes.add(startTime);
                    //}
                    //if (endTime != treeEndTime) {    
			obsvTimes.add(endTime);
                    //}
		} else {
                    obsvTimes.add(startTime);
                    double newTime = startTime;
                    while (newTime < endTime) {
			newTime = newTime + dt;
                        if (newTime <= endTime) {
                            obsvTimes.add(newTime);
                        }
                    }
		}
                
                //Add observation events to events array
                for (int times = 0; times < obsvTimes.size(); ++times) {
                    double thisTime = obsvTimes.get(times);
                    int eventLoc = eventTimes.indexOf(thisTime);
                    if (eventLoc == -1) {
                        Event thisEvent = new Event();
                        thisEvent.eventID = events.size();
                        thisEvent.absoluteTime = thisTime;
                        thisEvent.eventTypes.add(0);
                        thisEvent.nodes.add(-1);
                        events.add(thisEvent);
                        eventTimes.add(thisTime);
                    } else {
                        events.get(eventLoc).eventTypes.add(0);
                        events.get(eventLoc).nodes.add(-1); 
                    }
                }
                
                Collections.sort(events); //sort based on absoluteTimes
		
		//Set Z vectors
		int numLineages = 0;
                for (int entry = events.size() - 1; entry >= 0; --entry) {
                    
                     for (int subEntry = 0; subEntry < events.get(entry).eventTypes.size(); ++subEntry) {   
                        
                        int type = (events.get(entry).eventTypes.get(subEntry));
                        if (type == 2) {
                            ++numLineages; //sample event
                        }
                        if (type == 1) {
                            --numLineages; //coal event
                        }
                     }
                     
                     omegaEvents.add(events.get(entry).eventTypes);
                     absoluteTimes.add(events.get(entry).absoluteTime);
                     omegaLineages.add(numLineages);
                     nodePointers.add(events.get(entry).nodes);
                     
                }     
		
		absoluteTimes.reverse();
                omegaLineages.reverse();
		omegaLineages.remove(0); //first element corresponds to interval before t = 0;
                omegaLineages.add(0); //add zero for time interval after t = end, so ZVectors have same length;
		Collections.reverse(omegaEvents);
                Collections.reverse(nodePointers);
	} //END method
        
        public void setZVectorsWithTimeSeries(TreeNode[] tree, double startTimeIn, double endTimeIn, double dt, TimeSeries timeSeries)
	{ 
                startTime = startTimeIn; //this can change based on the position of the root
                endTime = endTimeIn; //this probably should never change
            
                omegaEvents = new ArrayList<ArrayList<Integer>>(); //type of event
                absoluteTimes = new DoubleArrayList(); //absolute times of events
                omegaLineages = new DoubleArrayList(); //number of lineages in the genealogy at this time
                nodePointers = new ArrayList<ArrayList<Integer>>(); //pointer to tree node for events in the genealogy
            
                DoubleArrayList sampleTimes = TreeUtils.getTipTimes(tree);
		ArrayList<Integer> sampleLocs = TreeUtils.getTipLocs(tree);
		DoubleArrayList coalTimes = TreeUtils.getInternalNodeTimes(tree);
		ArrayList<Integer> coalLocs = TreeUtils.getInternalNodeLocs(tree);
                DoubleArrayList timeSeriesTimes = timeSeries.times;
                
                ArrayList<Event> events = new ArrayList<Event>();
                DoubleArrayList eventTimes = new DoubleArrayList();
                
                //Add sampling events to events array
                for (int sample = 0; sample < sampleTimes.size(); ++sample) {
                    double thisTime = sampleTimes.get(sample);
                    int eventLoc = eventTimes.indexOf(thisTime);
                    if (eventLoc == -1) {
                        //new event time
                        Event thisEvent = new Event();
                        thisEvent.eventID = events.size();
                        thisEvent.absoluteTime = thisTime;
                        thisEvent.eventTypes.add(2);
                        thisEvent.nodes.add(sampleLocs.get(sample)); 
                        events.add(thisEvent);
                        eventTimes.add(thisTime);
                    } else {
                        //event time already exists
                        events.get(eventLoc).eventTypes.add(2);
                        events.get(eventLoc).nodes.add(sampleLocs.get(sample));
                    }      
                }
                
                //Add coalescence events to events array
                for (int coal = 0; coal < coalTimes.size(); ++coal) {
                    double thisTime = coalTimes.get(coal);
                    int eventLoc = eventTimes.indexOf(thisTime);
                    if (eventLoc == -1) {
                        Event thisEvent = new Event();
                        thisEvent.eventID = events.size();
                        thisEvent.absoluteTime = thisTime;
                        thisEvent.eventTypes.add(1);
                        thisEvent.nodes.add(coalLocs.get(coal));
                        events.add(thisEvent);
                        eventTimes.add(thisTime);
                    } else {
                        events.get(eventLoc).eventTypes.add(1);
                        events.get(eventLoc).nodes.add(coalLocs.get(coal)); 
                    }
                }
                
                
                //Get dt observation times;
		coalTimes.sort();
                sampleTimes.sort();
		double treeStartTime = coalTimes.get(0);
		double treeEndTime = sampleTimes.get(sampleTimes.size()-1);
		if (startTime == Double.NEGATIVE_INFINITY) {
			startTime = treeStartTime;
		}
		if (endTime == Double.POSITIVE_INFINITY) {
			endTime = treeEndTime;
		}
		
		DoubleArrayList obsvTimes = new DoubleArrayList(); //also includes dt times
		if (dt == Double.POSITIVE_INFINITY) { //startTime is root time
                    
                    //if (startTime != treeStartTime) {
                        obsvTimes.add(startTime);
                    //}
                    //if (endTime != treeEndTime) {    
			obsvTimes.add(endTime);
                    //}
		} else {
                    obsvTimes.add(startTime);
                    double newTime = startTime;
                    while (newTime < endTime) {
			newTime = newTime + dt;
                        if (newTime <= endTime) {
                            obsvTimes.add(newTime);
                        }
                    }
		}
                
                //Add observation events to events array
                for (int times = 0; times < obsvTimes.size(); ++times) {
                    double thisTime = obsvTimes.get(times);
                    int eventLoc = eventTimes.indexOf(thisTime);
                    if (eventLoc == -1) {
                        Event thisEvent = new Event();
                        thisEvent.eventID = events.size();
                        thisEvent.absoluteTime = thisTime;
                        thisEvent.eventTypes.add(0);
                        thisEvent.nodes.add(-1);
                        events.add(thisEvent);
                        eventTimes.add(thisTime);
                    } else {
                        events.get(eventLoc).eventTypes.add(0);
                        events.get(eventLoc).nodes.add(-1); 
                    }
                }
                
                //Add timeSeries events to events array
                for (int times = 0; times < timeSeriesTimes.size(); ++times) {
                    double thisTime = timeSeriesTimes.get(times);
                    int eventLoc = eventTimes.indexOf(thisTime);
                    if (eventLoc == -1) {
                        Event thisEvent = new Event();
                        thisEvent.eventID = events.size();
                        thisEvent.absoluteTime = thisTime;
                        thisEvent.eventTypes.add(3);
                        thisEvent.nodes.add(-1);
                        events.add(thisEvent);
                        eventTimes.add(thisTime);
                    } else {
                        events.get(eventLoc).eventTypes.add(3);
                        events.get(eventLoc).nodes.add(-1); 
                    }
                }
                
                Collections.sort(events); //sort based on absoluteTimes
		
		//Set Z vectors
		int numLineages = 0;
                for (int entry = events.size() - 1; entry >= 0; --entry) {
                    
                     for (int subEntry = 0; subEntry < events.get(entry).eventTypes.size(); ++subEntry) {   
                        
                        int type = (events.get(entry).eventTypes.get(subEntry));
                        if (type == 2) {
                            ++numLineages; //sample event
                        }
                        if (type == 1) {
                            --numLineages; //coal event
                        }
                     }
                     
                     omegaEvents.add(events.get(entry).eventTypes);
                     absoluteTimes.add(events.get(entry).absoluteTime);
                     omegaLineages.add(numLineages);
                     nodePointers.add(events.get(entry).nodes);
                     
                }     
		
		absoluteTimes.reverse();
                omegaLineages.reverse();
		omegaLineages.remove(0); //first element corresponds to interval before t = 0;
                omegaLineages.add(0); //add zero for time interval after t = end, so ZVectors have same length;
		Collections.reverse(omegaEvents);
                Collections.reverse(nodePointers);
	} //END method
        
        public void viewZVectors()
        {
            //Just for debugging
            for (int aT = 0; aT < omegaEvents.size(); ++aT) {
                for (int subEntry = 0; subEntry < omegaEvents.get(aT).size(); ++subEntry) {
                    System.out.println(absoluteTimes.get(aT) + "  " + omegaEvents.get(aT).get(subEntry) + "  " + nodePointers.get(aT).get(subEntry));
                }
                System.out.println("Lineages:" + " " + omegaLineages.get(aT));
                //System.out.println("Inactive lineages:" + " " + inactiveLineages.get(aT));
                System.out.println();
            }
            
        }//END method
	
}//END class
