import java.util.ArrayList;
import java.util.Collections;
import cern.colt.list.DoubleArrayList;
/**
 *
 * @author David
 */
public class NewickWriter {
    
    String treeString;

    public void treeToString(TreeNode[] tree) {
    
        //Put internal nodes in ArrayList
        int nodes = tree.length;
        int internalNodes = ((nodes + 1) / 2) - 1;
        int tipNodes = nodes - internalNodes;
        ArrayList<TreeNode> treeArray = new ArrayList<TreeNode>();
        for (int n = 0; n < nodes; ++n) {
            treeArray.add(tree[n]);
        }
        
        //Sort internal nodes by node heights
        Collections.sort(treeArray);
        
        //Post order tree traversal
        for (int node = nodes - 1; node >= 0; --node) {
            
            if (treeArray.get(node).nodeType == "tip") {
                //skip
            } else {
                int child1 = treeArray.get(node).childNodes[0];
                int child2 = treeArray.get(node).childNodes[1];
                String child1Name = tree[child1 - 1].nodeName;
                String child2Name = tree[child2 - 1].nodeName;
                String child1Length = Double.toString(tree[child1 - 1].distanceToParent);
                String child2Length = Double.toString(tree[child2 - 1].distanceToParent);
                String newName = ('(' + child1Name + ':' + child1Length + ',' + child2Name + ":" + child2Length + ')');
                treeArray.get(node).nodeName = newName;
            }
        }  
        treeString = treeArray.get(0).nodeName + ';';  
    }
    
    public void probTreeToString(TreeNode[] tree, LineageStateSamples samples) 
    {
        //This method only works if we're starting at the root, otherwise LineageEntropies will not have absoluteTimes before the startTime
        int totalClasses = 20;
        samples.discretizeSingleStatePosteriorProbs(totalClasses);
        
        int nodes = tree.length;
        int internalNodes = ((nodes + 1) / 2) - 1;
        int tipNodes = nodes - internalNodes;
        ArrayList<TreeNode> treeArray = new ArrayList<TreeNode>();
        for (int n = 0; n < nodes; ++n) {
            treeArray.add(tree[n]);
        }
        
        //Sort internal nodes by node heights
        Collections.sort(treeArray);
        
        //Post order tree traversal
        for (int node = nodes - 1; node >= 0; --node) {
            //System.out.println("Node =" + node);
            if (treeArray.get(node).nodeType == "tip") {
                //do nothing
            } else {
            
                int child1 = treeArray.get(node).childNodes[0];
                int child2 = treeArray.get(node).childNodes[1];
                String child1Name = tree[child1 - 1].nodeName;
                String child2Name = tree[child2 - 1].nodeName;

                //Get lineage posterior probs for child1
                DoubleArrayList lineageTimesChild1 = new DoubleArrayList();
                if (samples.absoluteTimes.size() >= (child1 - 1)) {
                    lineageTimesChild1 = samples.absoluteTimes.get(child1 - 1);
                } else {
                    System.out.println();
                }
                int totalTimesChild1 = lineageTimesChild1.size();
                lineageTimesChild1.sort();
                lineageTimesChild1.reverse();
                Collections.reverse(samples.discreteProbs.get(child1 - 1));
                int t = 0;
                double currTime = lineageTimesChild1.get(0); //moving backwards through time
                int currIndexLoc = samples.absoluteTimes.get(child1 - 1).indexOf(currTime);
                int currClass = samples.discreteProbs.get(child1 - 1).get(currIndexLoc);
                double segmentLength = 0.0;
                String child1String = "{";
                while (t < totalTimesChild1 - 1) {
                    double nextTime = lineageTimesChild1.get(t+1);
                    int nextIndexLoc = samples.absoluteTimes.get(child1 - 1).indexOf(nextTime);
                    int nextClass = samples.discreteProbs.get(child1 - 1).get(nextIndexLoc);
                    if ((t+1) != totalTimesChild1 - 1) {
                        if (nextClass == currClass) {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            currTime = nextTime;
                        } else {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            String segLengthString = Double.toString(segmentLength);
                            String segClassString = Integer.toString(currClass);
                            String stringForSeg = segClassString + ',' + segLengthString + ':'; 
                            child1String = child1String + stringForSeg;
                            segmentLength = 0.0;
                            currClass = nextClass;
                            currTime = nextTime;
                        }
                    } else {
                        double addTime = currTime - nextTime;
                        segmentLength += addTime;
                        String segLengthString = Double.toString(segmentLength);
                        String segClassString = Integer.toString(currClass);
                        String stringForSeg = segClassString + ',' + segLengthString + '}'; 
                        child1String = child1String + stringForSeg;    
                    }
                    t += 1; //increment t
                }

                //Get lineage entropies for child2
                DoubleArrayList lineageTimesChild2 = samples.absoluteTimes.get(child2 - 1);
                int totalTimesChild2 = lineageTimesChild2.size();
                lineageTimesChild2.sort();
                lineageTimesChild2.reverse();
                Collections.reverse(samples.discreteProbs.get(child2 - 1));
                t = 0;
                currTime = lineageTimesChild2.get(0); //moving backwards through time
                currIndexLoc = samples.absoluteTimes.get(child2 - 1).indexOf(currTime);
                currClass = samples.discreteProbs.get(child2 - 1).get(currIndexLoc);
                segmentLength = 0.0;
                String child2String = "{";
                while (t < totalTimesChild2 - 1) {
                    double nextTime = lineageTimesChild2.get(t+1);
                    int nextIndexLoc = samples.absoluteTimes.get(child2 - 1).indexOf(nextTime);
                    int nextClass = samples.discreteProbs.get(child2 - 1).get(nextIndexLoc);
                    if ((t+1) != totalTimesChild2 - 1) {
                        if (nextClass == currClass) {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            currTime = nextTime;
                        } else {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            String segLengthString = Double.toString(segmentLength);
                            String segClassString = Integer.toString(currClass);
                            String stringForSeg = segClassString + ',' + segLengthString + ':'; 
                            child2String = child2String + stringForSeg;
                            segmentLength = 0.0;
                            currClass = nextClass;
                            currTime = nextTime;
                        }
                    } else {
                        double addTime = currTime - nextTime;
                        segmentLength += addTime;
                        String segLengthString = Double.toString(segmentLength);
                        String segClassString = Integer.toString(currClass);
                        String stringForSeg = segClassString + ',' + segLengthString + '}'; 
                        child2String = child2String + stringForSeg;    
                    }
                    t += 1; //increment t
                }

                //Get new name for this node
                String newName = ('(' + child1Name + ':' + child1String + ',' + child2Name + ":" + child2String + ')');
                treeArray.get(node).nodeName = newName;
            } //end IF
           
        } //end FOR  
        treeString = treeArray.get(0).nodeName + ';';
        //System.out.println();
        
    }
    
    public void simMapTreeToString(ParticleLineageStateArray S, TreeNode[] tree, ZVectors dataZ) 
    {
        //This method only works if we're starting at the root, otherwise LineageEntropies will not have absoluteTimes before the startTime
        //int totalClasses = 20;
        //samples.discretizeSingleStatePosteriorProbs(totalClasses);
        
        int nodes = tree.length;
        int internalNodes = ((nodes + 1) / 2) - 1;
        int tipNodes = nodes - internalNodes;
        ArrayList<TreeNode> treeArray = new ArrayList<TreeNode>();
        for (int n = 0; n < nodes; ++n) {
            treeArray.add(tree[n]);
        }
        
        //Sort internal nodes by node heights
        Collections.sort(treeArray);
        
        //Post order tree traversal
        for (int node = nodes - 1; node >= 0; --node) {
            //System.out.println("Node =" + node);
            if (treeArray.get(node).nodeType == "tip") {
                //do nothing
            } else {
            
                int child1 = treeArray.get(node).childNodes[0];
                int child2 = treeArray.get(node).childNodes[1];
                String child1Name = tree[child1 - 1].nodeName;
                String child2Name = tree[child2 - 1].nodeName;

                //Get lineage states for child1
                double child1StartTime = treeArray.get(node).nodeHeight;
                int zIndexStart = dataZ.absoluteTimes.indexOf(child1StartTime);
                double child1EndTime = tree[child1 - 1].nodeHeight;
                int zIndexEnd = dataZ.absoluteTimes.indexOf(child1EndTime);
                DoubleArrayList lineageTimesChild1 = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zIndexStart, zIndexEnd);
                int totalTimesChild1 = lineageTimesChild1.size();
                lineageTimesChild1.sort();
                lineageTimesChild1.reverse();
                int t = 0;
                double currTime = lineageTimesChild1.get(0); //moving backwards through time
                int mapTimeIndex = S.absoluteTimes.indexOf(currTime);
                int currState = 0;
                if (mapTimeIndex >= 0) {
                    int currMapIndex = S.arrayMap.get(mapTimeIndex).indexOf(child1 - 1);;
                    currState = S.sampleArray.get(mapTimeIndex).get(currMapIndex);
                } else {
                    currState = 3;
                }
                double segmentLength = 0.0;
                String child1String = "{";
                while (t < totalTimesChild1 - 1) {
                    double nextTime = lineageTimesChild1.get(t+1);
                    mapTimeIndex--;
                    int nextState = 0;
                    if (mapTimeIndex >= 0) {
                        int nextMapIndex = S.arrayMap.get(mapTimeIndex).indexOf(child1 - 1);
                        nextState = S.sampleArray.get(mapTimeIndex).get(nextMapIndex);
                    } else {
                        nextState = 3;
                    }
                    if ((t+1) != totalTimesChild1 - 1) { //if not at lineage start time
                        if (nextState == currState) {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            currTime = nextTime;
                        } else {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            String segLengthString = Double.toString(segmentLength);
                            String segClassString = Integer.toString(currState + 1);
                            String stringForSeg = segClassString + ',' + segLengthString + ':'; 
                            child1String = child1String + stringForSeg;
                            segmentLength = 0.0;
                            currState = nextState;
                            currTime = nextTime;
                        }
                    } else {
                        double addTime = currTime - nextTime;
                        segmentLength += addTime;
                        String segLengthString = Double.toString(segmentLength);
                        String segClassString = Integer.toString(currState + 1);
                        String stringForSeg = segClassString + ',' + segLengthString + '}'; 
                        child1String = child1String + stringForSeg;    
                    }
                    t += 1; //increment t
                }

                //Get lineage states for child2
                double child2StartTime = treeArray.get(node).nodeHeight;
                zIndexStart = dataZ.absoluteTimes.indexOf(child2StartTime);
                double child2EndTime = tree[child2 - 1].nodeHeight;
                zIndexEnd = dataZ.absoluteTimes.indexOf(child2EndTime);
                DoubleArrayList lineageTimesChild2 = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zIndexStart, zIndexEnd);
                int totalTimesChild2 = lineageTimesChild2.size();
                lineageTimesChild2.sort();
                lineageTimesChild2.reverse();
                t = 0;
                currTime = lineageTimesChild2.get(0); //moving backwards through time
                mapTimeIndex = S.absoluteTimes.indexOf(currTime);
                currState = 0;
                if (mapTimeIndex >= 0) {
                    int currMapIndex = S.arrayMap.get(mapTimeIndex).indexOf(child2 - 1);;
                    currState = S.sampleArray.get(mapTimeIndex).get(currMapIndex);
                } else {
                    currState = 3;
                }
                segmentLength = 0.0;
                String child2String = "{";
                while (t < totalTimesChild2 - 1) {
                    double nextTime = lineageTimesChild2.get(t+1);
                    mapTimeIndex--;
                    int nextState = 0;
                    if (mapTimeIndex >= 0) {
                        int nextMapIndex = S.arrayMap.get(mapTimeIndex).indexOf(child2 - 1);
                        nextState = S.sampleArray.get(mapTimeIndex).get(nextMapIndex);
                    } else {
                        nextState = 3;
                    }
                    if ((t+1) != totalTimesChild2 - 1) { //if not at lineage start time
                        if (nextState == currState) {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            currTime = nextTime;
                        } else {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            String segLengthString = Double.toString(segmentLength);
                            String segClassString = Integer.toString(currState + 1);
                            String stringForSeg = segClassString + ',' + segLengthString + ':'; 
                            child2String = child2String + stringForSeg;
                            segmentLength = 0.0;
                            currState = nextState;
                            currTime = nextTime;
                        }
                    } else {
                        double addTime = currTime - nextTime;
                        segmentLength += addTime;
                        String segLengthString = Double.toString(segmentLength);
                        String segClassString = Integer.toString(currState + 1);
                        String stringForSeg = segClassString + ',' + segLengthString + '}'; 
                        child2String = child2String + stringForSeg;    
                    }
                    t += 1; //increment t
                }
                


                //Get new name for this node
                String newName = ('(' + child1Name + ':' + child1String + ',' + child2Name + ":" + child2String + ')');
                treeArray.get(node).nodeName = newName;
            } //end IF
           
        } //end FOR  
        treeString = treeArray.get(0).nodeName + ';';
        //System.out.println();
        
    }
    
    
    public void entropyTreeToString(TreeNode[] tree, LineageStateSamples samples) 
    {
        //This method only works if we're starting at the root, otherwise LineageEntropies will not have absoluteTimes before the startTime
        int totalClasses = 20;
        samples.discretizeEntropies(totalClasses);
        
        int nodes = tree.length;
        int internalNodes = ((nodes + 1) / 2) - 1;
        int tipNodes = nodes - internalNodes;
        ArrayList<TreeNode> treeArray = new ArrayList<TreeNode>();
        for (int n = 0; n < nodes; ++n) {
            treeArray.add(tree[n]);
        }
        
        //Sort internal nodes by node heights
        Collections.sort(treeArray);
        
        //Post order tree traversal
        for (int node = nodes - 1; node >= 0; --node) {
            //System.out.println("Node =" + node);
            if (treeArray.get(node).nodeType == "tip") {
                //do nothing
            } else {
            
                int child1 = treeArray.get(node).childNodes[0];
                int child2 = treeArray.get(node).childNodes[1];
                String child1Name = tree[child1 - 1].nodeName;
                String child2Name = tree[child2 - 1].nodeName;

                //Get lineage entropies for child1
                DoubleArrayList lineageTimesChild1 = samples.absoluteTimes.get(child1 - 1);
                int totalTimesChild1 = lineageTimesChild1.size();
                lineageTimesChild1.sort();
                lineageTimesChild1.reverse();
                Collections.reverse(samples.discreteEntropies.get(child1 - 1));
                int t = 0;
                double currTime = lineageTimesChild1.get(0); //moving backwards through time
                int currIndexLoc = samples.absoluteTimes.get(child1 - 1).indexOf(currTime);
                int currClass = samples.discreteEntropies.get(child1 - 1).get(currIndexLoc);
                double segmentLength = 0.0;
                String child1String = "{";
                while (t < totalTimesChild1 - 1) {
                    double nextTime = lineageTimesChild1.get(t+1);
                    int nextIndexLoc = samples.absoluteTimes.get(child1 - 1).indexOf(nextTime);
                    int nextClass = samples.discreteEntropies.get(child1 - 1).get(nextIndexLoc);
                    if ((t+1) != totalTimesChild1 - 1) {
                        if (nextClass == currClass) {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            currTime = nextTime;
                        } else {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            String segLengthString = Double.toString(segmentLength);
                            String segClassString = Integer.toString(currClass);
                            String stringForSeg = segClassString + ',' + segLengthString + ':'; 
                            child1String = child1String + stringForSeg;
                            segmentLength = 0.0;
                            currClass = nextClass;
                            currTime = nextTime;
                        }
                    } else {
                        double addTime = currTime - nextTime;
                        segmentLength += addTime;
                        String segLengthString = Double.toString(segmentLength);
                        String segClassString = Integer.toString(currClass);
                        String stringForSeg = segClassString + ',' + segLengthString + '}'; 
                        child1String = child1String + stringForSeg;    
                    }
                    t += 1; //increment t
                }

                //Get lineage entropies for child2
                DoubleArrayList lineageTimesChild2 = samples.absoluteTimes.get(child2 - 1);
                int totalTimesChild2 = lineageTimesChild2.size();
                lineageTimesChild2.sort();
                lineageTimesChild2.reverse();
                Collections.reverse(samples.discreteEntropies.get(child2 - 1));
                t = 0;
                currTime = lineageTimesChild2.get(0); //moving backwards through time
                currIndexLoc = samples.absoluteTimes.get(child2 - 1).indexOf(currTime);
                currClass = samples.discreteEntropies.get(child2 - 1).get(currIndexLoc);
                segmentLength = 0.0;
                String child2String = "{";
                while (t < totalTimesChild2 - 1) {
                    double nextTime = lineageTimesChild2.get(t+1);
                    int nextIndexLoc = samples.absoluteTimes.get(child2 - 1).indexOf(nextTime);
                    int nextClass = samples.discreteEntropies.get(child2 - 1).get(nextIndexLoc);
                    if ((t+1) != totalTimesChild2 - 1) {
                        if (nextClass == currClass) {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            currTime = nextTime;
                        } else {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            String segLengthString = Double.toString(segmentLength);
                            String segClassString = Integer.toString(currClass);
                            String stringForSeg = segClassString + ',' + segLengthString + ':'; 
                            child2String = child2String + stringForSeg;
                            segmentLength = 0.0;
                            currClass = nextClass;
                            currTime = nextTime;
                        }
                    } else {
                        double addTime = currTime - nextTime;
                        segmentLength += addTime;
                        String segLengthString = Double.toString(segmentLength);
                        String segClassString = Integer.toString(currClass);
                        String stringForSeg = segClassString + ',' + segLengthString + '}'; 
                        child2String = child2String + stringForSeg;    
                    }
                    t += 1; //increment t
                }

                //Get new name for this node
                String newName = ('(' + child1Name + ':' + child1String + ',' + child2Name + ":" + child2String + ')');
                treeArray.get(node).nodeName = newName;
            } //end IF
           
        } //end FOR  
        treeString = treeArray.get(0).nodeName + ';';
        //System.out.println();
        
    }
    
    public void oldEntropyTreeToString(TreeNode[] tree, LineageEntropies samples) 
    {
        //This method only works if we're starting at the root, otherwise LineageEntropies will not have absoluteTimes before the startTime
        int totalClasses = 20;
        samples.discretizeEntropies(totalClasses);
        
        int nodes = tree.length;
        int internalNodes = ((nodes + 1) / 2) - 1;
        int tipNodes = nodes - internalNodes;
        ArrayList<TreeNode> treeArray = new ArrayList<TreeNode>();
        for (int n = 0; n < nodes; ++n) {
            treeArray.add(tree[n]);
        }
        
        //Sort internal nodes by node heights
        Collections.sort(treeArray);
        
        //Post order tree traversal
        for (int node = nodes - 1; node >= 0; --node) {
            //System.out.println("Node =" + node);
            if (treeArray.get(node).nodeType == "tip") {
                //do nothing
            } else {
            
                int child1 = treeArray.get(node).childNodes[0];
                int child2 = treeArray.get(node).childNodes[1];
                String child1Name = tree[child1 - 1].nodeName;
                String child2Name = tree[child2 - 1].nodeName;

                //Get lineage entropies for child1
                DoubleArrayList lineageTimesChild1 = samples.absoluteTimes.get(child1 - 1);
                int totalTimesChild1 = lineageTimesChild1.size();
                lineageTimesChild1.sort();
                lineageTimesChild1.reverse();
                //Collections.reverse(samples.discreteEntropies.get(child1 - 1));
                int t = 0;
                double currTime = lineageTimesChild1.get(0); //moving backwards through time
                int currIndexLoc = samples.absoluteTimes.get(child1 - 1).indexOf(currTime);
                int currClass = samples.discreteEntropies.get(child1 - 1).get(currIndexLoc);
                double segmentLength = 0.0;
                String child1String = "{";
                while (t < totalTimesChild1 - 1) {
                    double nextTime = lineageTimesChild1.get(t+1);
                    int nextIndexLoc = samples.absoluteTimes.get(child1 - 1).indexOf(nextTime);
                    int nextClass = samples.discreteEntropies.get(child1 - 1).get(nextIndexLoc);
                    if ((t+1) != totalTimesChild1 - 1) {
                        if (nextClass == currClass) {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            currTime = nextTime;
                        } else {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            String segLengthString = Double.toString(segmentLength);
                            String segClassString = Integer.toString(currClass);
                            String stringForSeg = segClassString + ',' + segLengthString + ':'; 
                            child1String = child1String + stringForSeg;
                            segmentLength = 0.0;
                            currClass = nextClass;
                            currTime = nextTime;
                        }
                    } else {
                        double addTime = currTime - nextTime;
                        segmentLength += addTime;
                        String segLengthString = Double.toString(segmentLength);
                        String segClassString = Integer.toString(currClass);
                        String stringForSeg = segClassString + ',' + segLengthString + '}'; 
                        child1String = child1String + stringForSeg;    
                    }
                    t += 1; //increment t
                }

                //Get lineage entropies for child2
                DoubleArrayList lineageTimesChild2 = samples.absoluteTimes.get(child2 - 1);
                int totalTimesChild2 = lineageTimesChild2.size();
                lineageTimesChild2.sort();
                lineageTimesChild2.reverse();
                //Collections.reverse(samples.discreteEntropies.get(child2 - 1));
                t = 0;
                currTime = lineageTimesChild2.get(0); //moving backwards through time
                currIndexLoc = samples.absoluteTimes.get(child2 - 1).indexOf(currTime);
                currClass = samples.discreteEntropies.get(child2 - 1).get(currIndexLoc);
                segmentLength = 0.0;
                String child2String = "{";
                while (t < totalTimesChild2 - 1) {
                    double nextTime = lineageTimesChild2.get(t+1);
                    int nextIndexLoc = samples.absoluteTimes.get(child2 - 1).indexOf(nextTime);
                    int nextClass = samples.discreteEntropies.get(child2 - 1).get(nextIndexLoc);
                    if ((t+1) != totalTimesChild2 - 1) {
                        if (nextClass == currClass) {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            currTime = nextTime;
                        } else {
                            double addTime = currTime - nextTime;
                            segmentLength += addTime;
                            String segLengthString = Double.toString(segmentLength);
                            String segClassString = Integer.toString(currClass);
                            String stringForSeg = segClassString + ',' + segLengthString + ':'; 
                            child2String = child2String + stringForSeg;
                            segmentLength = 0.0;
                            currClass = nextClass;
                            currTime = nextTime;
                        }
                    } else {
                        double addTime = currTime - nextTime;
                        segmentLength += addTime;
                        String segLengthString = Double.toString(segmentLength);
                        String segClassString = Integer.toString(currClass);
                        String stringForSeg = segClassString + ',' + segLengthString + '}'; 
                        child2String = child2String + stringForSeg;    
                    }
                    t += 1; //increment t
                }

                //Get new name for this node
                String newName = ('(' + child1Name + ':' + child1String + ',' + child2Name + ":" + child2String + ')');
                treeArray.get(node).nodeName = newName;
            } //end IF
           
        } //end FOR  
        treeString = treeArray.get(0).nodeName + ';';
        //System.out.println();
        
    }
    
    
}
