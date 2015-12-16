import cern.colt.list.*;
import java.io.*;

/**
 *
 * @author David
 */
public class LineagesThroughTime {
    
    public void printLTT(ZVectors dataZ) throws IOException
    {
        SampleWriter sampleWriter = new SampleWriter();
        FileWriter linTTFileWriter = new FileWriter("PHYLterLTTData");
        PrintWriter linTTPrintWriter = new PrintWriter(linTTFileWriter);
        sampleWriter.samplesToString(dataZ.absoluteTimes);
        linTTPrintWriter.print(sampleWriter.sampleString); linTTPrintWriter.println();
        sampleWriter.samplesToString(dataZ.omegaLineages);
        linTTPrintWriter.print(sampleWriter.sampleString); linTTPrintWriter.println();
        
        linTTPrintWriter.flush();
        linTTPrintWriter.close();
        linTTFileWriter.close();
    }
    
    
}
