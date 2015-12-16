import java.io.*;

public class FastaReader
{

    public static SeqAlign getAlignFromFasta(String fastaFileName)
    {
        String line = null;
	StringBuilder seqLines = new StringBuilder();
	StringBuilder headerLines = new StringBuilder();
	String delim =  ">";
	String seqReadStart = "FALSE";
	try {
            File newFile = new File(fastaFileName);
            FileReader fileReader = new FileReader(newFile);
            BufferedReader reader = new BufferedReader(fileReader);
            while ((line = reader.readLine()) != null) {
                if (line.length() > 0) {
                    char firstChar = line.charAt(0);
                    if (firstChar == delim.charAt(0)) {
                        headerLines.append(line);
                        seqReadStart = "TRUE";
                    } else {
                        if (seqReadStart == "TRUE") {
                            seqLines.append(delim);
			}
			seqLines.append(line);
			seqReadStart = "FALSE";
                    }
		}
            }
            reader.close();
	} catch (Exception ex) {
            ex.printStackTrace();
	}
	String headerString = headerLines.toString();
	String seqString = seqLines.toString();
		
	headerString = headerString.substring(1,(headerString.length())); //remove opening '>'
	seqString = seqString.substring(1,(seqString.length()));
	String[] seqArray = seqString.split(delim);
	String[] headerArray = headerString.split(delim);
		
	SeqAlign alignment = new SeqAlign();
	alignment.makeAlignmentNoGaps(seqArray, headerArray);
	return alignment;
    }
}
