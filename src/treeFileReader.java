import java.io.*;

public class treeFileReader
{
	public static String readTree(String treeFileName)
	{
		String treeLine = null;
		StringBuilder allLines = new StringBuilder();
		try {
			File treeFile = new File(treeFileName);
			FileReader fileReader = new FileReader(treeFile);
			BufferedReader reader = new BufferedReader(fileReader);
			while ((treeLine = reader.readLine()) != null) {
				allLines.append(treeLine);
			}
			reader.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		return allLines.toString();
	}
}
