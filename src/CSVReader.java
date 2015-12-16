import java.io.*;
import java.util.Arrays;
import java.util.ArrayList;

/**
 *
 * @author David
 */
public class CSVReader {
    
    public ArrayList<ArrayList<String>> readMatrixFromCSV(String fileName) throws IOException 
    {    
        
        //String fileName = "thisFile.csv";
        BufferedReader reader = new BufferedReader(new FileReader(fileName));
        
        ArrayList<ArrayList<String>> data = new ArrayList<ArrayList<String>>();
        String dataRow = reader.readLine(); //Read first line
        int row = 0;
        
        while (dataRow != null) {
            String[] dataArray = dataRow.split(",");
            data.add(new ArrayList<String>());
            for (String item:dataArray) {
                data.get(row).add(item);
                //System.out.print(item + "\t");
            }
            //System.out.println();
            dataRow = reader.readLine();
            row++;
        }
        reader.close();
        
        return data;
        
    }
    
}
