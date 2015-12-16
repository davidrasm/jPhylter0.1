import java.io.*;
import cern.colt.matrix.*;
import java.util.ArrayList;

public class WriteMatrixToCSV
{
    
//DoubleMatrix2D dataMatrix;

	public void writeCSV(DoubleMatrix2D data, String fileName) throws IOException
	{
		FileWriter fw = new FileWriter(fileName);
		PrintWriter pw = new PrintWriter(fw);
		
		//Mock data
		//DoubleFactory2D factory2D;
		//factory2D = DoubleFactory2D.dense;
		//int size = 4;
		//DoubleMatrix2D dataMatrix = factory2D.make(size, size);
		//dataMatrix.assign(0.0);
		//System.out.println(dataMatrix);
		DoubleMatrix2D dataMatrix = data.copy();
                
		int nrows = dataMatrix.rows();
		int ncols = dataMatrix.columns();
		
		for (int j = 0; j < nrows; ++j) {
			for(int k = 0; k < (ncols - 1); ++k) {
				double currDouble = dataMatrix.get(j,k);
				String currString = Double.toString(currDouble);
				pw.print(currString);
				pw.print(",");
			}
			double currDouble = dataMatrix.get(j,(ncols-1));
			String currString = Double.toString(currDouble);
			pw.print(currString);
			pw.println();
		}
		
		pw.flush();
    	pw.close();
    	fw.close();
	}
        
    public static DoubleMatrix2D readMatrixFromCSV(String fileName) throws IOException 
    {    
        
        //String fileName = "thisFile.csv";
        BufferedReader reader = new BufferedReader(new FileReader(fileName));
        
        ArrayList<ArrayList<Double>> data = new ArrayList<ArrayList<Double>>();
        
        String dataRow = reader.readLine(); //Read first line
        int row = 0;
        while (dataRow != null) {
            String[] dataArray = dataRow.split(",");
            data.add(new ArrayList<Double>());
            for (String item:dataArray) {
                data.get(row).add(Double.valueOf(item));
                //System.out.print(item + "\t");
            }
            //System.out.println();
            dataRow = reader.readLine();
            row++;
        }
        reader.close();
        
        int rows = data.size();
        int cols = data.get(0).size();
        DoubleFactory2D factory2D = DoubleFactory2D.dense;
        DoubleMatrix2D dataMatrix = factory2D.make(rows,cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                dataMatrix.setQuick(i,j,data.get(i).get(j));
            }
        }
        
        return dataMatrix;
        
    }
        
}
