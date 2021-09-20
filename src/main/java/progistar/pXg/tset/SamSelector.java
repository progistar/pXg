package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class SamSelector {
	
	public static void main(String[] args) throws IOException {
		writeSam("C:\\Users\\progi\\Desktop\\Projects\\pXg\\chr1toy.sam", 2);
	}
	
	/**
	 * Select some lines of records from a sam file. <br>
	 * Use for testing purpose. <br>
	 * 
	 * @param samFilePath
	 * @param numOfRecords
	 */
	public static void writeSam (String samFilePath, int numOfRecords) {
		try {
			File samFile = new File(samFilePath);
			File outSamFile = new File(samFile.getAbsolutePath().replace(".sam", ".sub"+numOfRecords+".sam"));
			
			BufferedReader BR = new BufferedReader(new FileReader(samFile));
			BufferedWriter BW = new BufferedWriter(new FileWriter(outSamFile));
			String line = null;
			
			while((line = BR.readLine()) != null) {
				if(numOfRecords == 0) break;
				
				// meta data is not counted.
				if(!line.startsWith("@")) {
					if(!line.contains("HISEQ:111:C6ADFANXX:1:2101:10310:56482")) continue;
					numOfRecords--;
				}
				
				BW.append(line);
				BW.newLine();
			}
			
			BR.close();
			BW.close();
			
		}catch (IOException ioe) {
			
		}
	}

}
