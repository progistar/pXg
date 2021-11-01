package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class CSVtoTSV {

	public static void main(String[] args) throws IOException {
		String fileName = args[0];
		
		BufferedReader BR = new BufferedReader(new FileReader(fileName));
		BufferedWriter BW = new BufferedWriter(new FileWriter(fileName.replace(".csv", ".tsv")));
		
		String line = null;
		
		int printShort = 10;
		while((line = BR.readLine()) != null) {
			line = line.replace(",", "\t");
			BW.append(line);
			BW.newLine();
			
			if(printShort-- > 0) System.out.println(line);
		}
		
		
		BW.close();
		BR.close();
	}
}
