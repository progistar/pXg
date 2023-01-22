package progistar.pXg.data.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class pXgTSVParser {
	
	private pXgTSVParser () {
		// void
	}
	
	public static ArrayList<String> parseTSV (File file) {
		ArrayList<String> simpleRecords = new ArrayList<String>();
		try {
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			
			// read header
			String header = BR.readLine();
			
			// parsing records
			while((line = BR.readLine()) != null) {
				simpleRecords.add(line);
			}
			
			BR.close();
		}catch (IOException ioe) {
			
		}
		
		return simpleRecords;
	}
}
