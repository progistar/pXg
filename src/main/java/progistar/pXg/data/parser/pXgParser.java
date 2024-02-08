package progistar.pXg.data.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import progistar.pXg.data.pXgRecord;

public class pXgParser {

	public static String[] header = null;
	
	private pXgParser() {}
	
	public static ArrayList<pXgRecord> parse (File file) throws IOException {
		long startTime = System.currentTimeMillis();
		System.out.println("Parsing "+file.getName());
		
		ArrayList<pXgRecord> records = new ArrayList<pXgRecord>();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		pXgParser.header = BR.readLine().split("\t");
		String line = null;
		
		Hashtable<String, String> checkDuplicates = new Hashtable<String, String>();
		int totalPRSM = 0;
		int canonical = 0;
		int noncanonical = 0;
		while((line = BR.readLine()) != null) {
			totalPRSM++;
			pXgRecord record = new pXgRecord(line.split("\t"));
			String header = record.getHeader();
			
			if(checkDuplicates.get(header) == null) {
				records.add(record);
				checkDuplicates.put(header, "");
				
				if(record.isCanonical()) {
					canonical ++;
				} else {
					noncanonical ++;
				}
			}
		}
		
		BR.close();
		
		System.out.println("A total of "+totalPRSM+" PRSMs were parsed");
		System.out.println("A total of "+checkDuplicates.size()+" unique entries were saved (canonical: "+canonical+", non-canonical:"+noncanonical+")");
		
		long endTime = System.currentTimeMillis();
		System.out.println("Elapsed time: "+(endTime - startTime)/1000 +" sec");
		System.out.println();
		return records;
	}
}
