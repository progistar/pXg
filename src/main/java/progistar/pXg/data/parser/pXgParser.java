package progistar.pXg.data.parser;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
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
		
		if(Parameters.isStringent) {
			System.out.println("Exclude non-canonical peptides with FastaIDs...");
			ArrayList<pXgRecord> selectedRecords = new ArrayList<pXgRecord>();
			for(pXgRecord record : records) {
				if(record.isCanonical() || !record.hasFastaID()) {
					selectedRecords.add(record);
				}
			}
			
			int excluded = records.size() - selectedRecords.size();
			noncanonical -= excluded;
			System.out.println(excluded +" non-canonical peptides were excluded");
			records = selectedRecords;
		}
		
		System.out.println("A total of "+records.size()+" unique entries were saved (canonical: "+canonical+", non-canonical:"+noncanonical+")");
		
		long endTime = System.currentTimeMillis();
		System.out.println("Elapsed time: "+(endTime - startTime)/1000 +" sec");
		System.out.println();
		return records;
	}
	
	/**
	 * Merge all temporary pXg output files into a single one.<br>
	 * 
	 * 
	 * @param tmpRecords
	 * @param outputFilePath
	 */
	public static void writeMergedResult (ArrayList<pXgRecord>[] tmpRecords, String outputFilePath) {
		String[] samFileNames = new String[tmpRecords.length];
		Hashtable<String, Hashtable<String, pXgRecord[]>> spectrumToKey = new Hashtable<String, Hashtable<String, pXgRecord[]>>();
		
		for(int i=0; i<samFileNames.length; i++) {
			samFileNames[i] = Parameters.tmpOutputFilePaths[i].replace("."+Constants.UNIQUE_RUN_ID, "");
			
			ArrayList<pXgRecord> records = tmpRecords[i];
			
			for(int j=0; j<records.size(); j++) {
				pXgRecord record = records.get(j);
				String specId = record.getValueByFieldName("SpecID");
				Hashtable<String, pXgRecord[]> keyToRecords = spectrumToKey.get(specId);
				if(keyToRecords == null) {
					keyToRecords = new Hashtable<String, pXgRecord[]>();
					spectrumToKey.put(specId, keyToRecords);
				}
				
				String key = record.getID();
				
				pXgRecord[] recordArray = keyToRecords.get(key);
				if(recordArray == null) {
					recordArray = new pXgRecord[samFileNames.length];
					keyToRecords.put(key, recordArray);
				}
				recordArray[i] = record;
			}
		}
		ArrayList<String> finalResults = new ArrayList<String>();
		StringBuilder tmp = new StringBuilder();
		
		for(int i=0; i<pXgParser.header.length; i++) {
			if(i != 0) {
				tmp.append("\t");
			}
			tmp.append(pXgParser.header[i]);
		}
		
		if(samFileNames.length > 1) {
			for(int i=0; i<samFileNames.length; i++) {
				tmp.append("\t").append(new File(samFileNames[i]).getName());
			}
		}
		
		String header = tmp.toString();
		
		spectrumToKey.forEach((specId, keyToRecords)->{
			keyToRecords.forEach((key, array)->{
				StringBuilder counts = new StringBuilder();
				int mergedReads = 0;
				pXgRecord maxRecord = null;
				int maxReads = 0;
				
				for(int i=0; i<array.length; i++) {
					int reads = 0;
					if(array[i] != null) {
						reads = Integer.parseInt(array[i].getValueByFieldName("Reads"));
						mergedReads += reads;
						if(reads > maxReads) {
							maxRecord = array[i];
							maxReads = reads;
						}
					}
					counts.append("\t"+reads);
				}
				
				// change representative reads to merged reads.
				maxRecord.setValueByFieldName("Reads", mergedReads+"");
				
				tmp.setLength(0);
				tmp.append(maxRecord.toString()).append(counts.toString());
				finalResults.add(tmp.toString());
			});
		});
		
		try {
			BufferedWriter BW = new BufferedWriter(new FileWriter(outputFilePath));
			
			Collections.sort(finalResults, new Comparator<String>() {
				
				public int compare (String s1, String s2) {
					String specId1 = s1.split("\t")[0];
					String specId2 = s2.split("\t")[0];
					int cVal = specId1.compareTo(specId2);
					if(cVal > 0) {
						return 1;
					} else if(cVal < 0) {
						return -1;
					}
					return 0;
				}
				
			});
			
			BW.append(header);
			BW.newLine();
			for(String record : finalResults) {
				BW.append(record);
				BW.newLine();
			}
			
			BW.close();
		}catch(IOException ioe) {
			
		}
		
	}
}
