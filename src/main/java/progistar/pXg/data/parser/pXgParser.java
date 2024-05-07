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

	public static ArrayList<pXgRecord> parse (File file, boolean removeDuplicates) throws IOException {
		long startTime = System.currentTimeMillis();
		System.out.println("Parsing "+file.getName());

		ArrayList<pXgRecord> records = new ArrayList<>();

		BufferedReader BR = new BufferedReader(new FileReader(file));
		pXgParser.header = BR.readLine().split("\t");
		String line = null;

		Hashtable<String, String> checkDuplicates = new Hashtable<>();
		int totalPRSM = 0;
		int canonical = 0;
		int noncanonical = 0;
		while((line = BR.readLine()) != null) {
			totalPRSM++;
			pXgRecord record = new pXgRecord(line.split("\t"));
			String header = record.getHeader();

			if(!removeDuplicates || checkDuplicates.get(header) == null) {
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
			ArrayList<pXgRecord> selectedRecords = new ArrayList<>();
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
		
		if(removeDuplicates) {
			System.out.println("A total of "+records.size()+" unique entries were saved (canonical: "+canonical+", non-canonical:"+noncanonical+")");
		} else {
			System.out.println("A total of "+records.size()+" entries were saved (canonical: "+canonical+", non-canonical:"+noncanonical+")");
		}

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
		Hashtable<String, Hashtable<String, pXgRecord[]>> spectrumToKey = new Hashtable<>();

		for(int i=0; i<samFileNames.length; i++) {
			samFileNames[i] = Parameters.tmpOutputFilePaths[i].replace("."+Constants.UNIQUE_RUN_ID, "");

			ArrayList<pXgRecord> records = tmpRecords[i];

			for (pXgRecord record : records) {
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
		ArrayList<String> finalResults = new ArrayList<>();
		StringBuilder tmp = new StringBuilder();

		for(int i=0; i<pXgParser.header.length; i++) {
			if(i != 0) {
				tmp.append("\t");
			}
			tmp.append(pXgParser.header[i]);
		}

		if(samFileNames.length > 1) {
			for (String samFileName : samFileNames) {
				tmp.append("\t").append(new File(samFileName).getName());
			}
		}

		String header = tmp.toString();
		
		// The combination of genomic loci + strand + nucleotide sequence is mapping to unique gneomic ID.
		Hashtable<String, Integer> genomicIDMapper = new Hashtable<>();
		
		spectrumToKey.forEach((specId, keyToRecords)->{
			keyToRecords.forEach((key, array)->{
				StringBuilder counts = new StringBuilder();
				int mergedReads = 0;
				pXgRecord maxRecord = null;
				int maxReads = 0;

				for (pXgRecord element : array) {
					int reads = 0;
					if(element != null) {
						reads = Integer.parseInt(element.getValueByFieldName("Reads"));
						mergedReads += reads;
						if(reads > maxReads) {
							maxRecord = element;
							maxReads = reads;
						}
					}
					counts.append("\t"+reads);
				}

				Integer genomicID = genomicIDMapper.get(key);
				if(genomicID == null) {
					genomicID = genomicIDMapper.size()+1;
					genomicIDMapper.put(key, genomicID);
				}
				
				// change representative reads to merged reads.
				maxRecord.setValueByFieldName("Reads", mergedReads+"");
				maxRecord.setValueByFieldName("GenomicID", genomicID+"");
				
				tmp.setLength(0);
				tmp.append(maxRecord.toString()).append(counts.toString());
				finalResults.add(tmp.toString());
			});
		});

		try {
			BufferedWriter BW = new BufferedWriter(new FileWriter(outputFilePath));

			Collections.sort(finalResults, new Comparator<String>() {

				@Override
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
