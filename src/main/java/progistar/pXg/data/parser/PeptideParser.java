package progistar.pXg.data.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import progistar.pXg.constants.Parameters;
import progistar.pXg.constants.RunInfo;
import progistar.pXg.data.PBlock;
import progistar.pXg.data.PeptideAnnotation;

public class PeptideParser {

	private static Pattern	peptideRegExr;
	private static String[]	commentMarkers;
	
	private PeptideParser () {}
	
	/**
	 * The peptideFile must contain field (column names). <br>
	 * Header lines are optional and if exist, then the lines must be positioned on the top. <br>
	 * 
	 * 
	 * @param peptideFilePath
	 */
	public static void parseResult (String peptideFilePath) {
		System.out.print("Parsing peptide file: "+peptideFilePath);
		long startTime = System.currentTimeMillis();
		
		// set regular expressions
		peptideRegExr = Pattern.compile(Parameters.peptideParserRegExr);
		commentMarkers = Parameters.commentMarker.split("\\|");
		
		StringBuilder pSeq = new StringBuilder();
		try {
			File file = new File(peptideFilePath);
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;

			int recordCount = -1;
			while((line = BR.readLine()) != null) {
				// skip header marker
				// comment marker is not considered record.
				for(String headerMarker : commentMarkers) {
					if(line.startsWith(headerMarker)) continue;
				}
				
				// the first line after headers must be field line.
				if(recordCount == -1) {
					if(Parameters.sepType.equalsIgnoreCase("tsv")) {
						PeptideAnnotation.setFields(line.split("\t"));
					} else if(Parameters.sepType.equalsIgnoreCase("csv")) {
						PeptideAnnotation.setFields(line.split(","));
					}
				} 
				// record
				else {
					String[] record = null;
					if(Parameters.sepType.equalsIgnoreCase("tsv")) {
						record = line.split("\t");
					} else if(Parameters.sepType.equalsIgnoreCase("csv")) {
						record = line.split(",");
					}
					String peptide = record[Parameters.peptideColumnIndex];
					
					// find peptide strip sequence
					Matcher matcher = peptideRegExr.matcher(peptide);
					
					while(matcher.find()) {
						pSeq.append(matcher.group());
					}
					
					PBlock pBlock = new PBlock(record, pSeq.toString());
					
					PeptideAnnotation.pBlocks.add(pBlock);
					pSeq.setLength(0);
					
				}
				
				recordCount ++;
				
			}
			
			BR.close();
			
		}catch (IOException ioe) {
			
		}
		

		long endTime = System.currentTimeMillis();
		
		RunInfo.initialPeptideNum = PeptideAnnotation.getPeptideSize();
		RunInfo.initialScanNum = PeptideAnnotation.getScanSize();
		
		// filter PSMs by delta-rank
		PeptideAnnotation.filter();
		
		RunInfo.rankFilterPeptideNum1 = PeptideAnnotation.getPeptideSize();
		RunInfo.rankFilterScanNum1 = PeptideAnnotation.getScanSize();
		
		// length filter
		PeptideAnnotation.peptideLengthFilter();
		
		RunInfo.lengthFilterPeptideNum2 = PeptideAnnotation.getPeptideSize();
		RunInfo.lengthFilterScanNum2 = PeptideAnnotation.getScanSize();
		
		System.out.println("\tElapsed time: "+((endTime-startTime)/1000) + " sec");

		
		// build keyword-trie
		PeptideAnnotation.buildKeywordTrie();
	}
}
