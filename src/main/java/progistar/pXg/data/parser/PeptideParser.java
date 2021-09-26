package progistar.pXg.data.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import progistar.pXg.constants.Parameters;
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
					PeptideAnnotation.setFields(line.split("\t"));
				} 
				// record
				else {
					String[] record = line.split("\t");
					String peptide = record[Parameters.peptideColumnIndex];
					
					// find peptide strip sequence
					Matcher matcher = peptideRegExr.matcher(peptide);
					
					while(matcher.find()) {
						pSeq.append(matcher.group());
					}
					
					if(pSeq.length() >= Parameters.minPeptLen && pSeq.length() <= Parameters.maxPeptLen) {
						PBlock pBLock = new PBlock(record, pSeq.toString());
						PeptideAnnotation.pBlocks.add(pBLock);
					}
					pSeq.setLength(0);
					
				}
				
				recordCount ++;
				
			}
			
			BR.close();
			
		}catch (IOException ioe) {
			
		}
		

		long endTime = System.currentTimeMillis();
		
		// filter PSMs by delta-score threshold
		PeptideAnnotation.filter();
		System.out.println("\tElapsed time: "+((endTime-startTime)/1000) + " sec");

		
		// build keyword-trie
		PeptideAnnotation.buildKeywordTrie();
	}
}
