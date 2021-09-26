package progistar.pXg.data.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import progistar.pXg.constants.Parameters;

public class ParameterParser {

	public static void parseParams (String parameterFilePath) {
		File file = new File(parameterFilePath);
		
		try {
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			// TODO:
			// check all parameters are properly taken.
			// if indexed files already exist, then skip the indexing process. 
			
			while((line = BR.readLine()) != null) {
				// skip comment
				if(line.startsWith("#")) continue;
				if(line.length() == 0) continue;
				
				if(line.contains("=")) {
					String param = line.split("\\=")[0].trim();
					String value = line.split("\\=")[1].trim();
					
					if(param.equalsIgnoreCase(Parameters.GENOMIC_ANNOTATION_PATH)) {
						Parameters.genomicAnnotationFilePath = value;
					}
					
					else if(param.equalsIgnoreCase(Parameters.GENOMIC_SEQUENCE_PATH)) {
						Parameters.sequenceFilePath = value;
					}
					
					else if(param.equalsIgnoreCase(Parameters.PEPTIDE_ANNOTATION_PATH)) {
						Parameters.peptideFilePath = value;
					}
					
					else if(param.equalsIgnoreCase(Parameters.PEPTIDE_COLUMN_INDEX)) {
						// For the purpose of user-friendly, it takes one-based and converts to zero-based. 
						Parameters.peptideColumnIndex = Integer.parseInt(value) - 1;
					}
					
					else if(param.equalsIgnoreCase(Parameters.pParserRegExr)) {
						Parameters.peptideParserRegExr = value;
					}
					
					else if(param.equalsIgnoreCase(Parameters.cMarker)) {
						Parameters.commentMarker = value;
					}
					
					else if(param.equalsIgnoreCase(Parameters.GENOMIC_ANNOTATION_PARTITION_SIZE)) {
						Parameters.partitionSize = Integer.parseInt(value);
					}
					
					else if(param.equalsIgnoreCase(Parameters.GENOMIC_SEQUENCE_PARTITION_SIZE)) {
						Parameters.readSize = Integer.parseInt(value);
					}
					
					else if(param.equalsIgnoreCase(Parameters.numOfThreads)) {
						Parameters.nThreads = Integer.parseInt(value);
					}
					
					else if(param.equalsIgnoreCase(Parameters.OUTPUT_PATH)) {
						Parameters.outputFilePath = value;
					}
					
					else if(param.equalsIgnoreCase(Parameters.SCAN_COLUMN_INDICES)) {
						String[] indicies = value.split("\\,");
						int[] scanIndicies = new int[indicies.length];
						for(int i=0; i<indicies.length; i++) {
							// For the purpose of user-friendly, it takes one-based and converts to zero-based.
							scanIndicies[i] = Integer.parseInt(indicies[i])-1;
						}
						Parameters.scanColumnIndices = scanIndicies;
					}
					
					else if(param.equalsIgnoreCase(Parameters.SCORE_COLUMN_INDEX)) {
						// For the purpose of user-friendly, it takes one-based and converts to zero-based.
						Parameters.scoreColumnIndex = Integer.parseInt(value);
					}
					
					else if(param.equalsIgnoreCase(Parameters.DELTA_SCORE_THRESHOLD)) {
						Parameters.deltaScoreThreshold = Double.parseDouble(value);
					}
				}
				
			}
			
			BR.close();
		} catch (IOException ioe) {
			
		}
	}
	
	private static boolean isAlreadyIndexed () {
		boolean isAlreadyIndexed = false;
		
		File file = new File(Parameters.sequenceFilePath);
		
		
		
		return isAlreadyIndexed;
	}
}
