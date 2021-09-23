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
					
					if(param.equalsIgnoreCase(Parameters.gAnnPath)) {
						Parameters.genomicAnnotationFilePath = value;
					}
					
					else if(param.equalsIgnoreCase(Parameters.gSeqPath)) {
						Parameters.sequenceFilePath = value;
					}
					
					else if(param.equalsIgnoreCase(Parameters.pAnnPath)) {
						Parameters.peptideFilePath = value;
					}
					
					else if(param.equalsIgnoreCase(Parameters.pColumnIndex)) {
						// For the purpose of user-friendly, it takes one-based and converts to zero-based. 
						Parameters.peptideColumnIndex = Integer.parseInt(value) - 1;
					}
					
					else if(param.equalsIgnoreCase(Parameters.pParserRegExr)) {
						Parameters.peptideParserRegExr = value;
					}
					
					else if(param.equalsIgnoreCase(Parameters.cMarker)) {
						Parameters.commentMarker = value;
					}
					
					else if(param.equalsIgnoreCase(Parameters.gPartitionSize)) {
						Parameters.partitionSize = Integer.parseInt(value);
					}
					
					else if(param.equalsIgnoreCase(Parameters.gSeqReadSize)) {
						Parameters.readSize = Integer.parseInt(value);
					}
					
					else if(param.equalsIgnoreCase(Parameters.numOfThreads)) {
						Parameters.nThreads = Integer.parseInt(value);
					}
					
					else if(param.equalsIgnoreCase(Parameters.oPath)) {
						Parameters.outputFilePath = value;
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
