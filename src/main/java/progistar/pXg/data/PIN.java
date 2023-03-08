package progistar.pXg.data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import progistar.pXg.constants.Parameters;

public class PIN {

	private static String PIN_HEADER = "SpecID\tLabel\tScanNr\tMainScore\tLog2Reads";
	private static String[] pXgADDED_HEADERS = {"UniqueID", "Label"};
	private static String[] pXg_DEFAULT_FEATURES = {"Reads", "Rank", "InferredPeptide"};
	
	private PIN() {};
	
	public static void parseOutput () {
		try {
			File pXgOutput = new File(Parameters.outputFilePath);
			
			BufferedReader BR = new BufferedReader(new FileReader(pXgOutput));
			BufferedWriter BW = new BufferedWriter(new FileWriter(Parameters.pinFilePath));
			ArrayList<String> records = new ArrayList<String>();
			String[] headerFields = BR.readLine().split("\t");
			int indexShiftSize = pXgADDED_HEADERS.length;
			
			// find min-max charge
			String line = null;
			int minCharge = 100;
			int maxCharge = 0;
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				int charge = Integer.parseInt(fields[Parameters.chargeColumnIndex + indexShiftSize]);
				minCharge = Math.min(minCharge, charge);
				maxCharge = Math.max(maxCharge, charge);
				records.add(line);
			}
			
			BR.close();
			

			// find InferredPeptide index
			// find Rank index
			// find Reads
			int[] pXgDefaultFeatIdices = new int[pXg_DEFAULT_FEATURES.length];
			for(int fIdx=0; fIdx<pXg_DEFAULT_FEATURES.length; fIdx++) {
				for(int i=0; i<headerFields.length; i++) {
					if(headerFields[i].equalsIgnoreCase(pXg_DEFAULT_FEATURES[fIdx])) {
						pXgDefaultFeatIdices[fIdx] = i;
						break;
					}
				}
			}
			
			// add charge header
			for(int charge=minCharge; charge<=maxCharge; charge++) {
				PIN_HEADER += "\tCharge"+charge;
			}
			

			// find additional feature index
			if(Parameters.additionalFeatureIndices != null) {
				StringBuilder additionalFeatureHeader = new StringBuilder();
				for(int i=0; i<Parameters.additionalFeatureIndices.length; i++) {
					additionalFeatureHeader.append("\t").append(headerFields[Parameters.additionalFeatureIndices[i] + indexShiftSize]);
				}
				PIN_HEADER += additionalFeatureHeader.toString();
			}
			
			// last header
			PIN_HEADER += "\tLength\tRank\tPeptide\tProteins";
			
			BW.append(PIN_HEADER);
			BW.newLine();
			
			StringBuilder pinOutput = new StringBuilder();
			for(String record : records) {
				// init pinoutout string
				pinOutput.setLength(0);
				
				String[] fields = record.split("\t");
				
				String specId = fields[0];
				String label = fields[1];
				String scanNr = fields[Parameters.scanColumnIndex + indexShiftSize];
				String mainScore = fields[Parameters.scoreColumnIndex + indexShiftSize];
				String log2Reads = "" + Math.log(Double.parseDouble(fields[pXgDefaultFeatIdices[0]])+1)/Math.log(2);
				int charge = Integer.parseInt(fields[Parameters.chargeColumnIndex + indexShiftSize]);
				
				pinOutput.append(specId+"\t"+label+"\t"+scanNr+"\t"+mainScore+"\t"+log2Reads);
				// append charge
				for(int c=minCharge; c<=maxCharge; c++) {
					if(c == charge) {
						pinOutput.append("\t1");
					} else {
						pinOutput.append("\t0");
					}
				}
				
				String rank = fields[pXgDefaultFeatIdices[1]];
				String peptide = fields[pXgDefaultFeatIdices[2]];
				int length = peptide.length();
				
				pinOutput.append("\t").append(length);
				pinOutput.append("\t").append(rank);
				pinOutput.append("\t").append(peptide);
				// target or decoy
				if(label.equalsIgnoreCase("1")) {
					pinOutput.append("\tTarget");
				} else {
					pinOutput.append("\tXXX_Decoy");
				}
				
				BW.append(pinOutput.toString());
				BW.newLine();
			}
			
			BW.close();
			
		}catch(IOException ioe) {
			
		}
	}
}
