package progistar.pXg.data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;

public class PIN {

	/**
	 * 
	 * ScanNr is actually acting as scan index
	 */
	private static String PIN_HEADER = "SpecId\tLabel\tScanNr\tMainScore\tLog2Reads";
	private static String[] pXgADDED_HEADERS = {"UniqueID", "Label"};
	private static String[] pXg_DEFAULT_FEATURES = {"DeltaScore","Reads","MeanQScore", "InferredPeptide"};
	
	private PIN() {};
	
	/**
	 * Further analysis for group-specific characterization,
	 * Very specific level of events can be considered in future.
	 * 
	 * @param isStableMethod
	 */
	public static void parseOutput () {
		try {
			File pXgOutput = new File(Parameters.outputFilePath);
			ArrayList<String> pinRecords = new ArrayList<String>();
			BufferedReader BR = new BufferedReader(new FileReader(pXgOutput));
			BufferedWriter BW = new BufferedWriter(new FileWriter(Parameters.pinFilePath));
			Hashtable<String, Integer> specIDtoScanIdx = new Hashtable<String, Integer>();
			ArrayList<String> records = new ArrayList<String>();
			String[] headerFields = BR.readLine().split("\t");
			
			/// Find pXg default features
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
			
			int deltaScoreIdx = pXgDefaultFeatIdices[0];
			int readIdx = pXgDefaultFeatIdices[1];
			int meanQScoreIdx = pXgDefaultFeatIdices[2];
			int infPeptIdx = pXgDefaultFeatIdices[3];
			
			// to adjust index caused by appending "UniqueID" and "Label" to the original input,
			// the original index must be shifted by 2.
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
				
				// scan idx gen
				String specID = fields[0];
				if(specIDtoScanIdx.get(specID) == null) {
					specIDtoScanIdx.put(specID, specIDtoScanIdx.size()+1);
				}
			}
			
			BR.close();
			
			/******************8 HEADER maker 8*******************/
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
			PIN_HEADER += "\tDeltaScore\tMeanQScore\tPeptide\tProteins";
			
			pinRecords.add(PIN_HEADER);
			/******************8 Gen PIN 8*******************/
			
			StringBuilder pinOutput = new StringBuilder();
			for(String record : records) {
				// init pinoutout string
				pinOutput.setLength(0);
				
				String[] fields = record.split("\t");
				
				String specId = fields[0];
				String label = fields[1];
				String scanNr = specIDtoScanIdx.get(specId)+"";
				String mainScore = fields[Parameters.scoreColumnIndex + indexShiftSize];
				String log2Reads = "" + Math.log(Double.parseDouble(fields[readIdx])+1)/Math.log(2);
				
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
				
				// additional features
				if(Parameters.additionalFeatureIndices != null) {
					for(int i=0; i<Parameters.additionalFeatureIndices.length; i++) {
						pinOutput.append("\t").append(fields[Parameters.additionalFeatureIndices[i] + indexShiftSize]);
					}
				}
				
				// pXg default features
				String deltaScore = fields[deltaScoreIdx];
				String meanQScore = fields[meanQScoreIdx];
				String peptide = fields[infPeptIdx];
				
				pinOutput.append("\t").append(deltaScore);
				pinOutput.append("\t").append(meanQScore);
				pinOutput.append("\t").append(peptide);
				// target or decoy
				if(label.equalsIgnoreCase("1")) {
					pinOutput.append("\tTarget");
				} else {
					pinOutput.append("\trandom_Decoy");
				}
				
				pinRecords.add(pinOutput.toString());
			}
			
			
			for(int i=0; i<pinRecords.size(); i++) {
				BW.append(pinRecords.get(i));
				BW.newLine();
			}
			
			BW.close();
			
		}catch(IOException ioe) {
			
		}
	}
	
}
