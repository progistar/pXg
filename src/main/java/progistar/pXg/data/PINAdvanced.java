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

public class PINAdvanced {

	private static String PIN_HEADER = "SpecId\tLabel\tScanNr\tMainScore\tLog2Reads";
	private static String[] pXgADDED_HEADERS = {"UniqueID", "Label"};
	private static String[] pXg_DEFAULT_FEATURES = {"Mutations","Events","IsCanonical","DeltaScore","Reads", "Rank", "InferredPeptide"};
	
	private PINAdvanced() {};
	
	/**
	 * Further analysis for group-specific characterization,
	 * Very specific level of events can be considered in future.
	 * 
	 * @param isStableMethod
	 */
	public static void parseOutput (boolean isStableMethod) {
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
			
			int mutationIdx = pXgDefaultFeatIdices[0];
			int eventIdx = pXgDefaultFeatIdices[1];
			int isCanonicalIdx = pXgDefaultFeatIdices[2];
			int deltaScoreIdx = pXgDefaultFeatIdices[3];
			int readIdx = pXgDefaultFeatIdices[4];
			int rankIdx = pXgDefaultFeatIdices[5];
			int infPeptIdx = pXgDefaultFeatIdices[6];
			
			// to adjust index caused by appending "UniqueID" and "Label" to the original input,
			// the original index must be shifted by 2.
			int indexShiftSize = pXgADDED_HEADERS.length;
			
			// find min-max charge
			String line = null;
			int minCharge = 100;
			int maxCharge = 0;
			ArrayList<String> eventSegments = new ArrayList<String>();
			Hashtable<String, String> eventSegmentCheck = new Hashtable<String, String>();
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
				
				// figure out event segments
				// such as AS, sense, PC etc...
				String[] events = fields[eventIdx].split("\\|");
				for(String event : events) {
					String[] segments = event.split("\\;");
					
					for(String segment : segments) {
						if(eventSegmentCheck.get(segment) == null) {
							eventSegmentCheck.put(segment, "");
							
							// ban event
							// sense is duplicated information to asRNA
							if(!segment.equalsIgnoreCase(Constants.EVENT_SENSE)) {
								eventSegments.add(segment);
							}
						}
					}
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
			
			// mutation header
			// # of SNV and INDEL
			PIN_HEADER += "\tSNV\tINDEL";
			
			// event segment
			for(String eventSegment : eventSegments) {
				PIN_HEADER += "\t"+eventSegment;
			}
			
			// last header
			PIN_HEADER += "\tIsCanonical\tDeltaScore\tLength\tRank\tPeptide\tProteins";
			
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
				String[] mutations = fields[mutationIdx].split("\\|");
				String[] events = fields[eventIdx].split("\\|");
				String rank = fields[rankIdx];
				String peptide = fields[infPeptIdx];
				String deltaScore = fields[deltaScoreIdx];
				int length = peptide.length();
				
				/// mutation
				int snvCnt = 0;
				int indelCnt = 0;
				for(String mutation : mutations) {
					if(mutation.contains(">")) {
						snvCnt++;
					} else if(mutation.contains("del") || mutation.contains("ins")) {
						indelCnt++;
					}
				}
				pinOutput.append("\t").append(snvCnt);
				pinOutput.append("\t").append(indelCnt);
				
				// event
				Hashtable<String, String> thisEventSegments = new Hashtable<String, String>();
				for(String event : events) {
					String[] segments = event.split("\\;");
					for(String segment : segments) {
						thisEventSegments.put(segment, "");
					}
				}
				for(String eventSegment : eventSegments) {
					if(thisEventSegments.get(eventSegment) == null) {
						pinOutput.append("\t").append(0);
					} else {
						pinOutput.append("\t").append(1);
					}
				}
				
				// isCanonical
				if(fields[isCanonicalIdx].equalsIgnoreCase("TRUE")) {
					pinOutput.append("\t1");
				} else {
					pinOutput.append("\t-1");
				}
				
				pinOutput.append("\t").append(deltaScore);
				pinOutput.append("\t").append(length);
				pinOutput.append("\t").append(rank);
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
