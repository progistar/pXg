package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import progistar.pXg.constants.Parameters;

class HEADER {
	
	String[] thisHeader = null;
	double[] lengthReadSigma = new double[20];
	double[] lengthReadMean = new double[20];
	ArrayList[] lengthReads = new ArrayList[20];
	
	int specIdIdx;
	int labelIdx;
	int scanNrIdx;
	int mainScoreIdx;
	int log2ReadsIdx;
	int deltaScoreIdx;
	ArrayList<Integer> chargeIndices = new ArrayList<Integer>();
	int ppmIdx;
	int peptIdx;
	int protIdx;
	

	int snvIdx = -1;
	int indelIdx = -1;
	int fsIdx = -1;
	int utr5Idx = -1;
	int utr3Idx = -1;
	int asRNAIdx = -1;
	int ncRNAIdx = -1;
	int IGRIdx = -1;
	int IRIdx = -1;
	int ASIdx = -1;
	int unknownIdx = -1;
	int SAIdx = -1;
	int bestDeltaRTIdx = -1;
	int baIdx = -1;
	int avgDeltaRTIdx = -1;
	int lenIdx = -1;
	
	public String genRecord1 (String[] fields) {
		StringBuilder sb = new StringBuilder();
		String specId = fields[specIdIdx];
		String label = fields[labelIdx];
		String scanNr = fields[scanNrIdx];
		String mainScore = fields[mainScoreIdx];
		String log2Reads = fields[log2ReadsIdx];
		String ppm = fields[ppmIdx];
		String pept = fields[peptIdx];
		String prot = fields[protIdx];
		String deltaScore = fields[deltaScoreIdx];
		ArrayList<String> charges = new ArrayList<String>();
		for(int i=0; i<chargeIndices.size(); i++) {
			charges.add(fields[chargeIndices.get(i)]);
		}
		
		double crypticScore = 0;
		int snv = 0;
		if(snvIdx != -1) {
			snv = Integer.parseInt(fields[snvIdx]);
		}
		int indel = 0;
		if(indelIdx != -1) {
			indel = Integer.parseInt(fields[indelIdx]);
		}
		int fs = 0;
		if(fsIdx != -1) {
			fs = Integer.parseInt(fields[fsIdx]);
		}
		int utr5 = 0;
		if(utr5Idx != -1) {
			utr5 = Integer.parseInt(fields[utr5Idx]);
		}
		int utr3 = 0;
		if(utr3Idx != -1) {
			utr3 = Integer.parseInt(fields[utr3Idx]);
		}
		int asRNA = 0;
		if(asRNAIdx != -1) {
			asRNA = Integer.parseInt(fields[asRNAIdx]);
		}
		int IGR = 0;
		if(IGRIdx != -1) {
			IGR = Integer.parseInt(fields[IGRIdx]);
		}
		int IR = 0;
		if(IRIdx != -1) {
			IR = Integer.parseInt(fields[IRIdx]);
		}
		int ncRNA = 0;
		if(ncRNAIdx != -1) {
			ncRNA = Integer.parseInt(fields[ncRNAIdx]);
		}
		int AS = 0;
		if(ASIdx != -1) {
			AS = Integer.parseInt(fields[ASIdx]);
		}
		
		int unknown = 0;
		if(unknownIdx != -1) {
			unknown = Integer.parseInt(fields[unknownIdx]);
		}
		
		double sa = 0;
		if(SAIdx != -1) {
			sa = Double.parseDouble(fields[SAIdx]);
		}
		
		double bestDeltaRT = 0;
		if(bestDeltaRTIdx != -1) {
			bestDeltaRT = Math.abs(Double.parseDouble(fields[bestDeltaRTIdx]));
		}
		
		double baScore = 0;
		
		if(baIdx != -1) {
			baScore = Double.parseDouble(fields[baIdx]);
		}
		
		double avgDeltaRT = 0;
		if(avgDeltaRTIdx != -1) {
			avgDeltaRT = Double.parseDouble(fields[avgDeltaRTIdx]);
		}
		
		String length = fields[lenIdx];
		
		
		crypticScore = (double) (snv * (Math.log(Parameters.PENALTY_MUTATION+1)/Math.log(2)) +
				indel * (Math.log(Parameters.PENALTY_MUTATION+1)/Math.log(2)) +
				fs * (Math.log(Parameters.PENALTY_FS+1)/Math.log(2)) +
				utr5 * (Math.log(Parameters.PENALTY_5UTR+1)/Math.log(2)) +
				utr3 * (Math.log(Parameters.PENALTY_3UTR+1)/Math.log(2)) +
				asRNA * (Math.log(Parameters.PENALTY_asRNA+1)/Math.log(2)) +
				IGR * (Math.log(Parameters.PENALTY_IGR+1)/Math.log(2)) +
				IR * (Math.log(Parameters.PENALTY_IR+1)/Math.log(2)) +
				ncRNA * (Math.log(Parameters.PENALTY_ncRNA+1)/Math.log(2)) +
				AS * (Math.log(Parameters.PENALTY_AS+1)/Math.log(2)) +
				unknown * (Math.log(Parameters.PENALTY_UNMAP+1)/Math.log(2)));

		if(crypticScore > 0) {
			crypticScore = -1;
		} else {
			crypticScore = 1;
		}

		if(baScore < 0) {
			return null;
		}
		
		crypticScore = 0;
		length = "0";
		//baScore = 0;
		avgDeltaRT = 0;
		//log2Reads = "0";
		/*
		Double log2ReadsNorm = Double.parseDouble(log2Reads);
		if(lengthReadSigma[pept.length()] != 0) {
			log2ReadsNorm = (Double.parseDouble(log2Reads) - lengthReadMean[pept.length()]) / lengthReadSigma[pept.length()];
		}
		log2Reads = log2ReadsNorm + "";
		*/
		//bestDeltaRT = 0;
		//sa = 0;
		
		sb.append(specId).append("\t").append(label).append("\t").append(scanNr).append("\t").append(mainScore)
		.append("\t").append(deltaScore).append("\t").append(log2Reads).append("\t").append(ppm);
		for(int i=0; i<charges.size(); i++) {
			sb.append("\t").append(charges.get(i));
		}
		sb.append("\t").append(crypticScore)
		.append("\t").append(sa)
		.append("\t").append(bestDeltaRT)
		.append("\t").append(avgDeltaRT)
		.append("\t").append(baScore)
		.append("\t").append(length)
		.append("\t").append(pept).append("\t").append(prot);
		
		return sb.toString();
	}
	
	public String toHeader1 () {
		StringBuilder sb = new StringBuilder();
		sb.append(thisHeader[specIdIdx]).append("\t")
		.append(thisHeader[labelIdx]).append("\t")
		.append(thisHeader[scanNrIdx]).append("\t")
		.append(thisHeader[mainScoreIdx]).append("\t")
		.append(thisHeader[deltaScoreIdx]).append("\t")
		.append(thisHeader[log2ReadsIdx]).append("\t")
		.append(thisHeader[ppmIdx]);
		for(int i=0; i<chargeIndices.size(); i++) {
			sb.append("\t").append(thisHeader[chargeIndices.get(i)]);
		}
		sb.append("\t").append("IsCanonical")
		.append("\t").append(thisHeader[SAIdx])
		.append("\t").append(thisHeader[bestDeltaRTIdx])
		.append("\t").append(thisHeader[avgDeltaRTIdx])
		.append("\t").append(thisHeader[baIdx])
		.append("\t").append(thisHeader[lenIdx])
		.append("\t").append(thisHeader[peptIdx])
		.append("\t").append(thisHeader[protIdx]);
		
		return sb.toString();
	}
	
	public HEADER (String[] headers) {
		thisHeader = headers;
		for(int i=0; i<headers.length; i++) {
			if(headers[i].equalsIgnoreCase("SNV")) {
				snvIdx = i;
			} else if(headers[i].equalsIgnoreCase("INDEL")) {
				indelIdx = i;
			}  else if(headers[i].equalsIgnoreCase("IR")) {
				IRIdx = i;
			}  else if(headers[i].equalsIgnoreCase("asRNA")) {
				asRNAIdx = i;
			}  else if(headers[i].equalsIgnoreCase("IGR")) {
				IGRIdx = i;
			}  else if(headers[i].equalsIgnoreCase("FS")) {
				fsIdx = i;
			}  else if(headers[i].equalsIgnoreCase("3`-UTR")) {
				utr3Idx = i;
			}  else if(headers[i].equalsIgnoreCase("5`-UTR")) {
				utr5Idx = i;
			} else if(headers[i].equalsIgnoreCase("AS")) {
				ASIdx = i;
			} else if(headers[i].equalsIgnoreCase("unknown")) {
				unknownIdx = i;
			} else if(headers[i].equalsIgnoreCase("ncRNA")) {
				ncRNAIdx = i;
			} else if(headers[i].equalsIgnoreCase("specId")) {
				specIdIdx = i;
			} else if(headers[i].equalsIgnoreCase("Label")) {
				labelIdx = i;
			} else if(headers[i].equalsIgnoreCase("ScanNr")) {
				scanNrIdx = i;
			} else if(headers[i].equalsIgnoreCase("mainScore")) {
				mainScoreIdx = i;
			} else if(headers[i].equalsIgnoreCase("log2Reads")) {
				log2ReadsIdx = i;
			} else if(headers[i].equalsIgnoreCase("ppm")) {
				ppmIdx = i;
			} else if(headers[i].toLowerCase().startsWith("charge")) {
				chargeIndices.add(i);
			} else if(headers[i].equalsIgnoreCase("peptide")) {
				peptIdx = i;
			} else if(headers[i].equalsIgnoreCase("proteins")) {
				protIdx = i;
			} else if(headers[i].equalsIgnoreCase("deltaScore")) {
				deltaScoreIdx = i;
			} else if(headers[i].equalsIgnoreCase("SA")) {
				SAIdx = i;
			} else if(headers[i].equalsIgnoreCase("BestDeltaRT")) {
				bestDeltaRTIdx = i;
			} else if(headers[i].equalsIgnoreCase("Length")) {
				lenIdx = i;
			} else if(headers[i].equalsIgnoreCase("mLog2BestELRank")) {
				baIdx = i;
			}  else if(headers[i].equalsIgnoreCase("AvgDelta")) {
				avgDeltaRTIdx = i;
			}
			
		}
	}
}

public class _FeatComb {

	public static void main(String[] args) throws IOException {
		
		File[] files = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads_features").listFiles();
		for(File file : files) {
			if(file.getName().endsWith("pxg.predfeat.pin")) {
				System.out.println(file.getName());
				genTypeFeat(file);
			}
		}
		

	}
	
	public static void genTypeFeat (File file) throws IOException {

		BufferedReader BR = new BufferedReader(new FileReader(file));
		BufferedWriter BW = new BufferedWriter(new FileWriter(file.getAbsolutePath().replace(".pin", ".feat1.pin")));
		String line = null;
		String header = BR.readLine();
		
		String[] initHeaders = header.split("\t");
		HEADER header_ = new HEADER(initHeaders);
		
		BW.append(header_.toHeader1());
		BW.newLine();
		
		Hashtable<String, String> checkRedundantPeptide = new Hashtable<String, String>();
		ArrayList<String> records = new ArrayList<String>();
		while((line = BR.readLine()) != null) {
			records.add(line);
			String[] fields = line.split("\t");
			String peptide = fields[header_.peptIdx];
			int peptLen = peptide.length();
			Double log2Reads = Double.parseDouble(fields[header_.log2ReadsIdx]);
			
			if(checkRedundantPeptide.get(peptide) == null) {
				checkRedundantPeptide.put(peptide, "");
				if(header_.lengthReads[peptLen] == null) {
					header_.lengthReads[peptLen] = new ArrayList<Double>();
				}
				header_.lengthReads[peptLen].add(log2Reads);
			}
			
		}
		
		for(int i=0; i<header_.lengthReads.length; i++) {
			if(header_.lengthReads[i] == null) {
				continue;
			}
			if(header_.lengthReads[i].size() > 1) {
				for(int j=0; j<header_.lengthReads[i].size(); j++) {
					header_.lengthReadMean[i] += (Double) header_.lengthReads[i].get(j);
				}
				header_.lengthReadMean[i] /= (double)(header_.lengthReads[i].size());
				
				for(int j=0; j<header_.lengthReads[i].size(); j++) {
					header_.lengthReadSigma[i] += Math.pow(((Double) header_.lengthReads[i].get(j)) - header_.lengthReadMean[i],2);
				}
				
				header_.lengthReadSigma[i] /= (double)(header_.lengthReads[i].size()-1);
				header_.lengthReadSigma[i] = Math.sqrt(header_.lengthReadSigma[i]);
			}
		}
		
		
		for(int i=0; i<records.size(); i++) {
			String[] fields = records.get(i).split("\t");
			String record = header_.genRecord1(fields);
			if(record != null)
			{
				BW.append(record);
				BW.newLine();
			}
			
		}
		
		BW.close();
		BR.close();
	}
}
