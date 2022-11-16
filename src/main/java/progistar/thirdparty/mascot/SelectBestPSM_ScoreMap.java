package progistar.thirdparty.mascot;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

class MascotRecord_ implements Comparable<MascotRecord_> {
	double score = 0;
	String peptide = null;
	String protein = null;
	String scan = null;
	
	
	@Override
	public int compareTo(MascotRecord_ o) {
		if(this.score < o.score) {
			return 1;
		} else if(this.score > o.score) {
			return -1;
		}
		return 0;
	}
	
	public String toString() {
		StringBuilder str = new StringBuilder();
		str.append(scan+"\t"+protein+"\t"+peptide+"\t"+score);
		return str.toString();
	}
	
}

public class SelectBestPSM_ScoreMap {

	public static void main(String[] args) throws IOException {
		File file = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/LaumontResults/Reanalysis_Mascot_Reformat.tsv");
		BufferedReader BR = new BufferedReader(new FileReader(file));
		BR.readLine();
		String line = null;
		
		int protAccIdx = 1;
		int scoreIdx = 17;
		int sequenceIdx = 20;
		int modIdx = 22;
		int modPosIdx = 23;
		int scanIdx = 24;
		
		// PTM idx
		int deamiIdx = 1;
		int cysteinylIdx = 5;
		int phosphoSTIdx = 3;
		int phosphoYIdx = 4;
		int metIdx = 2;
		String deamiPTM = "Deamidated (NQ)";
		String cysteinylPTM = "Cysteinyl (C)";
		String phosphoSTPTM = "Phospho (ST)";
		String phosphoYPTM = "Phospho (Y)";
		String metPTM = "Oxidation (M)";
		
		
		String proteinStr = null;
		ArrayList<MascotRecord_> records = new ArrayList<MascotRecord_>();
		while((line = BR.readLine()) != null) {
			// skip unassigned spectra
			if(line.startsWith("Unassigned")) {
				continue;
			}
			String[] fields = line.split("\t");
			String thisProteinStr = fields[protAccIdx];
			String peptideSeq = fields[sequenceIdx];
			String scan = fields[scanIdx];
			String score = fields[scoreIdx];
			String modPos = fields[modPosIdx];
			if(modPos.length() != 0) {
				modPos = modPos.split("\\.")[1];
			}
			String mod = "";
			
			for(int i=0; i<modPos.length(); i++) {
				if(modPos.charAt(i)-'0' == deamiIdx) {
					mod += deamiPTM+"-"+peptideSeq.charAt(i)+(i+1)+"|";
				} else if(modPos.charAt(i)-'0' == cysteinylIdx) {
					mod += cysteinylPTM+"-"+peptideSeq.charAt(i)+(i+1)+"|";
				} else if(modPos.charAt(i)-'0' == phosphoSTIdx) {
					mod += phosphoSTPTM+"-"+peptideSeq.charAt(i)+(i+1)+"|";
				} else if(modPos.charAt(i)-'0' == phosphoYIdx) {
					mod += phosphoYPTM+"-"+peptideSeq.charAt(i)+(i+1)+"|";
				} else if(modPos.charAt(i)-'0' == metIdx) {
					mod += metPTM+"-"+peptideSeq.charAt(i)+(i+1)+"|";
				}
			}
			if(mod.length() > 0) {
				mod = mod.substring(0, mod.length()-1);
			}
			
			if(thisProteinStr.length() != 0) {
				proteinStr = thisProteinStr;
			}
			
			MascotRecord_ mr = new MascotRecord_();
			if(mod.length() != 0) {
				mr.peptide = peptideSeq + "|" +mod;
			} else {
				mr.peptide = peptideSeq;
			}
			mr.scan = scan;
			mr.score = Double.parseDouble(score);
			mr.protein = proteinStr;
			records.add(mr);
		}
		
		BR.close();
		
		// aggregate at peptide level
		Hashtable<String, ArrayList<MascotRecord_>> peptideLevel = new Hashtable<String, ArrayList<MascotRecord_>>();
		Hashtable<String, MascotRecord_> spectrumLevel = new Hashtable<String, MascotRecord_>();
		
		for(MascotRecord_ mr : records) {
			String peptide = mr.peptide;
			ArrayList<MascotRecord_> mrs = peptideLevel.get(peptide);
			if(mrs == null) {
				mrs = new ArrayList<MascotRecord_>();
			}
			mrs.add(mr);
			peptideLevel.put(peptide, mrs);
			spectrumLevel.put(mr.scan, mr);
		}
		
		// matching to Laumont report
		File lAllFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/LaumontResults/Laumont_All.tsv");
		int lSequenceIdx = 0;
		int lPTMIdx = 21;
		
		BR = new BufferedReader(new FileReader(lAllFile));
		String lHeader = BR.readLine(); //  skip header
		String[] nullStr = lHeader.split("\t");
		String nullRecord = "";
		lHeader = "";
		for(int i=0; i<nullStr.length; i++) {
			if(i == 0) {
				nullRecord += "0";
				lHeader += "Laumont_"+nullStr[i];
			} else {
				nullRecord += "\t0";	
				lHeader += "\tLaumont_"+nullStr[i];
			}
		}
		final String nullRecordFinal = nullRecord;
		Hashtable<String, String> laumontResMapper = new Hashtable<String, String>();
		while((line = BR.readLine()) != null) {
			line = line.replace("*", "");
			String[] fields = line.split("\t");
			String sequence = fields[lSequenceIdx];
			String mod = fields[lPTMIdx];
			if(!mod.equalsIgnoreCase("none")) {
				sequence += "|" +mod;
			}
			laumontResMapper.put(sequence, line);
		}
		
		BR.close();
		
		File pXgRes = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/LaumontResults/pXg_Subject1.tsv");
		BR = new BufferedReader(new FileReader(pXgRes));
		String pXgHeader = BR.readLine();
		nullStr = pXgHeader.split("\t");
		pXgHeader = "";
		nullRecord = "";
		for(int i=0; i<nullStr.length; i++) {
			if(i == 0) {
				nullRecord += "0";
				pXgHeader += "pXg_"+nullStr[i];
			} else {
				nullRecord += "\t0";	
				pXgHeader += "\tpXg_"+nullStr[i];
			}
		}
		final String nullRecordFinalpXg = nullRecord; 
		
		int pXgFileIdx = 1;
		int pXgScanIdx = 4;
		int pXgChargeIdx = 10;
		int peaksPeptIdx = 3;
		int pXgInferPeptIdx = 20;
		Hashtable<String, String> pXgPeptideMapper = new Hashtable<String, String>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			// spectrum
			//String keyS = fields[pXgFileIdx].split("\\.")[0]+"."+fields[pXgScanIdx].split("\\:")[1]+"."+fields[pXgScanIdx].split("\\:")[1]+"."+fields[pXgChargeIdx];
			// peptide
			String key = fields[pXgInferPeptIdx];
			pXgPeptideMapper.put(key, line);
		}
		
		BR.close();
		
		File pXgResPSMs = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/LaumontResults/pXg_Subject1_PSMs.tsv");
		BR = new BufferedReader(new FileReader(pXgResPSMs));
		BR.readLine();
		Hashtable<String, String> pXgScanMapper = new Hashtable<String, String>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			// spectrum
			String key = fields[pXgFileIdx].split("\\.")[0]+"."+fields[pXgScanIdx].split("\\:")[1]+"."+fields[pXgScanIdx].split("\\:")[1]+"."+fields[pXgChargeIdx];
			// peptide
			pXgScanMapper.put(key, line);
		}
		
		BR.close();
		
		
		File peaksRes = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/9.Unmodified_10ppm/PEAKS/S1.RAW.PEAKS.csv");
		Hashtable<String, String> peaksPeptides = new Hashtable<String, String>();
		BR = new BufferedReader(new FileReader(peaksRes));
		Hashtable<String, String> peaksMapper = new Hashtable<String, String>();
		String peaksHeader = BR.readLine().replace(",", "\t");
		nullStr = peaksHeader.split("\t");
		peaksHeader = "";
		nullRecord = "";
		for(int i=0; i<nullStr.length; i++) {
			if(i == 0) {
				nullRecord += "0";
				peaksHeader += "PEAKS_"+nullStr[i];
			} else {
				nullRecord += "\t0";
				peaksHeader += "\tPEAKS_"+nullStr[i];
			}
		}
		final String nullRecordFilePeaks = nullRecord;
		while((line = BR.readLine()) != null) {
			line = line.replace(",", "\t");
			String[] fields = line.split("\t");
			String key = fields[pXgFileIdx].split("\\.")[0]+"."+fields[pXgScanIdx].split("\\:")[1]+"."+fields[pXgScanIdx].split("\\:")[1]+"."+fields[pXgChargeIdx];
			if(peaksMapper.get(key) == null) {
				peaksMapper.put(key, line);
				peaksPeptides.put(fields[peaksPeptIdx], "");
			}
			
		}
		BR.close();
		
		System.out.println("Mascot_Scan\tMascot_Category\tMascot_Peptide\tMascotScore\t"+lHeader+"\t"+pXgHeader+"\t"+peaksHeader+"\tpXg_Match\tCategory");
		// idx 0: Peptide
		// idx 1: PSM
		int[] overlapBetweenMascotLaumont = new int[2];
		
		peptideLevel.forEach((key, mrs) -> {
			Collections.sort(mrs);
			String laumontReport = laumontResMapper.get(key);
			String scan = mrs.get(0).scan;
			boolean isLaumont = false;
			
			if(laumontReport == null) {
				laumontReport = nullRecordFinal;
			} else {
				isLaumont = true;
				overlapBetweenMascotLaumont[0]++;
				overlapBetweenMascotLaumont[1] += mrs.size();
			}
			
			boolean ispXg = false;
			String pXgReport = pXgPeptideMapper.get(key);
			String pXgScan = null;
			String pXgPept = "";
			if(pXgReport == null) {
				pXgReport = nullRecordFinalpXg;
			} else {
				String[] fields = pXgReport.split("\t");
				pXgScan = fields[pXgFileIdx].split("\\.")[0]+"."+fields[pXgScanIdx].split("\\:")[1]+"."+fields[pXgScanIdx].split("\\:")[1]+"."+fields[pXgChargeIdx];
				ispXg = true;
				pXgPept = fields[pXgInferPeptIdx].replace("I", "L");
			}
			
			
			String category = "?";
			if(ispXg && isLaumont) {
				category = "Overlap";
				//laumontResMapper.remove(key);
				pXgPeptideMapper.remove(key);
			} else if(!ispXg && isLaumont) {
				category = "Laumont";
				if(pXgScanMapper.get(scan) != null) {
					pXgReport = pXgScanMapper.get(scan);
					pXgPept = pXgReport.split("\t")[pXgInferPeptIdx].replace("I", "L");
				}
				//laumontResMapper.remove(key);
			}
			
			String peaksReport = peaksMapper.get(mrs.get(0).scan);
			if(peaksReport == null) {
				peaksReport = nullRecordFilePeaks;
			}
			
			String thisPept = key.split("\\|")[0].replace("I", "L");
			
			
			// is rotated?
			String peaksMark = "NA";
			int[] thisAACount = new int[26];
			int[] peaksAACount = new int[26];
			for(int i=0; i<thisPept.length(); i++) {
				thisAACount[thisPept.charAt(i) - 'A']++;
			}
			for(int i=0; i<pXgPept.length(); i++) {
				peaksAACount[pXgPept.charAt(i) - 'A']++;
			}
			
			boolean isSameAA = true;
			for(int i=0; i<26; i++) {
				if(thisAACount[i] != peaksAACount[i]) {
					isSameAA = false;
				}
			}
			
			if(thisPept.equalsIgnoreCase(pXgPept)) {
				// actually it is exact sequence match,
				// for the sake of simplicity,
				peaksMark = "Equivalent AA composition";
			} else if(isSameAA) {
				peaksMark = "Equivalent AA composition";
			} else {
				peaksMark = "Different AA composition";
			}
			
			if(!category.equalsIgnoreCase("?")) {
				String print = mrs.get(0).toString()+"\t"+laumontReport+"\t"+pXgReport+"\t"+peaksReport+"\t"+peaksMark;
				print = print.replace("PC;sense", "Coding");
				print = print.replace(";sense;", "");
				print = print.replace(";sense", "");
				print = print.replace("sense", "");
				System.out.println(print+"\t"+category);
			}
		});

		pXgPeptideMapper.forEach((key, pXgReport) -> {
			if(laumontResMapper.get(key) == null) {

				if(pXgReport == null) {
					System.out.println("ERROR--;");
					pXgReport = nullRecordFinalpXg;
				}
				
				String[] fields = pXgReport.split("\t");
				String pXgScan = fields[pXgFileIdx].split("\\.")[0]+"."+fields[pXgScanIdx].split("\\:")[1]+"."+fields[pXgScanIdx].split("\\:")[1]+"."+fields[pXgChargeIdx];
				String laumontReport = nullRecordFinal;
				String mascotReport = "0\t0\t0\t0";
				String peaksReport = peaksMapper.get(pXgScan);
				String peaksMark = "NA";
				if(spectrumLevel.get(pXgScan) != null) {
					MascotRecord_ m = spectrumLevel.get(pXgScan);
					mascotReport = m.toString();

					String thisPept = "";
					String pXgPept = fields[pXgInferPeptIdx].replace("I", "L");
					
					if(laumontResMapper.get(m.peptide) != null) {
						laumontReport = laumontResMapper.get(m.peptide);
						thisPept = m.peptide.split("\\|")[0].replace("I", "L");
					}
					
					// is rotated?
					
					int[] thisAACount = new int[26];
					int[] peaksAACount = new int[26];
					for(int i=0; i<thisPept.length(); i++) {
						thisAACount[thisPept.charAt(i) - 'A']++;
					}
					for(int i=0; i<pXgPept.length(); i++) {
						peaksAACount[pXgPept.charAt(i) - 'A']++;
					}
					
					boolean isSameAA = true;
					for(int i=0; i<26; i++) {
						if(thisAACount[i] != peaksAACount[i]) {
							isSameAA = false;
						}
					}
					
					if(thisPept.equalsIgnoreCase(pXgPept)) {
						// actually it is exact sequence match,
						// for the sake of simplicity,
						peaksMark = "Equivalent AA composition";
					} else if(isSameAA) {
						peaksMark = "Equivalent AA composition";
					} else {
						peaksMark = "Different AA composition";
					}
				}
				
				
				String print = mascotReport+"\t"+laumontReport+"\t"+pXgReport+"\t"+peaksReport+"\t"+peaksMark;
				print = print.replace("PC;sense", "Coding");
				print = print.replace(";sense;", "");
				print = print.replace(";sense", "");
				print = print.replace("sense", "");
				System.out.println(print+"\tpXg");
			}
			
			
		});
		
		
		System.out.println("Mascot summary");
		System.out.println("A total of assigned peptides: "+peptideLevel.size());
		System.out.println("A total of assigned PSMs: "+records.size());
		System.out.println("A total of assigned & laumont peptides: "+overlapBetweenMascotLaumont[0]);
		System.out.println("A total of assigned & laumont PSMs: "+overlapBetweenMascotLaumont[1]);
		System.out.println("PEAKS summary");
		System.out.println("A total of peptides: "+peaksPeptides.size());
	}
}
