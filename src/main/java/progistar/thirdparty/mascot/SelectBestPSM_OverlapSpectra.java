package progistar.thirdparty.mascot;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

class MascotRecord__ implements Comparable<MascotRecord__> {
	double score = 0;
	String peptide = null;
	String protein = null;
	String scan = null;
	
	
	@Override
	public int compareTo(MascotRecord__ o) {
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

public class SelectBestPSM_OverlapSpectra {

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
		ArrayList<MascotRecord__> records = new ArrayList<MascotRecord__>();
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
			
			MascotRecord__ mr = new MascotRecord__();
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
		Hashtable<String, ArrayList<MascotRecord__>> peptideLevel = new Hashtable<String, ArrayList<MascotRecord__>>();
		Hashtable<String, MascotRecord__> spectrumLevel = new Hashtable<String, MascotRecord__>();
		
		for(MascotRecord__ mr : records) {
			String peptide = mr.peptide;
			ArrayList<MascotRecord__> mrs = peptideLevel.get(peptide);
			if(mrs == null) {
				mrs = new ArrayList<MascotRecord__>();
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
			// peptide
			String key = fields[pXgInferPeptIdx];
			pXgPeptideMapper.put(key, line);
		}
		
		BR.close();
		
		File pXgResPSMs = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/LaumontResults/pXg_Subject1_PSMs_Binder.tsv");
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
		System.out.println("Mascot_Scan\tMascot_Category\tMascot_Peptide\tMascot_Score\t"+lHeader+"\t"+pXgHeader+"\tCategory");
		peptideLevel.forEach((key, mrs) -> {
			Collections.sort(mrs);
			String laumontReport = laumontResMapper.get(key);
			String scan = mrs.get(0).scan;
			if(laumontReport != null) {
				
				String pXgReport = pXgScanMapper.get(scan);
				if(pXgReport != null) {

					String print = mrs.get(0).toString()+"\t"+laumontReport+"\t"+pXgReport;
					print = print.replace("PC;sense", "Coding");
					print = print.replace(";sense;", "");
					print = print.replace(";sense", "");
					print = print.replace("sense", "");
					
					String category = "Different sequence";
					if(key.replace("I", "L").equalsIgnoreCase(pXgReport.split("\t")[pXgInferPeptIdx].replace("I", "L"))) {
						category = "Equivalent sequence";
					}
					
					System.out.println(print+"\t"+category);
				}
				
			}
			
			
		});
	}
}
