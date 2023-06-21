package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class Monkey {

	public static void main(String[] args) throws IOException {
		String[] rnas = {
				"GATTTGTATGGGCTTTCCCAGTAGGGA",
				"CAAGGCTTGGACGTCGTTAATAAGTTG",
				"GGCAGCGAGGGTTTTCTCATATGTAAG",
				"GTTCTAGCTCACCAACCGTTTAACTTG",
				"AAGGGTTGCAAGGTTGGGGTACAGTGA",
				"TCTCTTCTTGATGGGGTTAAAACCTTA",
				"CTTCTTGATGGGGTTAAAACCTTAGTC",
				"AGCCATAACTCCCGCGGGAGGTAGAAG",
				"TCTCTTCTTGATGGGGTTAAAACCTTAGTC",
				"GGCCCCTGTATCAAGAATGCCAGAAAAGAGCTT",
				
				"TAGGCTGGGGTTGGAGATATCGAGGGC",
				"TAGGCTGGGGTTGGAGATATTGAGGGC", // mutation 
				"TAGGGTTACAATCAGTTTTAAAGCTAC",
				"CACTCTGCTGATCTTCTCCAAAATATA",
				"ACCCTGTCTGACTACAACATCCGGATA"
		};
		String[] peptide = {
				"SLLGKPIQI",
				"QLINDVQAL",
				"LTYEKTLAA",
				"VLAHQPFNL",
				"SLYPNLATL",
				"SLLDGVKTL",
				"LLDGVKTLV",
				"LLPPAGVMA",
				"SLLDGVKTLV",
				"KLFSGILDTGA",
				
				"ALDISNPSL",
				"ALDISNPSL (N3D)",
				"VALKLIVTL",
				"YILEKISRV",
				"TLSDYNIRI"
		};
		boolean[] strands = {
				false,
				false,
				false,
				true,
				false,
				true,
				true,
				false,
				true,
				false,
				
				false,
				false,
				false,
				false,
				true
		};
		String[] regions = {
				"LTR1","IGR1" ,"gag protein",
				"IGR2", "pol protein", "IGR3", 
				"envelope glycoprotein", "IGR4", "LTR2"
		};
		int[] regionStarts = {
			1, 459, 644, 
			2873, 3891, 6288, 
			6317, 8045, 8417
		};
		
		int[] regionEnds = {
			458, 643, 2872, 
			3890, 6287, 6316, 
			8044, 8416, 8787
		};
		
		File genomeFile = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/SquirrelMonkey/MK561030.1.fa");
		BufferedReader BR = new BufferedReader(new FileReader(genomeFile));
		StringBuilder genome = new StringBuilder();
		BufferedWriter BW = new BufferedWriter(new FileWriter(genomeFile.getAbsolutePath().replace(".fa", ".bed")));
		
		String line = null;
		String header = BR.readLine();// skip header
		while((line = BR.readLine()) != null) {
			genome.append(line);
		}
		BR.close();
		for(int i=0; i<rnas.length; i++) {
			String rna = rnas[i];
			if(!strands[i]) {
				rna = revComp(rna);
			}
			String strand = strands[i] ? "+" : "-";
			int start = genome.indexOf(rna);
			int end = start + rna.length();
			String regionInfo = "";
			for(int j=0; j<regions.length; j++) {
				int rStart = regionStarts[j];
				int rEnd = regionEnds[j];
				if(rStart <= (start+1) && start <= rEnd) {
					if(regionInfo.length() != 0) {
						regionInfo += "&";
					}
					regionInfo += regions[j];
				}
			}
			if(start == -1) {
				System.out.println(peptide[i] +" is not found...");
			} else {
				BW.append("MK561030.1\t"+start+"\t"+end+"\t["+regionInfo+"]"+peptide[i]+"\t100\t"+strand);
				BW.newLine();
			}
		}
		
		BW.close();
	}
	
	public static String revComp (String rna ) {
		String nRNA = "";
		for(int i=0; i<rna.length(); i++) {
			char nt = rna.charAt(i);
			if(nt == 'A') {
				nRNA += "T";
			}else if(nt == 'C') {
				nRNA += "G";
			}else if(nt == 'T') {
				nRNA += "A";
			}else if(nt == 'G') {
				nRNA += "C";
			}
		}
		
		return (new StringBuilder(nRNA).reverse().toString());
	}
}
