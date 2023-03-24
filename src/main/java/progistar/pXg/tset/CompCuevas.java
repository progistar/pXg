package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

import progistar.thirdparty.netMHCpan.NetMHCpanParser;
import progistar.thirdparty.netMHCpan.NetMHCpanResult;

public class CompCuevas {

	public static void main(String[] args) throws IOException {
//		matchToCuevasResult();
//		matchToReferenceResult();
		
		appendBA();
		System.exit(1);
		
		File refFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/2204_Human_179Contam_UniProteome.fasta");
		File resFile = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/SUDHL4.cuevas.BA.tsv");
		
		
		BufferedReader BR = new BufferedReader(new FileReader(resFile));
		String line = null;
		
		ArrayList<String> peptides = new ArrayList<String>();
		ArrayList<String> results = new ArrayList<String>();
		System.out.println(BR.readLine()+"\tRefSeq"); // print header
		Hashtable<String, Integer> mapper = new Hashtable<String, Integer>();
		
		while((line = BR.readLine()) != null) {
			peptides.add(line.split("\t")[16].replace("I", "L"));
			results.add(line);
		}
		
		BR.close();
		
		Trie trie = Trie.builder().addKeywords(peptides).build();
		BR = new BufferedReader(new FileReader(refFile));
		
		StringBuilder sequence = new StringBuilder();
		while((line = BR.readLine()) != null) {
			if(line.startsWith(">")) {
				if(sequence.length() != 0) {
					String toSequence = sequence.toString().replace("I", "L");
					
					Collection<Emit> matches = trie.parseText(toSequence);
					for(Emit emit : matches) {
						
						Integer cnt = mapper.get(emit.getKeyword());
						if(cnt == null) {
							cnt = 0;
						}
						mapper.put(emit.getKeyword(), ++cnt);
					}
				}
				sequence.setLength(0);
			} else {
				sequence.append(line);
			}
		}
		
		BR.close();
		
		for(int i=0; i<results.size(); i++) {
			Integer cnt = mapper.get(peptides.get(i));
			if(cnt == null) {
				cnt = 0;
			}
			System.out.println(results.get(i)+"\t"+cnt);
		}
		
	}
	
	public static void appendBA () throws IOException {
		File resFile = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/DOHH2.cuevas.BA.tsv");
		File netMHCpan = new File("/Users/gistar/tools/netMHCpan4.1/netMHCpan-4.1/peptide_DOHH2.cuevas.xls");
		NetMHCpanResult netMHCpanResult = NetMHCpanParser.parseNetMHCpan(netMHCpan.getAbsolutePath());
		
		BufferedReader BR = new BufferedReader(new FileReader(resFile));
		String line = null;
		
		String header = BR.readLine();
		System.out.println(header+"\t"+netMHCpanResult.getHeader());
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			System.out.println(line+"\t"+netMHCpanResult.getHLATyping(fields[16]));
		}
		
		BR.close();
	}
	
	public static void matchToCuevasResult () throws IOException {
		File refFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/2204_Human_179Contam_UniProteome.fasta");
		File resFile = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/SUDHL4.cuevas.8-15.tsv");
		File givenCanonicalFile = new File("/Users/gistar/projects/pXg/PreviousStudyResource/Immunopeptidome/SuDHL4/protein_origin_SUDHL4_Canonical.tsv");
		File givenNoncanonicalFile = new File("/Users/gistar/projects/pXg/PreviousStudyResource/Immunopeptidome/SuDHL4/protein_origin_SUDHL4_Non_Canonical.tsv");
		
		BufferedReader BR = new BufferedReader(new FileReader(givenCanonicalFile));
		String line = null;
		BR.readLine(); // skip header
		
		Hashtable<String, String> peptideList = new Hashtable<String, String>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String[] peptides = fields[2].split("\\,");
			
			for(String peptide : peptides) {
				peptideList.put(peptide, "canonical");
			}
		}
		
		BR.close();
		
		BR = new BufferedReader(new FileReader(givenNoncanonicalFile));
		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String[] peptides = fields[2].split("\\,");
			
			for(String peptide : peptides) {
				peptideList.put(peptide, "noncanonical");
			}
		}
		BR.close();
		
		Hashtable<String, Integer> mapper = new Hashtable<String, Integer>();
		ArrayList<String> compPeptides = new ArrayList<String>();
		
		
		BR = new BufferedReader(new FileReader(resFile));
		String header = BR.readLine()+"\tType\tRefSeq";
		
		System.out.println(header);
		Hashtable<String, String> overlapList = new Hashtable<String, String>();
		ArrayList<String> outputs = new ArrayList<String>();
		while((line = BR.readLine()) != null ) {
			String[] fields = line.split("\t");
			String peptide = fields[16];
			String type = peptideList.get(peptide);
			if(type == null) {
				type = "NA";
			} else {
				overlapList.put(peptide, "");
			}
			
			outputs.add(line+"\t"+type);
			compPeptides.add(peptide.replace("I", "L"));
		}
		
		if(overlapList.size() != peptideList.size()) {
			peptideList.forEach((peptide, type)->{
				if(overlapList.get(peptide) == null ) {
					System.out.println(peptide+"\t"+type);
				}
			});
		}
		BR.close();
		
		Trie trie = Trie.builder().addKeywords(compPeptides).build();
		BR = new BufferedReader(new FileReader(refFile));
		
		StringBuilder sequence = new StringBuilder();
		while((line = BR.readLine()) != null) {
			if(line.startsWith(">")) {
				if(sequence.length() != 0) {
					String toSequence = sequence.toString().replace("I", "L");
					
					Collection<Emit> matches = trie.parseText(toSequence);
					for(Emit emit : matches) {
						
						Integer cnt = mapper.get(emit.getKeyword());
						if(cnt == null) {
							cnt = 0;
						}
						mapper.put(emit.getKeyword(), ++cnt);
					}
				}
				sequence.setLength(0);
			} else {
				sequence.append(line);
			}
		}
		
		BR.close();
		
		for(String output : outputs) {
			String[] fields = output.split("\t");
			String peptide = fields[16];
			
			Integer cnt = mapper.get(peptide.replace("I", "L"));
			if(cnt == null) {
				cnt = 0;
			}
			System.out.println(output+"\t"+cnt);
		}
	}
	
	public static void matchToReferenceResult () throws IOException {
		File refFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/2204_Human_179Contam_UniProteome.fasta");
		File givenNoncanonicalFile = new File("/Users/gistar/projects/pXg/PreviousStudyResource/Immunopeptidome/SuDHL4/protein_origin_SUDHL4_Non_Canonical.tsv");
		
		BufferedReader BR = null;
		String line = null;
		ArrayList<String> peptideList = new ArrayList<String>();
		ArrayList<String> records = new ArrayList<String>();
		
		BR = new BufferedReader(new FileReader(givenNoncanonicalFile));
		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String[] peptides = fields[2].split("\\,");
			for(String peptide : peptides) {
				records.add(peptide+"\tNoncanonical");
				peptideList.add(peptide.replace("I", "L"));
			}
		}
		BR.close();
		
		Hashtable<String, Integer> mapper = new Hashtable<String, Integer>();
		
		Trie trie = Trie.builder().addKeywords(peptideList).build();
		BR = new BufferedReader(new FileReader(refFile));
		
		StringBuilder sequence = new StringBuilder();
		while((line = BR.readLine()) != null) {
			if(line.startsWith(">")) {
				if(sequence.length() != 0) {
					String toSequence = sequence.toString().replace("I", "L");
					
					Collection<Emit> matches = trie.parseText(toSequence);
					for(Emit emit : matches) {
						
						Integer cnt = mapper.get(emit.getKeyword());
						if(cnt == null) {
							cnt = 0;
						}
						mapper.put(emit.getKeyword(), ++cnt);
					}
				}
				sequence.setLength(0);
			} else {
				sequence.append(line);
			}
		}
		
		BR.close();
		
		System.out.println("Peptide\tClass\tRefSeq");
		for(int i=0; i<records.size(); i++) {
			String[] fields = records.get(i).split("\t");
			String peptide = fields[0];
			Integer cnt = mapper.get(peptide.replace("I", "L"));
			if(cnt == null) {
				cnt = 0;
			}
			System.out.println(records.get(i)+"\t"+cnt);
		}
	}
}
