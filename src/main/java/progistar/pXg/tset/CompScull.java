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

public class CompScull {

	public static void main(String[] args) throws IOException {
//		extractSequence();
		referenceAppender();
	}
	
	public static void extractSequence () throws IOException {
		File file = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/THP1_3.scull.tsv");
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		String header = BR.readLine();
		System.out.println(header+"\tSequence");
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String sequence = fields[0].replaceAll("[+-0123456789*\\(\\)]", "");
			System.out.println(line+"\t"+sequence);
		}
		
		BR.close();
	}
	
	public static void referenceAppender() throws IOException {

		File refFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/2204_Human_179Contam_UniProteome.fasta");
		File resFile = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/THP1_3.scull.tsv");
		
		
		BufferedReader BR = new BufferedReader(new FileReader(resFile));
		String line = null;
		
		ArrayList<String> peptides = new ArrayList<String>();
		ArrayList<String> results = new ArrayList<String>();
		System.out.println(BR.readLine()+"\tRefSeq"); // print header
		Hashtable<String, Integer> mapper = new Hashtable<String, Integer>();
		
		while((line = BR.readLine()) != null) {
			peptides.add(line.split("\t")[17].replace("I", "L"));
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
}
