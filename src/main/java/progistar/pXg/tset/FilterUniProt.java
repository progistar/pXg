package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

public class FilterUniProt {

	public static void main(String[] args) throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader("C:\\Users\\progi\\Desktop\\Projects\\pXg\\Laumont_NatCommun2016\\Results\\crypticPeptides.txt"));
		String line = null;
		
		ArrayList<String> sequences = new ArrayList<String>();
		ArrayList<String> oriSequences = new ArrayList<String>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[0];
			
			oriSequences.add(peptide);
			sequences.add(peptide.replace("I", "L"));
		}
		
		Trie trie = Trie.builder().addKeywords(sequences).build();
		
		BR.close();
		
		BR = new BufferedReader(new FileReader("C:\\Bioinformatics\\0.Databases\\2.HumanProteins\\uniProt_proteome_iso_202102.fasta"));
		StringBuilder protein = new StringBuilder();
		ArrayList<String> proteins = new ArrayList<String>();
		while((line = BR.readLine()) != null) {
			if(line.startsWith(">")) {
				if(protein.length() != 0) {
					proteins.add(protein.toString());
					protein.setLength(0);
				}
				continue;
			} else {
				protein.append(line.replace("I", "L"));
			}
		}
		
		BR.close();
		
		Hashtable<String, String> finds = new Hashtable<String, String>();
		for(String p : proteins) {
			Collection<Emit> emits = trie.parseText(p);
			
			for(Emit emit : emits) {
				finds.put(emit.getKeyword(), "");
			}
		}
		
		for(int i=0; i<oriSequences.size(); i++) {
			String ilSeq = sequences.get(i);
			if(finds.get(ilSeq) != null) {
				System.out.println(oriSequences.get(i)+"\tUniProt");
			} else {
				System.out.println(oriSequences.get(i)+"\tNovel");
			}
		}
	}
}
