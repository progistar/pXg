package progistar.thridparty.statistics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import progistar.pXg.utils.Codon;

public class DatabaseSize {

	public static boolean isProt = false;
	
	public static void main(String[] args) throws IOException {
		//String databaseName = "/Users/gistar/resources/Sequences/gencode.v38.pc_translations.fa";
		String databaseName = "/Users/gistar/resources/Sequences/gencode.v38.transcripts.fa";
		
		File dbFile = new File(databaseName);
		
		BufferedReader BR = new BufferedReader(new FileReader(dbFile));
		
		ArrayList<String> entries = new ArrayList<String>();
		StringBuilder protein = new StringBuilder();
		String line = null;
		
		while((line = BR.readLine()) != null) {
			if(line.startsWith(">")) {
				if(protein.length() != 0) {
					entries.add(protein.toString());
				}
				protein.setLength(0);
			} else {
				protein.append(line);
			}
		}
		
		BR.close();
		
		if(!isProt) {
			ArrayList<String> changedEntries = new ArrayList<String>();
			for(String nt : entries) {
				changedEntries.add(translation(nt, 0));
				changedEntries.add(translation(nt, 1));
				changedEntries.add(translation(nt, 2));
				changedEntries.add(reverseComplementTranslation(nt, 0));
				changedEntries.add(reverseComplementTranslation(nt, 1));
				changedEntries.add(reverseComplementTranslation(nt, 2));
			}
			
			entries = changedEntries;
		}
		
		System.out.println("a total of proteins: "+entries.size());
		// calculate
		Hashtable<String, String> lengthCounts = new Hashtable<String, String>();
		for(String p : entries) {
			p = p.replace("I", "L");
			int length = p.length();
			for(int i=8; i<=length; i++) {
				String peptide = p.substring(i-8, i);
				lengthCounts.put(peptide, "");
			}
		}
		System.out.println("Length 8 : "+lengthCounts.size());
	}
	
	/**
	 * frame is a start position. This is zero-base.
	 * 
	 * @param nucleotides
	 * @param frame
	 * @return
	 */
	public static String translation (String nucleotides, int frame) {
		StringBuilder peptides = new StringBuilder();
		for(int position=frame; position<nucleotides.length()-2; position+=3) {
			char aa = Codon.nuclToAmino(nucleotides.substring(position,position+3));
			peptides.append(aa);
		}
		return peptides.toString();
	}
	
	public static String reverseComplementTranslation (String nucleotides, int frame) {
		StringBuilder peptides = new StringBuilder();
		StringBuilder reverseComplementNTs = new StringBuilder(nucleotides);
		for(int i=0; i<nucleotides.length(); i++) {
			switch(reverseComplementNTs.charAt(i)) {
				case 'A': reverseComplementNTs.setCharAt(i, 'T'); break;
				case 'C': reverseComplementNTs.setCharAt(i, 'G'); break;
				case 'T': reverseComplementNTs.setCharAt(i, 'A'); break;
				case 'G': reverseComplementNTs.setCharAt(i, 'C'); break;
				default : break;
			}
		}
		
		nucleotides = reverseComplementNTs.reverse().toString();
		for(int position=frame; position<nucleotides.length()-2; position+=3) {
			char aa = Codon.nuclToAmino(nucleotides.substring(position,position+3));
			peptides.append(aa);
		}
		return peptides.toString();
	}
}
