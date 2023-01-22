package progistar.pXg.train;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

import progistar.pXg.constants.Parameters;
import progistar.pXg.constants.RunInfo;
import progistar.pXg.data.GenomicSequence;
import progistar.pXg.data.parser.SamParser;
import progistar.pXg.data.parser.pXgTSVParser;
import progistar.pXg.utils.Codon;

public class MakePin {

	
	public static void main(String[] args) throws IOException {
		long startTime = System.currentTimeMillis();
		Parameters.leucineIsIsoleucine = false;
		SamParser.ready("/Users/gistar/projects/pXg/Laumont_NatCommun2016/S1.RAW.PEAKS.DecoyOut.pxg.sam");
		
		// for nucleotide sequence trie
		ArrayList<String> records = pXgTSVParser.parseTSV(new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/S1.RAW.PEAKS.DecoyOut.pxg"));
		ArrayList<String> sequences = new ArrayList<String>();
		Hashtable<String, String> checkDuplicatedSequences = new Hashtable<String, String>();
		Hashtable<String, String> leftSeqs = new Hashtable<String, String>();
		Hashtable<String, String> rightSeqs = new Hashtable<String, String>();
		
		for(String record : records) {
			String[] fields = record.split("\t");
			if(fields[23].equalsIgnoreCase("-")) continue;
			if(!fields[37].equalsIgnoreCase("true")) continue;
			
			String peptide = fields[21];
			if(checkDuplicatedSequences.get(fields[24]) == null && checkDuplicatedSequences.get(peptide) == null) {
				sequences.add(fields[24]);
				checkDuplicatedSequences.put(fields[24], "");
				checkDuplicatedSequences.put(peptide, "");
			}
		}
		
		
		// peptide with genomic locus should be used for key in hashtable.
		Trie trie = Trie.builder().addKeywords(sequences).build();
		// loading Codon.
		Codon.mapping();
		
		while(!SamParser.isEndOfFile()) {
			ArrayList<GenomicSequence> genomicSequences = SamParser.parseSam(Parameters.readSize);
			for(GenomicSequence genomicSequence : genomicSequences) {
				String target = genomicSequence.getNucleotideString();
				Collection<Emit> emits = trie.parseText(target);
				
				for(Emit emit : emits) {
					int start = emit.getStart();
					int end = emit.getEnd()+1;
					
					int prefixIdx = Math.max(0, start - 3);
					int suffixIdx = Math.min(target.length(), end + 3);
					
					if(leftSeqs.get(emit.getKeyword()) == null ||
							leftSeqs.get(emit.getKeyword()).length() < target.substring(prefixIdx, start).length()) {
						leftSeqs.put(emit.getKeyword(), target.substring(prefixIdx, start));
					}
					
					if(rightSeqs.get(emit.getKeyword()) == null ||
							rightSeqs.get(emit.getKeyword()).length() < target.substring(end, suffixIdx).length()) {
						rightSeqs.put(emit.getKeyword(), target.substring(end, suffixIdx));
					}
				}
			}
		}
		
		long[][] leftPairs = new long[26][26];
		long[][] rightPairs = new long[26][26];
		for(String sequence : sequences) {
			if(leftSeqs.get(sequence) == null || rightSeqs.get(sequence) == null) {
				
			} else {
				if(leftSeqs.get(sequence).length() == 3 && rightSeqs.get(sequence).length() == 3) {
					String peptide = GenomicSequence.translation(sequence, 0);
					String leftAA = GenomicSequence.translation(leftSeqs.get(sequence), 0);
					String rightAA = GenomicSequence.translation(rightSeqs.get(sequence), 0);
//					System.out.println(leftSeqs.get(sequence)+"\t"+rightSeqs.get(sequence));
//					System.out.println(leftAA+"\t"+rightAA);
					leftPairs[leftAA.charAt(0)-'A'][peptide.charAt(0)-'A']++;
					rightPairs[peptide.charAt(peptide.length()-1)-'A'][rightAA.charAt(0)-'A']++;
				}
			}
		}
		
		System.out.print("AA");
		for(int i=0; i<26; i++) {
			if((char)(i+'A') == 'B' || (char)(i+'A') == 'J' || (char)(i+'A') == 'O' || (char)(i+'A') == 'U' || (char)(i+'A') == 'Z') {
				continue;
			}
			System.out.print("\t"+(char)(i+'A'));
		}
		System.out.println();
		for(int i=0; i<rightPairs.length; i++) {
			if((char)(i+'A') == 'B' || (char)(i+'A') == 'J' || (char)(i+'A') == 'O' || (char)(i+'A') == 'U' || (char)(i+'A') == 'Z') {
				continue;
			}
			System.out.print((char)(i+'A'));
			for(int j=0; j<rightPairs[i].length; j++) {
				
				if((char)(j+'A') == 'B' || (char)(j+'A') == 'J' || (char)(j+'A') == 'O' || (char)(j+'A') == 'U' || (char)(j+'A') == 'Z') {
					continue;
				}
				
				System.out.print("\t"+rightPairs[i][j]);
			}
			System.out.println();
		}
		
		
		SamParser.finish();
    	
    	long endTime = System.currentTimeMillis();
    	System.out.println("\tTotal Elapsed time: "+((endTime-startTime)/1000) + " sec with "+RunInfo.totalProcessedPeptides +" peptides and " +RunInfo.totalProcessedReads+" reads");
	}
}
