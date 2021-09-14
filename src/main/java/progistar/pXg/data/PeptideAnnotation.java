package progistar.pXg.data;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

import progistar.pXg.decoy.Decoy;

public class PeptideAnnotation {

	
	private static String[]			fields;
	private static Trie				trie;
	public static ArrayList<PBlock>	pBlocks = new ArrayList<PBlock>();
	public static Hashtable<String, Integer>	peptideIndexer	=	new Hashtable<String, Integer>();
	public static Hashtable<Integer, String>	indexedPeptide	=	new Hashtable<Integer, String>();
	
	private PeptideAnnotation() {}
	
	public static void setFields (String[] fields) {
		PeptideAnnotation.fields = fields;
	}
	
	public static void buildKeywordTrie () {
		System.out.println("Enumerate peptide sequences...");
		ArrayList<String> sequences = enumerateSequence();
		System.out.println("A total of "+sequences.size()+" peptides was detected without duplications.");
		System.out.println("Build keyword trie...");
		trie = Trie.builder().addKeywords(sequences).build();
		System.out.println("Done!");
	}
	
	public static ArrayList<Output> find (GenomicSequence gSeq) {
		ArrayList<Output> outputs = new ArrayList<Output>();
		
		for(int i=0; i<3; i++) {
			Collection<Emit> emits = trie.parseText(gSeq.getForwardStrandTranslation(i));
			
			for(Emit emit : emits) {
				// convert peptide index to nucleotide index
				int start = emit.getStart() * 3 + i;
				int end = (emit.getEnd()+1) * 3 + i - 1;
				Output output = new Output(gSeq, peptideIndexer.get(emit.getKeyword()), start, end, true);
				outputs.add(output);
			}
		}
		
		for(int i=0; i<3; i++) {
			int ntLen = gSeq.getNucleotideString().length();
			Collection<Emit> emits = trie.parseText(gSeq.getReverseStrandTranslation(i));
			
			for(Emit emit : emits) {
				// convert peptide index to nucleotide index
				int start = emit.getStart() * 3 + i;
				int end = (emit.getEnd()+1) * 3 + i - 1;
			
				// convert reverse index to forward index
				int tmp = start;
				start = ntLen - end - 1;
				end = ntLen - tmp - 1;
				
				Output output = new Output(gSeq, peptideIndexer.get(emit.getKeyword()), start, end, false);
				outputs.add(output);
			}
		}
		
		
		
		return outputs;
	}
	
	/**
	 * Get non-duplicated peptide sequences from PeptideAnnotation.<br>
	 * 
	 * 
	 * @return
	 */
	private static ArrayList<String> enumerateSequence () {
		ArrayList<String> sequences = new ArrayList<String>();
		
		// put peptide sequences into checks
		pBlocks.forEach(pBlock ->
			{
				if(peptideIndexer.get(pBlock.getPeptideSequence()) == null) {
					indexedPeptide.put(peptideIndexer.size()+1, pBlock.getPeptideSequence());
					peptideIndexer.put(pBlock.getPeptideSequence(), peptideIndexer.size()+1);
					
				}
			}
		);
		
		// add to ArrayList without duplications
		peptideIndexer.forEach((k, v) -> { sequences.add(k); });
		
		return sequences;
	}
	
}
