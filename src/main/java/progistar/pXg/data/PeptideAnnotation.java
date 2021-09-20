package progistar.pXg.data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

import progistar.pXg.constants.Constants;

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
	
	
	public static void parseTmpResults (File[] files) {
		try {
			
			// IL replaced peptide to PBlocks
			Hashtable<String, ArrayList<PBlock>> pBlockMapper = new Hashtable<String, ArrayList<PBlock>>();
			pBlocks.forEach(pBlock -> {
				String pPeptide = pBlock.getPeptideSequence(); // applied IL-replacement option.
				ArrayList<PBlock> mPBlocks = pBlockMapper.get(pPeptide);
				if(mPBlocks == null) mPBlocks = new ArrayList<PBlock>();
				mPBlocks.add(pBlock);
				
				pBlockMapper.put(pPeptide, mPBlocks);
			});
			
			Hashtable<String, XBlock> xBlocks = new Hashtable<String, XBlock>();
			int maxTargetCount = 0;
			int maxDecoyCount = 0;
			
			for(File file : files) {
				System.out.println("parsing "+file.getName()+" ...");
				BufferedReader BR = new BufferedReader(new FileReader(file));
				String line = null;
				
				String uniqueID = null;
				boolean isDecoy = false;
				
				while((line = BR.readLine()) != null) {
					String[] field = line.split("\t");
					if(field[0].equalsIgnoreCase(Constants.OUTPUT_G_UNIQUE_ID)) {
						uniqueID = field[1];
						// decoy decision
						if(uniqueID.startsWith("XXX")) {
							isDecoy = true;
						} else {
							isDecoy = false;
						}
					} else if(field[0].equalsIgnoreCase(Constants.OUTPUT_G_PEPTIDE)) {
						// skip unknown region
						if(field[2].contains("?")) continue;
						
						String key = field[3] + "_" + field[2];
						
						XBlock xBlock = xBlocks.get(key);
						if(xBlock == null) {
							xBlock = new XBlock();
							xBlock.genomicSequence = field[3];
							xBlock.genomicLocus = field[2];
							xBlocks.put(key, xBlock);
						}
						
						if(isDecoy) {
							xBlock.decoyReadCount++;
							maxDecoyCount = Math.max(maxDecoyCount, xBlock.decoyReadCount);
						} else {
							xBlock.targetReadCount++;
							maxTargetCount = Math.max(maxTargetCount, xBlock.targetReadCount);
						}
					}
				}
				
				
				BR.close();
			}
			
			int[] targetCounts = new int[maxTargetCount+1];
			int[] decoyCounts = new int[maxDecoyCount+1];
			
			System.out.println(maxTargetCount+" and "+maxDecoyCount);
			
			xBlocks.forEach((key, xBlock) -> {
				if(xBlock.targetReadCount != 0) targetCounts[xBlock.targetReadCount]++;
				if(xBlock.decoyReadCount != 0) decoyCounts[xBlock.decoyReadCount]++;
			});
			
			BufferedWriter BW = new BufferedWriter(new FileWriter("readmapping.stat"));
			
			int maxIndex = Math.max(maxTargetCount+1, maxDecoyCount+1);
			
			BW.append("Count\tTarget\tDecoy");
			BW.newLine();
			for(int i=1; i<maxIndex; i++) {
				int targetCount = 0;
				int decoyCount = 0;
				
				if(targetCounts.length > i) targetCount = targetCounts[i];
				if(decoyCounts.length > i) decoyCount = decoyCounts[i];
				
				BW.append(i+"\t"+targetCount+"\t"+decoyCount);
				BW.newLine();
			}
			
			BW.close();
			
		}catch(IOException ioe) {
			
		}
	}
}
