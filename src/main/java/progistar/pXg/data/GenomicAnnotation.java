package progistar.pXg.data;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;

public class GenomicAnnotation {
	
	// index id (tBlockID) to transcript block
	private Hashtable<Integer, TBlock> idToTx = new Hashtable<Integer, TBlock>();
	// txid to index id
	private Hashtable<String, Integer> txToId = new Hashtable<String, Integer>();
	// transcripts
	private ArrayList<TBlock> tBlocks = new ArrayList<TBlock>();
	
	public void putTBlock 
	(         byte chrIndex, boolean strand, int start, int end,
			  String transcriptID, String transcriptName, String transcriptType, 
			  String geneID,       String geneName,       String geneType) {
		int id = idToTx.size();
		
		if(idToTx.get(id) == null) {
			TBlock tBlock = new TBlock(id, chrIndex, strand, start, end, 
					transcriptID, transcriptName, transcriptType, geneID, geneName, geneType);
			tBlock.transcriptID = transcriptID;
			
			idToTx.put(id, tBlock);
			txToId.put(transcriptID, id);
			tBlocks.add(tBlock);
		}
	}
	
	/**
	 * Retrieve a TBlock matched to ID <br>
	 * 
	 * @param id
	 * @return
	 */
	public TBlock getTBlockByID (int id) {
		return idToTx.get(id);
	}
	
	/**
	 * Retrieve a TBlock matched to transcriptID <br> 
	 * 
	 * @param transcriptID
	 * @return
	 */
	public TBlock getTBlockByTID (String transcriptID) {
		return idToTx.get(txToId.get(transcriptID));
	}
	
	/**
	 * Assign types of annotation blocks in each transcript. <br>
	 * And then, sort transcripts by its start position (increasing order). <br>
	 * 
	 * 
	 */
	public void assignTypesInTBlocks () {
		// assign intra-transcript blocks
		// such as CDS, NCDS, UTR, INTRON
		if(idToTx.size() != 0) {
			idToTx.forEach((k, v) -> {((TBlock)v).assignBlockTypes();});
		}
		
		// sort transcripts by ascending order of its start position.
		Collections.sort(this.tBlocks);
	}
	
	/**
	 *  Use chrIndex defined by IndexConvertor Class. <br>
	 * The start position is 1-based and also inclusive. <br>
	 * 
	 * @param chrIndex
	 * @param start
	 * @param end
	 * @return
	 */
	public int[][] getIndexingBlocks (byte chrIndex, int start, int end) {
		System.out.print("Get indexing tBlocks from the annotation... ("+start+"-"+end+")");
		long startTime = System.currentTimeMillis();
		
		int size = end - start + 1;
		
		int[][] genomicAnnotationIndex 	= new int[size][];
		int[]   annotationIndexSize		= new int[size];
		
		// calculate annotation size of each index
		int sizeOfTBlocks = this.tBlocks.size();
		for(int i=0; i<sizeOfTBlocks; i++) {
			TBlock tBlock = this.tBlocks.get(i);
			
			// transcript start and end locus
			int tStart = tBlock.start;
			int tEnd = tBlock.end;
			
			// chr selection
			if(chrIndex != tBlock.chrIndex) continue;
			// out of range
			if(tStart > end) continue;
			
			for(int j=tStart; j<=tEnd; j++) {
				int relPos = j - start;
				
				if(relPos < 0 || relPos >= size) continue;
				
				annotationIndexSize[relPos] ++;
			}
		}
		
		for(int i=0; i<genomicAnnotationIndex.length; i++) {
			genomicAnnotationIndex[i] = new int[annotationIndexSize[i]];
		}
		
		for(int i=0; i<sizeOfTBlocks; i++) {
			TBlock tBlock = this.tBlocks.get(i);
			
			// transcript start and end locus
			int tStart = tBlock.start;
			int tEnd = tBlock.end;
			
			// chr selection
			if(chrIndex != tBlock.chrIndex) continue;
			// out of range
			if(tStart > end) continue;
			
			for(int j=tStart; j<=tEnd; j++) {
				int relPos = j - start;
				
				if(relPos < 0 || relPos >= size) continue;
				
				genomicAnnotationIndex[relPos][--annotationIndexSize[relPos]] = tBlock.tBlockID;
			}
		}
		
		long endTime = System.currentTimeMillis();
		System.out.println("\tElapsed time: "+((endTime-startTime)/1000) + " sec");
		
		return genomicAnnotationIndex;
	}
}
