package progistar.pXg.processor;

import java.util.Hashtable;
import java.util.Iterator;

import progistar.pXg.constants.Constants;
import progistar.pXg.data.Cigar;
import progistar.pXg.data.GenomicAnnotation;
import progistar.pXg.data.GenomicSequence;
import progistar.pXg.data.PeptideAnnotation;
import progistar.pXg.data.TBlock;

public class Mapper {

	/**
	 * genomicSequence must reside on genomicAnnotationIndex! <br>
	 * 
	 * @param genomicSequence
	 * @param start
	 * @param genomicAnnotationIndex
	 * @param genomicAnnotation
	 */
	public static void gMap (GenomicSequence genomicSequence, int start, 
			int[][] genomicAnnotationIndex, GenomicAnnotation genomicAnnotation) {
		int gStartPos = genomicSequence.startPosition;
		int gRelPos = gStartPos - start;
		
		// check if genomicSequence is part of genomicAnnotationIndex.
		assert gRelPos >= 0 && gRelPos < genomicAnnotationIndex.length;
		
		Hashtable<Integer, Boolean> txdIndexer = new Hashtable<Integer, Boolean>();
		
		int[] txdIndices = null;
		// aggregate matched tBlokc IDs
		for(Cigar cigar : genomicSequence.cigars) {
			char op = cigar.operation;
			
			switch (op) {
	    	case 'M': // match or mismatch
	    		for(int i=0; i<cigar.relativePositions.length; i++) {
	    			txdIndices = genomicAnnotationIndex[gRelPos + cigar.relativePositions[i]];
	    			for(int txdIndex : txdIndices) txdIndexer.put(txdIndex, true);
	    		}
	    		break;
	    		
	    	case 'I': // insertion
	    		for(int i=0; i<cigar.relativePositions.length; i++) {
	    			txdIndices = genomicAnnotationIndex[gRelPos + cigar.relativePositions[i]];
	    			for(int txdIndex : txdIndices) txdIndexer.put(txdIndex, true);
	    		}
	    		break;
	    		
	    	case 'D': // deletion
	    		break;	
	    		
	    	case 'N': // skip (ex> exon junction)
	    		break;
	    		
	    	case '*': // unmapped
	    		// NOTHING TO DO
	    		break;
	    		
	    	case 'S': // soft-clip
	    		break;
    		default :
    			break;
	    	}
		}
		
		// matched transcripts...
		// if the size of matched transcripts is zero, then this is intergenic mapping or unmapped
		if(txdIndexer.size() == 0) {
			genomicSequence.setNonTranscripts();
		}
		else {
			// enumerate txd index into tBlocks
			TBlock[] tBlocks = new TBlock[txdIndexer.size()];
			Iterator<Integer> indicies = (Iterator<Integer>) txdIndexer.keys();
			int index = 0;
			while(indicies.hasNext()) {
				int txdIndex = indicies.next();
				tBlocks[index++] = genomicAnnotation.getTBlockByID(txdIndex);
			}
			
			genomicSequence.matchedTxds = tBlocks.length;
			genomicSequence.tBlocks = tBlocks;
			
			for(Cigar cigar : genomicSequence.cigars) {
				char op = cigar.operation;
				
				switch (op) {
		    	case 'M': // match or mismatch
		    		cigar.annotations = new char[cigar.relativePositions.length][tBlocks.length];
		    		
		    		for(int i=0; i<cigar.annotations.length; i++) {
		    			for(int j=0; j<tBlocks.length; j++) {
		    				cigar.annotations[i][j] = tBlocks[j].getRegionMark(cigar.relativePositions[i] + genomicSequence.startPosition);
		    			}
		    		}
		    		
		    		break;
		    		
		    	case 'I': // insertion
		    		cigar.annotations = new char[cigar.relativePositions.length][tBlocks.length];
		    		
		    		for(int i=0; i<cigar.annotations.length; i++) {
		    			for(int j=0; j<tBlocks.length; j++) {
		    				cigar.annotations[i][j] = tBlocks[j].getRegionMark(cigar.relativePositions[i] + genomicSequence.startPosition);
		    			}
		    		}
		    		
		    		break;
		    		
		    	case 'D': // deletion
		    		break;	
		    		
		    	case 'N': // skip (ex> exon junction)
		    		break;
		    		
		    	case 'S': // soft-clip
		    		cigar.annotations = new char[cigar.relativePositions.length][tBlocks.length];
		    		
		    		for(int i=0; i<cigar.annotations.length; i++) {
		    			for(int j=0; j<tBlocks.length; j++) {
		    				cigar.annotations[i][j] = Constants.MARK_SOFTCLIP;
		    			}
		    		}
		    		
		    		break;	
		    	
	    		default :
	    			break;
		    	}
			}
		}
		
	}
}
