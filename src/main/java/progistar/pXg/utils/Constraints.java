package progistar.pXg.utils;

import java.util.ArrayList;

import progistar.pXg.data.GenomicSequence;

public class Constraints {

	/**
	 * Check if GenomicSequences are sorted by descending order.<br>
	 * 
	 * @param gSeqs
	 * @return
	 */
	public static boolean isSortedSAM (ArrayList<GenomicSequence> gSeqs) {
		boolean isSorted = true;
		int size = gSeqs.size();
		
		for(int i=1; i<size; i++) {
			GenomicSequence pGSeq = gSeqs.get(i-1);
			GenomicSequence nGSeq = gSeqs.get(i);
			
			// assume that same chromosomes are clustered at least..!
			if(pGSeq.chrIndex == nGSeq.chrIndex) {
				if(pGSeq.startPosition > nGSeq.startPosition) isSorted = false;
			}
		}
		
		return isSorted;
	}
}
