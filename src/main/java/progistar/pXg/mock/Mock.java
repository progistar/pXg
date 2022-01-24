package progistar.pXg.mock;

import java.util.ArrayList;

import progistar.pXg.constants.Constants;
import progistar.pXg.data.Cigar;
import progistar.pXg.data.GenomicSequence;

public class Mock  {

	/**
	 * Make gSeq with mock property.<br>
	 * 
	 * 
	 * @param gSeq
	 * @param mockMethod
	 * @return
	 */
	public static GenomicSequence makeMockRead(GenomicSequence gSeq, byte mockMethod) {
		assert mockMethod != Constants.MOCK_NONE;
		
		ArrayList<Cigar> cigars = new ArrayList<Cigar>();
		
		StringBuilder revNucleotides = new StringBuilder(gSeq.getNucleotideString());
		revNucleotides = revNucleotides.reverse();
		
		// revCigar
		// containing reverse nucleotides/relativePositions.
		// the whole Cigars are reversed into single Cigar
		
		Cigar revCigar = new Cigar(revNucleotides.length(), 'M');
		revCigar.relativePositions = new int[revNucleotides.length()];
		
		revCigar.nucleotides = revNucleotides.toString();
		
		int index = 0;
		int size = gSeq.cigars.size();
		for(int i=size-1; i>=0; i--) {
			Cigar cigar = gSeq.cigars.get(i);
			// skip zero-size nucleotide
			if(cigar.nucleotides.length() == 0) continue;
			
			int relSize = cigar.relativePositions.length;
			
			for(int j=relSize-1; j>=0; j--) {
				revCigar.relativePositions[index] = cigar.relativePositions[j];
				index++;
			}
			
		}
		cigars.add(revCigar);
		
		GenomicSequence gMSeq = new GenomicSequence("XXX_"+gSeq.uniqueID, gSeq.chrIndex, gSeq.startPosition, cigars, null);
		
		return gMSeq;
	}

}
