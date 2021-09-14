package progistar.pXg.decoy;

import java.util.ArrayList;
import java.util.Random;

import progistar.pXg.constants.Constants;
import progistar.pXg.data.Cigar;
import progistar.pXg.data.GenomicSequence;

public class Decoy  {

	/**
	 * Make gSeq with decoy property.<br>
	 * 
	 * 
	 * @param gSeq
	 * @param decoyMethod
	 * @return
	 */
	public static GenomicSequence makeDecoy(GenomicSequence gSeq, byte decoyMethod) {
		// if decoy is not used, this method MUST BE CALLED by others.
		assert decoyMethod != Constants.DECOY_NONE;
		
		ArrayList<Cigar> cigars = new ArrayList<Cigar>();
		
		for(Cigar cigar : cigars) {
			Cigar revCigar = new Cigar(cigar);
			
			
			// make reverse sequence
			if(decoyMethod == Constants.DECOY_REVERSE) {
				revCigar.nucleotides = getReverseSequence(cigar.nucleotides);
			} else if(decoyMethod == Constants.DECOY_SHUFFLE) {
				revCigar.nucleotides = getShuffleSequence(cigar.nucleotides);
			}
			
			cigars.add(revCigar);
		}
		
		GenomicSequence reverseSeq = new GenomicSequence("XXX_"+gSeq.uniqueID, gSeq.chrIndex, gSeq.startPosition, cigars);
		
		return reverseSeq;
	}
	
	private static String getReverseSequence (String nucleotide) {
		StringBuilder revNucleotides = new StringBuilder();
		revNucleotides.append(nucleotide);
		
		return revNucleotides.reverse().toString();
	}
	
	private static String getShuffleSequence (String nucleotide) {
		StringBuilder shfNucleotides = new StringBuilder();
		shfNucleotides.append(nucleotide);
		
		// Do shuffle
		Random rnd = new Random();
		int nucleotideLength = nucleotide.length();
		
		for(int i=0; i<nucleotideLength; i++) {
			int position = rnd.nextInt(nucleotideLength);
			char n1 = shfNucleotides.charAt(i);
			char n2 = shfNucleotides.charAt(position);
			
			// swap
			shfNucleotides.setCharAt(i, n2);
			shfNucleotides.setCharAt(position, n1);
			
		}
		
		return shfNucleotides.toString();
	}

}
