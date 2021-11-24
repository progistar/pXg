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
				cigar.relativePositions[index] = cigar.relativePositions[j];
				index++;
			}
		}
		cigars.add(revCigar);
		
		GenomicSequence gDSeq = new GenomicSequence("XXX_"+gSeq.uniqueID, gSeq.chrIndex, gSeq.startPosition, cigars, null);
		
		return gDSeq;
	}
	
	private static String getReverseSequence (String nucleotide) {
		StringBuilder revNucleotides = new StringBuilder(nucleotide);
		return revNucleotides.reverse().toString();
	}
	
	private static String getShuffleSequence (String nucleotide) {
		StringBuilder shfNucleotides = new StringBuilder(nucleotide);
		
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
