package progistar.pXg.mock;

import java.util.ArrayList;

import progistar.pXg.constants.Constants;
import progistar.pXg.data.Cigar;
import progistar.pXg.data.GenomicSequence;

/**
 * @deprecated: it is no longer used to estimate target/decoy match
 *
 * @author gistar
 *
 */
public class Mock  {

	private static final int PSD_REV_SIZE = 3;

	/**
	 * @deprecated
	 * Make gSeq with mock property.<br>
	 *
	 *
	 * @param gSeq
	 * @param mockMethod
	 * @return
	 */
	@Deprecated
	public static GenomicSequence makeMockRead(GenomicSequence gSeq, byte mockMethod) {
		assert mockMethod != Constants.MOCK_NONE;

		ArrayList<Cigar> cigars = new ArrayList<>();

		StringBuilder revNucleotides = null;
		String originSequence = gSeq.getNucleotideString();

		if(mockMethod == Constants.MOCK_REVERSE) {
			revNucleotides = new StringBuilder(originSequence);
			revNucleotides = revNucleotides.reverse();
		} else if(mockMethod == Constants.MOCK_PSD_REVERSE) {
			// "Pseudo Reverse"
			// Reverse original sequence codon by codon.
			// ex> GAA TGA GGA CAG GGG => GGG CAG GGA TGA GAA
			int len = originSequence.length();
			revNucleotides = new StringBuilder();
			for(int i=len; i>0; i-=PSD_REV_SIZE) {
				if(i-PSD_REV_SIZE >= 0) {
					revNucleotides.append(originSequence.substring(i-PSD_REV_SIZE, i));
				} else {
					revNucleotides.append(originSequence.substring(0, i));
				}
			}
		}

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
			if(cigar.nucleotides.length() == 0) {
				continue;
			}

			int relSize = cigar.relativePositions.length;

			for(int j=relSize-1; j>=0; j-=PSD_REV_SIZE) {
				if(j-PSD_REV_SIZE >= 0) {
					for(int k=j-PSD_REV_SIZE+1; k<=j; k++) {
						revCigar.relativePositions[index++] = cigar.relativePositions[k];
					}
				} else {
					for(int k=0; k<=j; k++) {
						revCigar.relativePositions[index++] = cigar.relativePositions[k];
					}
				}
			}

		}
		cigars.add(revCigar);

		GenomicSequence gMSeq = new GenomicSequence("XXX_"+gSeq.uniqueID, gSeq.chrIndex, gSeq.startPosition, cigars, null, -1);

		return gMSeq;
	}

}
