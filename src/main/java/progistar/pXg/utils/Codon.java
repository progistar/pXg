package progistar.pXg.utils;

import java.util.Arrays;

import progistar.pXg.constants.Constants;

public class Codon {
	private static final int nucleoIndexes = 8;
	private static String AminoToNuclArray[][];
	private static char NuclToAminoArray[][][];
	private static char ReversedNuclToAminoArray[][][];
	private static char ReversedComplementNuclToAminoArray[][][];
	private static boolean setOkay = false;
	private static char aminoAcids[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
	
	private static String nucleotides[][] = {
			/*A*/ {"GCT", "GCC", "GCA", "GCG"},
					{},
			/*C*/ {"TGT", "TGC"},
			/*D*/ {"GAT", "GAC"},
			/*E*/ {"GAA", "GAG"},
			/*F*/ {"TTT", "TTC"},
			/*G*/ {"GGT", "GGC", "GGA", "GGG"},
			/*H*/ {"CAT", "CAC"},
			/*I*/ {"ATT", "ATC", "ATA"},
					{},
			/*K*/ {"AAA", "AAG"},
			/*L*/ {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"},
			/*M*/ {"ATG"},
			/*N*/ {"AAT", "AAC"},
					{},
			/*P*/ {"CCT", "CCC", "CCA", "CCG"},
			/*Q*/ {"CAA", "CAG"},
			/*R*/ {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"},
			/*S*/ {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"},
			/*T*/ {"ACT", "ACC", "ACA", "ACG"},
					{},
			/*V*/ {"GTT", "GTA", "GTC", "GTG"},
			/*W*/ {"TGG"},
					{},
			/*Y*/ {"TAT", "TAC"},
					{}
	};
	
	
	public static void mapping() {
		AminoToNuclArray = new String[26][];
		NuclToAminoArray = new char[nucleoIndexes][nucleoIndexes][nucleoIndexes];
		ReversedNuclToAminoArray = new char[nucleoIndexes][nucleoIndexes][nucleoIndexes];
		ReversedComplementNuclToAminoArray = new char[nucleoIndexes][nucleoIndexes][nucleoIndexes];
		
		for(int ntPos = 0; ntPos<nucleoIndexes; ntPos++){
			for(int ntPos_ = 0; ntPos_<nucleoIndexes; ntPos_++){
				Arrays.fill(NuclToAminoArray[ntPos][ntPos_], 'X');
				Arrays.fill(ReversedNuclToAminoArray[ntPos][ntPos_], 'X');
				Arrays.fill(ReversedComplementNuclToAminoArray[ntPos][ntPos_], 'X');
			}
		}
		
		for(char AA : aminoAcids){
			AminoToNuclArray[AA - 'A'] = new String[nucleotides[AA -'A'].length];
			for(int ntPos = 0; ntPos<nucleotides[AA -'A'].length; ntPos++){
				AminoToNuclArray[AA - 'A'][ntPos] = nucleotides[AA -'A'][ntPos];
				NuclToAminoArray
				[nucleotides[AA -'A'][ntPos].charAt(0) & 7]
				[nucleotides[AA -'A'][ntPos].charAt(1) & 7]
				[nucleotides[AA -'A'][ntPos].charAt(2) & 7]
						= AA;
				
				ReversedNuclToAminoArray
				[nucleotides[AA -'A'][ntPos].charAt(2) & 7]
				[nucleotides[AA -'A'][ntPos].charAt(1) & 7]
				[nucleotides[AA -'A'][ntPos].charAt(0) & 7]
						= AA;
				
				char[] RCNts = new char[3];
				for(int nt=0; nt<3; nt++){
					switch(nucleotides[AA -'A'][ntPos].charAt(nt)){
					case 'A': RCNts[2-nt] = 'T'; break;
					case 'C': RCNts[2-nt] = 'G'; break;
					case 'T': RCNts[2-nt] = 'A'; break;
					case 'G': RCNts[2-nt] = 'C'; break;
					}
				}
				
				ReversedComplementNuclToAminoArray
				[RCNts[0] & 7]
				[RCNts[1] & 7]
				[RCNts[2] & 7]
						= AA;
			}
		}
		
		setOkay = true;
	}
	
	
	public static Character nuclToAmino (String nucleotides){
		if(!setOkay) mapping();
		if(nucleotides.length() != 3) return 'X';
		return NuclToAminoArray[nucleotides.charAt(0) & 7][nucleotides.charAt(1) & 7][nucleotides.charAt(2) & 7];
	}
	
	public static String[] aminoToNucl (String amino){
		if(!setOkay) mapping();
		return AminoToNuclArray[amino.charAt(0)-'A'];
	}
	/**
	 * Decide AA-level region annotation.<br>
	 * Note that the worst nucleotide is a representative to AA.<br> 
	 * 
	 * 
	 * @param nt1
	 * @param nt2
	 * @param nt3
	 * @return
	 */
	public static char getAARegion (char nt1, char nt2, char nt3) {
		char aaRegion = Constants.MARK_CDS;
		
		// worst annotation has the highest priority.
		// soft-clip is the worst case... will be discarded!
		if(nt1 == Constants.MARK_UNMAPPED || nt2 == Constants.MARK_UNMAPPED || nt3 == Constants.MARK_UNMAPPED) aaRegion = Constants.MARK_UNMAPPED;
		else if(nt1 == Constants.MARK_SOFTCLIP || nt2 == Constants.MARK_SOFTCLIP || nt3 == Constants.MARK_SOFTCLIP) aaRegion = Constants.MARK_SOFTCLIP;
		else if(nt1 == Constants.MARK_INTERGENIC || nt2 == Constants.MARK_INTERGENIC || nt3 == Constants.MARK_INTERGENIC) aaRegion = Constants.MARK_INTERGENIC;
		else if(nt1 == Constants.MARK_INTRON || nt2 == Constants.MARK_INTRON || nt3 == Constants.MARK_INTRON) aaRegion = Constants.MARK_INTRON;
		else if(nt1 == Constants.MARK_NCDS || nt2 == Constants.MARK_NCDS || nt3 == Constants.MARK_NCDS) aaRegion = Constants.MARK_NCDS;
		else if(nt1 == Constants.MARK_3UTR || nt2 == Constants.MARK_3UTR || nt3 == Constants.MARK_3UTR) aaRegion = Constants.MARK_3UTR;
		else if(nt1 == Constants.MARK_5UTR || nt2 == Constants.MARK_5UTR || nt3 == Constants.MARK_5UTR) aaRegion = Constants.MARK_5UTR;
		
		return aaRegion;
	}
}
