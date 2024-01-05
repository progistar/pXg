package progistar.pXg.data;

import java.util.ArrayList;

import progistar.pXg.constants.Constants;
import progistar.pXg.utils.Codon;
import progistar.pXg.utils.IndexConvertor;

public class Output {

	// positions in NGS-read
	// zero-based
	// [startPos, endPos] both are inclusive.
	// forward-oriented position.
	// it means that it is not care about strand.
	public int startPosInNGS;
	public int endPosInNGS;
	
	// if the mapped region resides on junction,
	// there are multiple genomic positions, possible.
	public ArrayList<Integer> startGenomicPositions;
	public ArrayList<Integer> endGenomicPositions;
	
	public int peptideIndex;
	public boolean strand;

	public GenomicSequence gSeq;
	
	public boolean isTarget = true;
	
	public Output (GenomicSequence gSeq, int peptideIndex, int startPos, int endPos, boolean strand, boolean isTarget) {
		this.gSeq = gSeq;
		this.peptideIndex = peptideIndex;
		this.startPosInNGS = startPos;
		this.endPosInNGS = endPos;
		this.strand = strand;
		this.isTarget = isTarget;
	}
	
	public String getPeptide () {
		return PeptideAnnotation.indexedPeptide.get(this.peptideIndex);
	}
	
	public String getMutationStatus () {
		String mutationStatus = "-";
		ArrayList<Mutation> mutations = this.getMutations();
		
		if(mutations.size() != 0) {

			// '-' indicates insertion.
			// reference sequence is "lower case".
			// to translate them, they are converted properly.
			String referenceSequence = this.getMatchedRefNucleotide().replace("-", "").toUpperCase();
			String observedSequence = this.getMatchedNucleotide();
			
			String refPeptide = null;
			String obsPeptide = null;
			
			if(this.strand) {
				refPeptide = GenomicSequence.translation(referenceSequence, 0);
				obsPeptide = GenomicSequence.translation(observedSequence, 0);
			} else {
				refPeptide = GenomicSequence.reverseComplementTranslation(referenceSequence, 0);
				obsPeptide = GenomicSequence.reverseComplementTranslation(observedSequence, 0);
			}
			
			if(refPeptide.equalsIgnoreCase(obsPeptide)) {
				mutationStatus = Constants.MUTATION_SILENT;
			} else {
				mutationStatus = Constants.MUTATION_MISSENSE;
			}
			
			// stop mark
			if(obsPeptide.contains("X")) {
				mutationStatus += "|"+Constants.MUTATION_STOPLOSS;
			}// we cannot determine stop gain. because if it is, then the peptide cannot be detected from MS/MS.
			
			boolean isINDEL = false;
			for(Mutation mutation : mutations) {
				if(mutation.type == Constants.INS || mutation.type == Constants.DEL) {
					isINDEL = true;
				}
			}
			
			if(isINDEL) {
				mutationStatus += "|"+Constants.MUTATION_INDELS;
			}
			
		}
		
		
		return mutationStatus;
	}
	
	public String getMatchedNucleotide () {
		String nucleotide = this.gSeq.getNucleotideString();
		return nucleotide.substring(this.startPosInNGS, this.endPosInNGS+1);
	}
	
	// retrieved reference sequences
	public String getMatchedRefNucleotide () {
		StringBuilder matchedNucleotide = new StringBuilder(this.getMatchedNucleotide());
		StringBuilder refNucleotide = new StringBuilder();
		ArrayList<Mutation> mutations = this.getMutations();
		
		int relPos = 0;
		for(int i=0; i<matchedNucleotide.length(); i++) {
			
			Mutation insOrDelMutation = null;
			for(Mutation mutation : mutations) {
				int mPos = mutation.relPos - this.startPosInNGS;
				
				if(mPos == relPos) {
					if(mutation.type == Constants.SNP) {
						matchedNucleotide.setCharAt(i, mutation.refSeq.toLowerCase().charAt(0));
					} else {
						insOrDelMutation = mutation;
					}
				}
			}
			
			if(insOrDelMutation == null) {
				refNucleotide.append(matchedNucleotide.charAt(i));
			} else {
				if(insOrDelMutation.type == Constants.INS) {
					int insSize = insOrDelMutation.altSeq.length();
					i += insSize - 1;
					for(int j=0; j<insSize; j++) {
						refNucleotide.append("-");
					}
				} else if(insOrDelMutation.type == Constants.DEL) {
					refNucleotide.append(insOrDelMutation.refSeq.toLowerCase());
					refNucleotide.append(matchedNucleotide.charAt(i));
				}
			}
			
			
			relPos++;
		}
		
		return refNucleotide.toString();
	}
	
	public ArrayList<Mutation> getMutations () {
		return this.gSeq.getMutationsByPositionInNGS(this.startPosInNGS, this.endPosInNGS);
	}
	
	/**
	 * return frame annotation such as: <br>
	 * IN_FRAME, NO_FRAME, OUT_OF_FRAME. <br>
	 * 
	 * @param transcriptNum
	 * @return
	 */
	public char getFrame (int transcriptNum) {
		TBlock tBlock = gSeq.tBlocks[transcriptNum];
		// this is intergenic
		// or non-coding
		if(tBlock == null || tBlock.transcriptCodingType == Constants.NON_CODING_TRANSCRIPT) return Constants.NO_FRAME;
		// anti-sense
		if(this.strand != tBlock.strand) return Constants.NO_FRAME;
		
		int size = this.startGenomicPositions.size();

		int genomicSize = 0;
		for(int i=0; i<size; i++) {
			int startPos = this.startGenomicPositions.get(i);
			int endPos = this.endGenomicPositions.get(i);
			genomicSize += endPos - startPos + 1;
		}
		byte[] frames = new byte[genomicSize];
		
		int fIndex = 0;
		for(int i=0; i<size; i++) {
			int startPos = this.startGenomicPositions.get(i);
			int endPos = this.endGenomicPositions.get(i);
			
			for(int gPos=startPos; gPos<=endPos; gPos++) {
				byte frame = tBlock.getFrameMark(gPos);
				frames[fIndex++] = frame;
				
				// If FRAME_X, it means no meaningful information about frame.
				if(frame == Constants.FRAME_X) return Constants.NO_FRAME;
			}
		}
		
		byte targetFrame = Constants.FRAME_0;
		for(int i=0; i<frames.length; i++) {
			if(frames[i] != targetFrame) {
				return Constants.OUT_OF_FRAME;
			} else {
				targetFrame++;
				if(targetFrame > Constants.FRAME_2) targetFrame = Constants.FRAME_0;
			}
		}
		
		// starting from FRAME_0 and ending with FRAME2+1 => FRAME_0.
		// INDELs can ruin this FRAME rule.
		if(targetFrame != Constants.FRAME_0) {
			return Constants.OUT_OF_FRAME;
		}
		
		return Constants.IN_FRAME;
		
	}
	
	/**
	 * return where the mapping is alternative splicing or canonical form <br>
	 * MARK_ALT, MARK_CAN. <br>
	 * 
	 * @param transcriptNum
	 * @return
	 */
	public char getAS (int transcriptNum) {
		TBlock tBlock = gSeq.tBlocks[transcriptNum];
		// this is intergenic
		// in this case, we do not judge where the intergenic is alternative spliced or not
		// just pretending to canonical form
		if(tBlock == null) return Constants.MARK_CA;
		// with soft-clip
		// actually, soft-clip is not allowed.
		if(this.startGenomicPositions.isEmpty() || this.endGenomicPositions.isEmpty()) return Constants.MARK_CA;
		
		// note that
		// genomic size cannot be inferred from peptide length in case of INDELs.
		return tBlock.isAS(this.startGenomicPositions, this.endGenomicPositions);
		
	}
	
	public String getLocus () {
		// unknown locus
		if(this.startGenomicPositions.isEmpty() || this.endGenomicPositions.isEmpty()) {
			return IndexConvertor.indexToChr(gSeq.chrIndex)+":?";
		}
		
		String locus = "";
		for(int i=0; i<startGenomicPositions.size(); i++) {
			if(i!=0) locus += "|";
			locus += IndexConvertor.indexToChr(gSeq.chrIndex) +":" +startGenomicPositions.get(i)+"-"+endGenomicPositions.get(i);
		}
		return locus;
	}
	
	/**
	 * 
	 * Genomic region matched to peptide.<br>
	 * 
	 * @param transcriptNum
	 * @return
	 */
	public String getGenomicRegion (int transcriptNum) {
		StringBuilder genomicRegion = new StringBuilder();
		int relPos = 0;
		for(Cigar cigar : gSeq.cigars) {
			// append sequence
			if(cigar.operation == 'M' || cigar.operation == 'I' || cigar.operation == '*' || cigar.operation == 'S') {
				for(int i=0; i<cigar.annotations.length; i++) {
					if(this.startPosInNGS <= relPos && relPos <= this.endPosInNGS) {
						genomicRegion.append(cigar.annotations[i][transcriptNum]);
					}
					relPos++;
				}
			}
		}
		
		return genomicRegion.toString();
	}
	
	/**
	 * 
	 * get genomic regional information at AA level. <br>
	 * 
	 * @param transcriptNum
	 * @return
	 */
	public String getAARegionAnnotation (int transcriptNum) {
		String genomicRegion = getGenomicRegion(transcriptNum);
		
		StringBuilder aaRegionAnnotation = new StringBuilder();
		for(int i=0; i<genomicRegion.length(); i+=3) {
			char ntMark1 = genomicRegion.charAt(i);
			char ntMark2 = genomicRegion.charAt(i+1);
			char ntMark3 = genomicRegion.charAt(i+2);
			
			aaRegionAnnotation.append(Codon.getAARegion(ntMark1, ntMark2, ntMark3));
		}
		
		
		if(!strand) {
			aaRegionAnnotation = aaRegionAnnotation.reverse();
		}
		// region shrinkage
		StringBuilder shortAnnotation = new StringBuilder();
		int count = 1;
		char mark = aaRegionAnnotation.charAt(0);
		for(int i=1; i<aaRegionAnnotation.length(); i++) {
			char nextMark = aaRegionAnnotation.charAt(i);
			if(mark == nextMark) count++;
			else {
				shortAnnotation.append(count).append(mark);
				mark = nextMark;
				count = 1;
			}
		}
		shortAnnotation.append(count).append(mark);
		
		return shortAnnotation.toString();
		
	}
	/**
	 * Genomic positions of mapped peptides <br>
	 *
	 * 
	 */
	public void mapGenomicAnnotation () {
		// genomic position setting
		this.startGenomicPositions = new ArrayList<Integer>();
		this.endGenomicPositions = new ArrayList<Integer>();
		
		int relPos = 0;
		
		int startGenomicPosition = -1;
		int endGenomicPosition = -1;
		for(Cigar cigar : gSeq.cigars) {
			// if the cigar operation is N,
			// there is a junction in the NGS-read.
			if(cigar.operation == 'N') {
				if(startGenomicPosition != -1 && endGenomicPosition != -1) {
					this.startGenomicPositions.add(startGenomicPosition);
					this.endGenomicPositions.add(endGenomicPosition);
				}
				
				startGenomicPosition = -1;
				endGenomicPosition = -1;
			}
			
			if(cigar.operation == 'M' || cigar.operation == 'I' || cigar.operation == '*' || cigar.operation == 'S') {
				for(int i=0; i<cigar.annotations.length; i++) {
					if(this.startPosInNGS <= relPos && relPos <= this.endPosInNGS) {
						if(startGenomicPosition == -1) startGenomicPosition = cigar.relativePositions[i] + gSeq.startPosition;
						if(endGenomicPosition == -1) endGenomicPosition = cigar.relativePositions[i] + gSeq.startPosition;
						
						startGenomicPosition = Math.min(startGenomicPosition, cigar.relativePositions[i] + gSeq.startPosition);
						endGenomicPosition = Math.max(endGenomicPosition, cigar.relativePositions[i] + gSeq.startPosition);
					}
					relPos++;
				}
			}
		}
		
		if(startGenomicPosition != -1 && endGenomicPosition != -1) {
			this.startGenomicPositions.add(startGenomicPosition);
			this.endGenomicPositions.add(endGenomicPosition);
		}
	}
}

