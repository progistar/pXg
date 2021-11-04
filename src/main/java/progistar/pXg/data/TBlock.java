package progistar.pXg.data;

import java.util.ArrayList;
import java.util.Collections;

import progistar.pXg.constants.Constants;

/**
 * Transcript information block <br>
 * 
 * @author gistar
 *
 */
public class TBlock implements Comparable<TBlock> {
	
	public int tBlockID;
	
	public int chrIndex;
	
	public String transcriptID;
	public String transcriptName;
	public String transcriptType;

	public String geneID;
	public String geneName;
	public String geneType;
	
	public int start;
	public int end;
	public boolean strand;
	
	public byte transcriptCodingType;
	public ArrayList<ABlock> aBlocks = new ArrayList<ABlock>();
	
	public TBlock(int tBlockID, int chrIndex, boolean strand, int start, int end,
			      String transcriptID, String transcriptName, String transcriptType, 
				  String geneID,       String geneName,       String geneType) {
		super();
		this.tBlockID = tBlockID;
		this.chrIndex = chrIndex;
		this.strand = strand;
		this.start = start;
		this.end = end;
		
		this.transcriptID = transcriptID;
		this.transcriptName = transcriptName;
		this.transcriptType = transcriptType;
		this.geneID = geneID;
		this.geneName = geneName;
		this.geneType = geneType;
	}
	
	public void assignBlockTypes () {
		// sort the blocks
		Collections.sort(this.aBlocks);
		
		byte[] blockTypes = new byte[this.end - this.start + 1];
		// the default value of blockTypes, filled with zeros... ( == Constants.INTRON )
		
		int size = this.aBlocks.size();
		boolean isNCDS = true;
		
		for(int i=0; i<size; i++) {
			ABlock aBlock = this.aBlocks.get(i);
			
			for(int idx = aBlock.start; idx <= aBlock.end; idx++) {
				int relIdx = idx - this.start;
				// the first a.Block.feature is CDS or EXON only.
				// CDS is greater score than EXON
				// see Constants Class
				blockTypes[relIdx] = aBlock.feature > blockTypes[relIdx] ? aBlock.feature : blockTypes[relIdx];
				if(blockTypes[relIdx] == Constants.CDS) isNCDS = false;
			}
		}
		
		// if this is non-coding transcript then all exons are interpreted as NCDS (non coding sequence)
		if(isNCDS) {
			// transcript coding type. it simply classifies coding or not.
			transcriptCodingType = Constants.NON_CODING_TRANSCRIPT;
			
			for(int idx=0; idx<blockTypes.length; idx++) {
				if( blockTypes[idx] != Constants.INTRON ) blockTypes[idx] = Constants.NCDS;
			}
			
			// non-coding transcript is assigned "NO_FRAME"
			// there is nothing to do because the default value is NO_FRAME
		} 
		// else if, this is coding transcript. In this case, exons are interpreted as UTR.
		else {
			// transcript coding type. it simply classifies coding or not.
			transcriptCodingType = Constants.CODING_TRANSCRIPT;
			
			boolean isUTR5 = this.strand;
			
			// block type assignment
			for(int idx=0; idx<blockTypes.length; idx++) {
				if( blockTypes[idx] == Constants.EXON && isUTR5 ) blockTypes[idx] = Constants.UTR5;
				else if( blockTypes[idx] == Constants.EXON && !isUTR5 ) blockTypes[idx] = Constants.UTR3;
				else if( blockTypes[idx] == Constants.CDS ) isUTR5 = !this.strand;
			}
		}
		
		// reconstruct aBlocks
		this.aBlocks.clear();
		int start = this.start;
		int end = this.start;
		ABlock aBlock = null;
		
		for(int idx=0; idx<blockTypes.length-1; idx++) {
			// type junction site
			if(blockTypes[idx] != blockTypes[idx+1]) {
				aBlock = new ABlock();
				end = this.start + idx;
				
				aBlock.feature = blockTypes[idx];
				aBlock.transcriptIndex = this.tBlockID;
				aBlock.strand = this.strand;
				aBlock.start = start;
				aBlock.end = end;
				
				start = end + 1;
				
				this.aBlocks.add(aBlock);
			}
		}
		
		// last aBlock.
		aBlock = new ABlock();
		aBlock.feature = blockTypes[blockTypes.length-1];
		aBlock.transcriptIndex = this.tBlockID;
		aBlock.strand = this.strand;
		aBlock.start = start;
		aBlock.end = this.end;
		
		this.aBlocks.add(aBlock);
	}
	
	/**
	 * Take a genomic position as one-based. <br>
	 * 
	 * @param pos
	 * @return
	 */
	public char getRegionMark (int pos) {
		char mark = Constants.MARK_INTERGENIC;
		
		for(ABlock aBlock : this.aBlocks) {
			if(aBlock.start <= pos && aBlock.end >= pos) {
				switch(aBlock.feature) {
				case Constants.CDS :
					mark = Constants.MARK_CDS;
					break;
				case Constants.UTR5 :
					mark = Constants.MARK_5UTR;
					break;
				case Constants.UTR3 :
					mark = Constants.MARK_3UTR;
					break;
				case Constants.INTRON :
					mark = Constants.MARK_INTRON;
					break;
				case Constants.NCDS :
					mark = Constants.MARK_NCDS;
					break;
				default :
					break;
				}
			}
			
			if(mark != Constants.MARK_INTERGENIC) break;
		}
		
		return mark;
	}
	
	/**
	 * currently, if not CDS, return FRAME_X. <br>
	 * If CDS, it returns the frameMark such as: <br>
	 * FRAME_0, FRAME_1, FRAME2. <br>
	 * 
	 * @param pos
	 * @return
	 */
	public byte getFrameMark (int pos) {
		byte mark = Constants.FRAME_X;

		int cds = 0;
		for(ABlock aBlock : this.aBlocks) {
			// inclusive
			if(aBlock.start <= pos && aBlock.end >= pos) {
				if(aBlock.feature == Constants.CDS) {
					cds += (pos - aBlock.start);
					mark = (byte) ( cds % 3);
				}
				
				// return NO_FRAME
				break;
			} 
			// exclusive
			else {
				if(aBlock.feature == Constants.CDS) {
					cds += (aBlock.end - aBlock.start + 1);
				}
			}
		}
		
		return mark;
	}
	
	/**
	 * Taking startLoci and endLoci of a region of mapped peptide. <br>
	 * return MARK_AS or MARK_CA. <br>
	 * 
	 * 
	 * @param startLoci
	 * @param endLoci
	 * @return
	 */
	public char isAS (ArrayList<Integer> startLoci, ArrayList<Integer> endLoci) {
		char isAS = Constants.MARK_AS;
		
		/**
		 * Think!
		 * 
		 * EXON1 - INTRON1 - EXON2 - INTRON2 - EXON3
		 * 
		 * 1) EXON1 - INTRON1 => is consecutive mapping?
		 * 2) EXON1 - EXON2 => is consecutive mapping?
		 * 3) EXON1 - EXON3 => okay. it must be AS
		 * 4) INTRON1 - INTRON2 => okay. it must be AS
		 * 
		 */
		
		// Consecutive mapping must be canonical form (no AS).
		if(startLoci.size() == 1) return Constants.MARK_CA;
		
		// non-consecutive mapping
		// this implies that 
		// 1) EXON1 - INTRON1 => okay. it must be AS
		// 2) EXON1 - EXON2 => both ends are consecutive?
		// 3) EXON1 - EXON3 => okay. it must be AS
		// 4) INTRON1 - INTRON2 => okay. it must be AS
		// So! we just care about case 2. If not, it must be AS
		int size = startLoci.size();
		// storing aBlock number for each start and end loci
		// if the locus resides on intergenic, store 0.
		
		boolean isCA = true;
		
		int direction = 0;
		int prevIdx = 0;
		for(int i=0; i<size; i++) {
			int start = startLoci.get(i);
			int end = endLoci.get(i);
			
			int idx = 0;
			for(ABlock aBlock : this.aBlocks) {
				if(i%2 == 0) {
					if(aBlock.end == end) {
						if(direction != 0) {
							if(idx - prevIdx != 2 || direction == 1) {
								isCA = false;
							}
						}
						
						prevIdx = idx;
						direction = 1;
						break;
					}
				} else {
					if(aBlock.start == start) {
						if(direction != 0) {
							if(idx - prevIdx != 2 || direction == -1) {
								isCA = false;
							}
						}
						
						prevIdx = idx;
						direction = -1;
						break;
					}
				}
				idx++;
			}
			
			// there is no matched regions in the transcript.
			// it implies that one of mapped region fully resides on intergenic region.
			// it cannot be the case of canonical form because exon junctions to intergenic region.
			if(idx == this.aBlocks.size()) {
				isCA = false;
			}
		}
		
		if(isCA) {
			isAS = Constants.MARK_CA;
		}
		
		return isAS;
	}

	public int compareTo(TBlock o) {
		if(this.start < o.start) {
			return -1;
		} else if(this.start > o.start) {
			return 1;
		}
		
		return 0;
	}
	
}
