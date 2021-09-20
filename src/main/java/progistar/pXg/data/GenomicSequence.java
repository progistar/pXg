package progistar.pXg.data;

import java.util.ArrayList;

import progistar.pXg.constants.Constants;
import progistar.pXg.utils.Codon;
import progistar.pXg.utils.IndexConvertor;

/**
 * Maintain sequence read information: <br>
 * UniqueID: qName <br>
 * chrIndex: index of chr string, which is an auto-increment byte generated from IndexConvertor <br>
 * startPosition: start position of mapped read. 1-based index <br>
 * ArrayList<Cigar>: list of cigars. <br>
 * 
 * @author gistar
 *
 */
public class GenomicSequence {

	public String uniqueID;
	public byte chrIndex;
	public int startPosition;
	public int endPosition;
	public ArrayList<Cigar> cigars;
	
	
	// annotation
	// if there is no model on the mapped region, then it means this read maps on intergenic region only.
	// => in this case, tBlocks has only one element with null.
	// otherwise, there is at least one annotated model.
	// => in this case, tBlocks has annotated transcript models and fill with their object.
	// The assignment, see Mapper Class
	public TBlock[]	tBlocks; // we can get strand from tBlock!
	public int 		matchedTxds = 1; // in the case of only mapping to intergenic, the number is 1.
	
	// default is three frame translation.
	// it has a specific translation frame as if it can infer from annotation.
	/**
	 * 
	 * Deprecated:: Assume that known peptides are filtered in advance. (known peptides are not our interest)
	 * 
	 * Frame decision problem.
	 * If we can find in-frame, then we should give an option of in-frame translation.
	 * To do so, we MUST decide when/how/what in-frame. 
	 * 
	 * Non-coding: three-frame
	 * == if the transcript is non-coding such as pseudogene, lncRNA something like that.
	 * 
	 * Coding: three-frame or in-frame
	 * == in-frame: if we can reasonably infer in-frame...
	 * ==== Of course, at least, the transcript is coding.
	 * ==== Assume that translation gets started from up-stream.
	 * ==== However, read can be aligned on UTRs or Intergenic as the first region.
	 * ==== In above case, frame follows the CDS closest to up-stream.
	 * 
	 */
	//public byte[] 	transFrames; // byte[sizeOfTranscripts]
	// In case of intergenic (it implies that this sequence cannot be explained by annotations), transFrames = new byte[1].
	
	public GenomicSequence (String uniqueID, byte chrIndex, int startPosition, ArrayList<Cigar> cigars) {
		this.uniqueID = uniqueID;
		this.chrIndex = chrIndex;
		this.startPosition = startPosition;
		this.cigars = cigars;
		
		this.endPosition = startPosition;
		for(Cigar cigar : this.cigars) {
			char op = cigar.operation;
			
			switch (op) {
	    	case 'M': case 'I':// match or mismatch or Insertion
	    		this.endPosition = Math.max(this.endPosition, this.endPosition+cigar.relativePositions[cigar.relativePositions.length-1]);
	    		break;
    		default :
    			break;
	    	}
		}
	}
	public String getLocus () {
		return IndexConvertor.indexToChr(chrIndex) +":" +this.startPosition+"-"+this.endPosition;
	}
	
	public String getNucleotideString () {
		StringBuilder nucleotides = new StringBuilder();
		for(Cigar cigar : cigars) {
			// append sequence
			nucleotides.append(cigar.nucleotides);
		}
		return nucleotides.toString();
	}
	
	public String getGenomieRegion (int transcriptNum) {
		StringBuilder genomicRegion = new StringBuilder();
		for(Cigar cigar : cigars) {
			// append sequence
			if(cigar.operation == 'S' || cigar.operation == 'M' || cigar.operation == 'I') {
				for(int i=0; i<cigar.annotations.length; i++) {
					genomicRegion.append(cigar.annotations[i][transcriptNum]);
				}
			}
		}
		return genomicRegion.toString();
	}
	
	public String getForwardStrandTranslation (int frame) {
		return this.translation(getNucleotideString(), frame);
	}
	
	public String getReverseStrandTranslation (int frame) {
		return this.reverseComplementTranslation(this.getNucleotideString(), frame);
	}
	/**
	 * frame is a start position. This is zero-base.
	 * 
	 * @param nucleotides
	 * @param frame
	 * @return
	 */
	public String translation (String nucleotides, int frame) {
		StringBuilder peptides = new StringBuilder();
		for(int position=frame; position<nucleotides.length()-2; position+=3) {
			char aa = Codon.nuclToAmino(nucleotides.substring(position,position+3));
			peptides.append(aa);
		}
		return peptides.toString();
	}
	
	public String reverseComplementTranslation (String nucleotides, int frame) {
		StringBuilder peptides = new StringBuilder();
		StringBuilder reverseComplementNTs = new StringBuilder(nucleotides);
		for(int i=0; i<nucleotides.length(); i++) {
			switch(reverseComplementNTs.charAt(i)) {
				case 'A': reverseComplementNTs.setCharAt(i, 'T'); break;
				case 'C': reverseComplementNTs.setCharAt(i, 'G'); break;
				case 'T': reverseComplementNTs.setCharAt(i, 'A'); break;
				case 'G': reverseComplementNTs.setCharAt(i, 'C'); break;
				default : break;
			}
		}
		
		nucleotides = reverseComplementNTs.reverse().toString();
		for(int position=frame; position<nucleotides.length()-2; position+=3) {
			char aa = Codon.nuclToAmino(nucleotides.substring(position,position+3));
			peptides.append(aa);
		}
		return peptides.toString();
	}
}
