package progistar.pXg.data;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
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

	/**
	 * Note that about MD tag.
	 * MD tag only presents about SNP and DEL reference sequences.
	 * So, basically Cigar 'M' string is only considered (in case of DEL, we can infer from 'M').
	 * Therefore, MD tag must be resolved by Cigar 'M' only! (do not use other Cigars to resolve MD tag). 
	 * 
	 * 
	 */
	private Pattern EACH_MD_REGEX = Pattern.compile("(([0-9]+)|([A-Z]+|\\^[A-Z]+))");
	
	public String uniqueID;
	public int chrIndex;
	public int startPosition;
	public int endPosition;
	public String mdString;
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
	
	public GenomicSequence (String uniqueID, int chrIndex, int startPosition, ArrayList<Cigar> cigars, String mdStr) {
		this.uniqueID = uniqueID;
		this.chrIndex = chrIndex;
		this.startPosition = startPosition;
		this.cigars = cigars;
		this.mdString = mdStr;
		this.endPosition = startPosition;
		
		for(Cigar cigar : this.cigars) {
			char op = cigar.operation;
			
			switch (op) {
	    	case 'M': case 'I':// match or mismatch or Insertion
	    		this.endPosition = Math.max(this.endPosition, this.startPosition+cigar.relativePositions[cigar.relativePositions.length-1]);
	    		this.endPosition = Math.max(this.endPosition, this.startPosition+cigar.relativePositions[0]);
	    		break;
	    	case '*': // for unmapped
	    		this.startPosition = 1;
	    		this.endPosition = this.startPosition + cigar.markerSize - 1;
	    		break;
    		default :
    			break;
	    	}
		}
	}
	
	/**
	 * Check unmapped status by cigar string.<br>
	 * If there is a cigar with '*', it will return 'false'.<br>
	 * Currently, we do not generate mock reads of unmapped reads. <br>
	 * 
	 * @return
	 */
	public boolean isMapped () {
		assert this.cigars.size() != 0;
		
		if(this.cigars.get(0).operation == '*') {
			return false;
		} else {
			return true;
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
	/**
	 * [start, end] zero-based. <br>
	 * 
	 * @param start
	 * @param end
	 * @return
	 */
	public ArrayList<Mutation> getMutationsByPositionInNGS (int start, int end) {
		ArrayList<Mutation> allMutations = new ArrayList<Mutation>();
		ArrayList<Mutation> inMutations = new ArrayList<Mutation>();
		
		// decoy... has no md string.
		// OR there is no available MD tag: because of alignment option.
		// In the case of STAR2, add --outSAMattributes MD
		if(mdString == null) return inMutations;
		
		// MD parsing
		Matcher mdMatcher = EACH_MD_REGEX.matcher(mdString);
		
		int mRelPos = 0;
		while(mdMatcher.find()) {
			String md = mdMatcher.group();
			char sign = md.charAt(0);
			
			// match size
			if(Character.isDigit(sign)) {
				mRelPos += Integer.parseInt(md);
			} 
			// nt change
			else if(Character.isAlphabetic(sign)) {
				for(int i=0; i<md.length(); i++) {
					mRelPos++;
					Mutation mutation = new Mutation();
					mutation.relPos = mRelPos - 1; // to zero-based
					mutation.refSeq = md.charAt(i)+"";
					allMutations.add(mutation);
				}
			} 
			// deletion sequence
			else if(sign == '^') {
				Mutation mutation = new Mutation();
				mutation.relPos = mRelPos; // to zero-based
				mutation.refSeq = md.substring(1);
				allMutations.add(mutation);
			}
		}
		
		int relPos = 0;
		mRelPos = 0;
		for(Cigar cigar : this.cigars) {
			if(cigar.operation == 'M') {
				for(int i=0; i<cigar.relativePositions.length; i++) {
					if(start <= relPos && relPos <= end) {
						
						for(int j=0; j<allMutations.size(); j++) {
							if(allMutations.get(j).relPos == mRelPos) {
								allMutations.get(j).altSeq = cigar.nucleotides.charAt(i) +"";
								allMutations.get(j).chrIndex = this.chrIndex;
								allMutations.get(j).genomicPosition = this.startPosition + cigar.relativePositions[i];
								allMutations.get(j).type = Constants.SNP;
								inMutations.add(allMutations.get(j));
								
								allMutations.remove(j);
							}
						}
						
					}
					mRelPos++;
					relPos++;
				}
			} else if(cigar.operation == 'I') {
				boolean isIncluded = false;
				for(int i=0; i<cigar.relativePositions.length; i++) {
					if(!isIncluded) {
						if(start <= relPos && relPos <= end) {
							Mutation mutation = new Mutation();
							mutation.altSeq = cigar.nucleotides;
							mutation.chrIndex = this.chrIndex;
							// the relative position of insertion is shifted by + 1 when parsing Cigar.
							// this is tiny issue, so just shift by -1.
							mutation.genomicPosition = this.startPosition + cigar.relativePositions[0] -1;
							mutation.type = Constants.INS;
							inMutations.add(mutation);
							isIncluded = true;
						}
					}
					relPos++;
				}
			} else if(cigar.operation == 'D') {
				for(int i=0; i<cigar.relativePositions.length; i++) {
					if(start <= relPos && relPos <= end) {
						for(int j=0; j<allMutations.size(); j++) {
							if(allMutations.get(j).relPos == mRelPos) {
								allMutations.get(j).chrIndex = this.chrIndex;
								allMutations.get(j).genomicPosition = this.startPosition + cigar.relativePositions[0];
								allMutations.get(j).type = Constants.DEL;
								
								inMutations.add(allMutations.get(j));
								allMutations.remove(j);
							}
						}
					}
				}
			}
		}
		
		return inMutations;
	}
	
	public String getGenomieRegion (int transcriptNum) {
		StringBuilder genomicRegion = new StringBuilder();
		for(Cigar cigar : cigars) {
			// append sequence
			if(cigar.operation == 'M' || cigar.operation == 'I') {
				for(int i=0; i<cigar.annotations.length; i++) {
					genomicRegion.append(cigar.annotations[i][transcriptNum]);
				}
			}
		}
		return genomicRegion.toString();
	}
	
	public String getForwardStrandTranslation (int frame) {
		if(Parameters.leucineIsIsoleucine) {
			return GenomicSequence.translation(getNucleotideString(), frame).replace("I", "L");
		} else {
			return GenomicSequence.translation(getNucleotideString(), frame);
		}
	}
	
	public String getReverseStrandTranslation (int frame) {
		if(Parameters.leucineIsIsoleucine) {
			return GenomicSequence.reverseComplementTranslation(this.getNucleotideString(), frame).replace("I", "L");
		} else {
			return GenomicSequence.reverseComplementTranslation(this.getNucleotideString(), frame);
		}
	}
	/**
	 * frame is a start position. This is zero-base.
	 * 
	 * @param nucleotides
	 * @param frame
	 * @return
	 */
	public static String translation (String nucleotides, int frame) {
		StringBuilder peptides = new StringBuilder();
		for(int position=frame; position<nucleotides.length()-2; position+=3) {
			char aa = Codon.nuclToAmino(nucleotides.substring(position,position+3));
			peptides.append(aa);
		}
		return peptides.toString();
	}
	
	public static String reverseComplementTranslation (String nucleotides, int frame) {
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
