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
	public int 		matchedTxds = 0; // in the case of only mapping to intergenic, the number is 1.
	
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
	
	/*
	public String getGenomicInformation () {
		
	}
	*/
	
	public String getFastaHeader () {
		return ">"+this.uniqueID;
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
	 * key: this.uniqueID. <br>
	 * header: <br>
	 * 1) this.uniqueID. <br>
	 * 2) chrName. <br>
	 * TODO genomic position will be multiple..?
	 * 3) genomic positions. <br>
	 * 4) number of contents. <br>
	 * contents: starts with @ <br>
	 * 1) @ENSG_ENST.  <br>
	 * If the locus falls in intergenic region, @Intergenic is annotated <br>
	 * 2) cigar annotations. <br>
	 * 
	 * 
	 * @return
	 */
	public String getGenomicFeature () {
		StringBuilder genomicFeature = new StringBuilder();
		String newLine = System.lineSeparator();
		
		genomicFeature.append(">"+this.uniqueID+"_"+IndexConvertor.indexToChr(this.chrIndex)+":"+this.startPosition+"-"+this.endPosition+"_"+matchedTxds);
		genomicFeature.append(newLine);
		
		for(int i=0; i<matchedTxds; i++) {
			if(tBlocks[i] == null) {
				genomicFeature.append("@Intergenic");
			} else {
				genomicFeature.append("@"+tBlocks[i].geneID+"_"+tBlocks[i].geneName+"_"+tBlocks[i].transcriptID+"_"+tBlocks[i].transcriptType);
			}
			genomicFeature.append(newLine);
			
			for(Cigar cigar : cigars) {
				if(cigar.operation == 'S' || cigar.operation == 'M' || cigar.operation == 'I') {
					for(int j=0; j<cigar.annotations.length; j++) {
						genomicFeature.append(cigar.annotations[j][i]);
					}
				}
			}
			genomicFeature.append(newLine);
		}
		
		return genomicFeature.toString();
	}

	/**
	 * To see what's going on <br>
	 * This is not associated with any functional effects.
	 */
	public void toPrint() {
		
		// TEST Peptide
		StringBuilder nucleotides = new StringBuilder();
		String targetPeptides = "QRFRAGPNM";
		
//		System.out.println(uniqueID+"\t"+startPosition+"-"+endPosition+"\t"+matchedTxds);
		for(Cigar cigar : cigars) {
//			System.out.print(cigar.nucleotides);
			
			// append sequence
			nucleotides.append(cigar.nucleotides);
			
		}
		
		boolean isFind = false;
		for(int i=0; i<3; i++) {
			
			String translatedSequences = translation(nucleotides.toString(), i, false);
			if(translatedSequences.contains(targetPeptides)) {
				isFind = true;
				break;
			}
			translatedSequences = reverseComplementTranslation(nucleotides.toString(), i, false);
			if(translatedSequences.contains(targetPeptides)) {
				isFind = true;
				break;
			}
		}
		
		if(!isFind) return;
		System.out.println(nucleotides.toString());
		System.out.println(uniqueID+"\t"+startPosition+"-"+endPosition+"\t"+matchedTxds);
		
		for(Cigar cigar : cigars) {
			System.out.print(cigar.markerSize+""+cigar.operation);
		}
		System.out.println();
		
		boolean isUTR = false;
		System.out.println(matchedTxds+" is matched Txds");
		for(int i=0; i<matchedTxds; i++) {
			for(Cigar cigar : cigars) {
				if(cigar.operation == 'S' || cigar.operation == 'M' || cigar.operation == 'I') {
					for(int j=0; j<cigar.annotations.length; j++) {
						System.out.print(cigar.annotations[j][i]);
						if(cigar.annotations[j][i] == Constants.MARK_UTR3 || cigar.annotations[j][i] == Constants.MARK_UTR5) isUTR = true;
					}
				}
			}
			
			if(tBlocks[i] == null) {
				System.out.println("\tIntergenic Only");
			} else {
				System.out.println("\t"+tBlocks[i].transcriptID);
			}
		}
		if(isUTR) System.out.println("ABOVE IS UTRs!");
		else System.out.println();
	}
	
	/**
	 * frame is a start position. This is zero-base.
	 * 
	 * @param nucleotides
	 * @param frame
	 * @return
	 */
	public String translation (String nucleotides, int frame, boolean isTerminatedAtStopCodon) {
		StringBuilder peptides = new StringBuilder();
		for(int position=frame; position<nucleotides.length()-2; position+=3) {
			char aa = Codon.nuclToAmino(nucleotides.substring(position,position+3));
			if(isTerminatedAtStopCodon && aa == 'X') break;
			peptides.append(aa);
		}
		return peptides.toString();
	}
	
	public String reverseComplementTranslation (String nucleotides, int frame, boolean isTerminatedAtStopCodon) {
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
			if(isTerminatedAtStopCodon && aa == 'X') break;
			peptides.append(aa);
		}
		return peptides.toString();
	}
	
	public String complementTranslation (String nucleotides, int frame) {
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
		
		nucleotides = reverseComplementNTs.toString();
		for(int position=frame; position<nucleotides.length()-2; position+=3) {
			peptides.append(Codon.nuclToAmino(nucleotides.substring(position,position+3)));
		}
		return peptides.toString();
	}
}
