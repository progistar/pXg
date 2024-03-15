package progistar.pXg.data;

import java.util.ArrayList;
import java.util.Hashtable;

import javax.swing.plaf.basic.BasicScrollBarUI;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.data.parser.ResultParser;
import progistar.pXg.utils.ENSTMapper;
import progistar.pXg.utils.Priority;

public class XBlock {
	// Note that only unmapped reads store sequenceID.
	// Mapped reads do not need this information
	// ==> ** Fixed: All reads store sequenceID. ** <==
	public String sequenceID		=	null;
	
	public int targetReadCount		=	0;
	public int mockReadCount		=	0;
	public char strand				=	'+';
	public String genomicLocus		=	null;
	public String mutations 		=	null;
	
	public int genomicSequenceIdx	=	-1;
	public int leftFlankSequenceIdx	=	-1;
	public int rightFlankSequenceIdx =	-1;
	
	public int referenceSequenceIdx	=	-1;
	public int leftFlankRefSequenceIdx	=	-1;
	public int rightFlankRefSequenceIdx	=	-1;
	
	public String mutationStatus	=	null;
	public String peptideSequence	=	null;
	public String tAnnotations		=	null; // transcript and additional annotations
	public String[] fastaIDs		=	new String[0];
	public String fullReadSequence	=	null; // for unmapped read
	public double bestRegionPriority 	= 	Double.MAX_VALUE;
	public double targetQScore		=	0;
	public double decoyQScore		=	0;
	
	public String[] exonLenghts = null;
	public String[] percentFullDistances	=	null;
	public String[] percentExonDistances	=	null;
	public String[] percentCDSDistances	=	null;
	public String[] fromStartDistances	=	null;
	public String[] fromStopDistances	=	null;
	
	// with the same key value block
	// the block will be null if there is no next sibling.
	public ArrayList<XBlock> siblingXBlocks	= new ArrayList<XBlock>();
	
	/**
	 * Jan. 26, 2024
	 * Get consensus XBlock based on left and right flank sequences.
	 * 
	 * @return
	 */
	public XBlock getConsensusSequenceXBlock () {
		XBlock consensusXBlock = this;
		ArrayList<XBlock> list = new ArrayList<XBlock>();
		list.add(this);
		list.addAll(this.siblingXBlocks);
		
		int maxLeftSize = 0;
		int maxRightSize = 0;
		for(XBlock xBlock : list) {
			String leftFlankSequence = Global.SEQUENCE_ARRAYLIST.get(xBlock.leftFlankRefSequenceIdx);
			String rightFlankSequence = Global.SEQUENCE_ARRAYLIST.get(xBlock.rightFlankSequenceIdx);
			
			maxLeftSize = Math.max(leftFlankSequence.length(), maxLeftSize);
			maxRightSize = Math.max(rightFlankSequence.length(), maxRightSize);
		}
		
		if(maxLeftSize > Parameters.maxFlankNSize) {
			maxLeftSize = Parameters.maxFlankNSize;
		}
		if(maxRightSize > Parameters.maxFlankNSize) {
			maxRightSize = Parameters.maxFlankNSize;
		}
		
		int[][] leftConsensusScore = new int[maxLeftSize][4];
		int[][] rightConsensusScore = new int[maxRightSize][4];
		
		
		// calculate score matrix
		for(XBlock xBlock : list) {
			// left
			String sequence = Global.SEQUENCE_ARRAYLIST.get(xBlock.leftFlankRefSequenceIdx);
			int idx = maxLeftSize-1;
			for(int i=sequence.length()-1; i>=0; i--) {
				char nt = sequence.charAt(i);
				int ntIdx = -1;
				if(nt == 'A') ntIdx = 0;
				else if(nt == 'C') ntIdx = 1;
				else if(nt == 'T') ntIdx = 2;
				else if(nt == 'G') ntIdx = 3;
				if(ntIdx != -1) {
					leftConsensusScore[idx][ntIdx] ++;
				}
				idx--;
			}
			
			// right
			sequence = Global.SEQUENCE_ARRAYLIST.get(xBlock.rightFlankSequenceIdx);
			for(int i=0; i<sequence.length(); i++) {
				char nt = sequence.charAt(i);
				int ntIdx = -1;
				if(nt == 'A') ntIdx = 0;
				else if(nt == 'C') ntIdx = 1;
				else if(nt == 'T') ntIdx = 2;
				else if(nt == 'G') ntIdx = 3;
				if(ntIdx != -1) {
					rightConsensusScore[i][ntIdx] ++;
				}
			}
		}
		
		int bestScore = 0;
		for(XBlock xBlock : list) {
			int score = 0;
			// left
			String sequence = Global.SEQUENCE_ARRAYLIST.get(xBlock.leftFlankRefSequenceIdx);
			int idx = maxLeftSize-1;
			for(int i=sequence.length()-1; i>=0; i--) {
				char nt = sequence.charAt(i);
				int ntIdx = -1;
				if(nt == 'A') ntIdx = 0;
				else if(nt == 'C') ntIdx = 1;
				else if(nt == 'T') ntIdx = 2;
				else if(nt == 'G') ntIdx = 3;
				if(ntIdx != -1) {
					score += leftConsensusScore[idx][ntIdx];
				}
				idx--;
			}
			
			// right
			sequence = Global.SEQUENCE_ARRAYLIST.get(xBlock.rightFlankSequenceIdx);
			for(int i=0; i<sequence.length(); i++) {
				char nt = sequence.charAt(i);
				int ntIdx = -1;
				if(nt == 'A') ntIdx = 0;
				else if(nt == 'C') ntIdx = 1;
				else if(nt == 'T') ntIdx = 2;
				else if(nt == 'G') ntIdx = 3;
				if(ntIdx != -1) {
					score += rightConsensusScore[i][ntIdx];
				}
			}
			
			if(score > bestScore) {
				consensusXBlock = xBlock;
				bestScore = score;
			}
		}
		
		return consensusXBlock;
	}
	
	public boolean isCannonical () {
		String events = toEvents().get("key");
		boolean isCannonical = false;
		
		// wildtype
		if(events.equalsIgnoreCase(Constants.EVENT_PROTEINCODING) && 
				!mutationStatus.equalsIgnoreCase(Constants.MUTATION_ALTERED)) {
			// PC
			isCannonical = true;
		}
		
		return isCannonical;
	}
	
	public boolean isMappedAmbiguous () {
		Hashtable<String, String> events = toEvents();
		String event = events.get("key");
		if(event.contains(Constants.EVENT_UNKNOWN)) {
			return false;
		}
		return true;
	}
	
	public String getKey () {
		String genomicSequence = Global.SEQUENCE_ARRAYLIST.get(this.genomicSequenceIdx);
		return genomicSequence+"_"+this.genomicLocus;
	}
	
	/**
	 * return record. <br>
	 * 
	 */
	public String toString (byte psmStatus) {
		XBlock consensusXBlock = this.getConsensusSequenceXBlock();
		Hashtable<String, String> geneIDs = toGeneIDs();
		Hashtable<String, String> geneNames = toGeneNames();
		Hashtable<String, String> events = toEvents();
		Hashtable<String, String> fastaIDs = toFastaIDs();
		
		int transCount = 0;
		String[] transcripts = tAnnotations.split("\\|");
		for(String transcript : transcripts) {
			if(!transcript.startsWith(Constants.EVENT_INTERGENIC) && 
					!transcript.startsWith(Constants.EVENT_UNKNOWN)) {
				transCount++;
			}
		}
		
		String lfs = Global.SEQUENCE_ARRAYLIST.get(consensusXBlock.leftFlankSequenceIdx).length() == 0 ? 
				"-" : Global.SEQUENCE_ARRAYLIST.get(consensusXBlock.leftFlankSequenceIdx);
		
		String rfs = Global.SEQUENCE_ARRAYLIST.get(consensusXBlock.rightFlankSequenceIdx).length() == 0 ? 
				"-" : Global.SEQUENCE_ARRAYLIST.get(consensusXBlock.rightFlankSequenceIdx);
		
		String lfsRef = Global.SEQUENCE_ARRAYLIST.get(consensusXBlock.leftFlankRefSequenceIdx).length() == 0 ? 
				"-" : Global.SEQUENCE_ARRAYLIST.get(consensusXBlock.leftFlankRefSequenceIdx);
		
		String rfsRef = Global.SEQUENCE_ARRAYLIST.get(consensusXBlock.rightFlankRefSequenceIdx).length() == 0 ? 
				"-" : Global.SEQUENCE_ARRAYLIST.get(consensusXBlock.rightFlankRefSequenceIdx);
		
		String genomicSequence = Global.SEQUENCE_ARRAYLIST.get(this.genomicSequenceIdx);
		String referenceSequence = Global.SEQUENCE_ARRAYLIST.get(this.referenceSequenceIdx);
		
		if(psmStatus == Constants.PSM_STATUS_DECOY) {
			
			double qScore = this.decoyQScore;
			if(Parameters.PHRED_CAL.equalsIgnoreCase(Constants.CAL_PHRED_AVG)) {
				qScore /= (double)mockReadCount;
			}
			
			return new StringBuilder(peptideSequence).reverse().toString() +"\t"+genomicLocus+"\t"
					+strand
					+"\t"+lfs
					+"\t"+genomicSequence
					+"\t"+rfs
					+"\t"+lfsRef
					+"\t"+referenceSequence
					+"\t"+rfsRef
					+"\t"+mutations
					+"\t"+mutationStatus
					+"\t"+tAnnotations+"\t"+transCount
					+"\t"+geneIDs.get("key")+"\t"+geneIDs.get("count")
					+"\t"+geneNames.get("key")+"\t"+geneNames.get("count")
					+"\t"+this.toDist()
					+"\t"+events.get("key")+"\t"+events.get("count")
					+"\t"+fastaIDs.get("key")+"\t"+fastaIDs.get("count")
					+"\t"+mockReadCount+"\t"+(qScore);
		} else {
			
			double qScore = this.targetQScore;
			if(Parameters.PHRED_CAL.equalsIgnoreCase(Constants.CAL_PHRED_AVG)) {
				qScore /= (double)targetReadCount;
			}
			
			return peptideSequence +"\t"+genomicLocus+"\t"
					+strand
					+"\t"+lfs
					+"\t"+genomicSequence
					+"\t"+rfs
					+"\t"+lfsRef
					+"\t"+referenceSequence
					+"\t"+rfsRef
					+"\t"+mutations
					+"\t"+mutationStatus
					+"\t"+tAnnotations+"\t"+transCount
					+"\t"+geneIDs.get("key")+"\t"+geneIDs.get("count")
					+"\t"+geneNames.get("key")+"\t"+geneNames.get("count")
					+"\t"+this.toDist()
					+"\t"+events.get("key")+"\t"+events.get("count")
					+"\t"+fastaIDs.get("key")+"\t"+fastaIDs.get("count")
					+"\t"+targetReadCount+"\t"+(qScore);
		}
	}
	
	/**
	 * Check whether there is matched fasta entries.<br>
	 * 
	 * 
	 * @return
	 */
	public boolean isFastaAssigned () {
		if(fastaIDs == null || this.fastaIDs.length == 0) {
			return false;
		}
		return true;
	}
	
	private String toDist () {
		int longestIdx = 0;
		int longestVal = 0;
		for(int i=0; i<this.exonLenghts.length; i++) {
			if(this.exonLenghts[i].equalsIgnoreCase("-")) continue;
			
			int exonLength = Integer.parseInt(this.exonLenghts[i]);
			if(exonLength > longestVal) {
				longestIdx = i;
				longestVal = exonLength;
			}
		}
		
		return percentFullDistances[longestIdx] + "\t" + 
		percentExonDistances[longestIdx] + "\t" +
		percentCDSDistances[longestIdx] + "\t" +
		fromStartDistances[longestIdx] + "\t" +
		fromStopDistances[longestIdx];
	}
	
	/**
	 * @deprecated
	 * Do not print PSM without read mapping.<br>
	 * -- The original usage is to print out PSMs with no reads.<br>
	 * @param fastaIDs
	 * @return
	 */
	private static String toFastaIDs (String[] fastaIDs) {
		if(fastaIDs == null || fastaIDs.length == 0) {
			return Constants.ID_NULL;
		}
		
		StringBuilder fasta= new StringBuilder();
		Hashtable<String, Boolean> isDuplicated = new Hashtable<String, Boolean>();
		
		for(String fastaID : fastaIDs) {
			if(isDuplicated.get(fastaID) == null) {
				fasta.append("|").append(fastaID);
				isDuplicated.put(fastaID, true);
			}
		}
		
		return fasta.substring(1).toString();
	}
	
	/**
	 * The return value is hashtable consisting of "key" and "count"<br>
	 * "key" represents the annotation of mapped IDs. <br>
	 * "count" represents the number of IDs. <br>
	 * 
	 * @return
	 */
	public Hashtable<String, String> toFastaIDs () {
		Hashtable<String, String> mapper = new Hashtable<String, String>();
		
		if(isFastaAssigned()) {
			StringBuilder fasta= new StringBuilder();
			Hashtable<String, Boolean> isDuplicated = new Hashtable<String, Boolean>();
			
			for(String fastaID : this.fastaIDs) {
				if(isDuplicated.get(fastaID) == null) {
					fasta.append("|").append(fastaID);
					isDuplicated.put(fastaID, true);
				}
			}
			String[] annotations = fasta.substring(1).toString().split("\\|");
			mapper.put("count", annotations.length+"");
			
			StringBuilder annotation = new StringBuilder();
			for(int i=0; i<annotations.length; i++) {
				if( i == Parameters.maxProteinOut) {
					annotation.append("|").append("...");
					break;
				}
				annotation.append("|").append(annotations[i]);
			}
			mapper.put("key", annotation.toString().substring(1));
		} else {
			// default behavior when there is no IDs.
			mapper.put("key", Constants.ID_NULL);
			mapper.put("count", "0");
		}
		
		
		return mapper;
	}
	
	/**
	 * The return value is hashtable consisting of "key" and "count"<br>
	 * "key" represents the annotation of matched events. <br>
	 * "count" represents the number of events. <br>
	 * 
	 * @param transcriptIDs
	 * @return
	 */
	public Hashtable<String, String> toEvents () {
		Hashtable<String, String> mapper = new Hashtable<String, String>();
		StringBuilder events= new StringBuilder();
		String[] genes = tAnnotations.split("\\|");
		
		Hashtable<String, Boolean> isDuplicated = new Hashtable<String, Boolean>();
		for(String gene : genes) {
			String event = Priority.getRegionEvent(gene);
			if(isDuplicated.get(event) == null) {
				events.append("|").append(event);
				isDuplicated.put(event, true);
			}
		}
		
		// there is no matched events? it is impossible.
		// if then, it must be treated as an error.
		assert events.length() != 0;
		String annotation = events.substring(1).toString();
		mapper.put("key", annotation);
		mapper.put("count", annotation.split("\\|").length+"");
			
		return mapper;
	}
	
	/**
	 * The return value is hashtable consisting of "key" and "count"<br>
	 * "key" represents the annotation of matched gene names. <br>
	 * "count" represents the number of gene names. <br>
	 * 
	 * @return
	 */
	private Hashtable<String, String> toGeneNames () {
		Hashtable<String, String> mapper = new Hashtable<String, String>();
		StringBuilder gRegions = new StringBuilder();
		String[] tRegions = tAnnotations.split("\\|");
		
		Hashtable<String, Boolean> isDuplicated = new Hashtable<String, Boolean>();
		for(String tRegion : tRegions) {
			String transcriptID = tRegion.split("\\(")[0];
			String geneName = ENSTMapper.getGeneNamebyENST(transcriptID);
			
			if(isDuplicated.get(geneName) == null && !geneName.equalsIgnoreCase(Constants.ID_NULL)) {
				isDuplicated.put(geneName, true);
				gRegions.append("|").append(geneName);
			}
		}
		
		// if either unknown or intergenic, then there is no matched geneNames.
		// 
		if(gRegions.length() == 0) {
			mapper.put("key", Constants.ID_NULL);
			mapper.put("count", "0");
		} else {
			String annotation = gRegions.substring(1).toString();
			mapper.put("key", annotation);
			mapper.put("count", annotation.split("\\|").length+"");
		}
		
		return mapper;
	}
	
	/**
	 * The return value is hashtable consisting of "key" and "count"<br>
	 * "key" represents the annotation of matched gene IDs. <br>
	 * "count" represents the number of gene IDs. <br>
	 * 
	 * @return
	 */
	private Hashtable<String, String> toGeneIDs () {
		Hashtable<String, String> mapper = new Hashtable<String, String>();
		StringBuilder gRegions = new StringBuilder();
		String[] tRegions = tAnnotations.split("\\|");
		
		Hashtable<String, Boolean> isDuplicated = new Hashtable<String, Boolean>();
		for(String tRegion : tRegions) {
			String transcriptID = tRegion.split("\\(")[0];
			String geneID = ENSTMapper.getENSGbyENST(transcriptID);
			
			if(isDuplicated.get(geneID) == null && !geneID.equalsIgnoreCase(Constants.ID_NULL)) {
				isDuplicated.put(geneID, true);
				gRegions.append("|").append(geneID);
			}
		}
		
		// if either unknown or intergenic, then there is no matched gene IDs.
		// 
		if(gRegions.length() == 0) {
			mapper.put("key", Constants.ID_NULL);
			mapper.put("count", "0");
		} else {
			String annotation = gRegions.substring(1).toString();
			mapper.put("key", annotation);
			mapper.put("count", annotation.split("\\|").length+"");
		}
		
		return mapper;
	}
	
	/**
	 * filter regions by priority. <br>
	 * and also set the best priority (smaller is better). <br>
	 * See Priority.java. <br>
	 */
	public void filterRegions () {
		String[] regions = tAnnotations.split("\\|");
		double[] penalties = new double[regions.length];
		
		for(int i=0; i<regions.length; i++) {
			penalties[i] = Priority.getRegionPenalty(regions[i], mutations);
			bestRegionPriority = Math.min(penalties[i], bestRegionPriority);
		}
		
		StringBuilder filteredAnnotations = new StringBuilder();
		
		Boolean[] bans = new Boolean[regions.length];
		for(int i=0; i<regions.length; i++) {
			if(penalties[i] == bestRegionPriority) {
				filteredAnnotations.append("|").append(regions[i]);
				bans[i] = false;
			} else {
				bans[i] = true;
			}
		}
		
		this.exonLenghts = ResultParser.getWithoutBanList(this.exonLenghts, bans);
		this.percentFullDistances = ResultParser.getWithoutBanList(this.percentFullDistances, bans);
		this.percentExonDistances = ResultParser.getWithoutBanList(this.percentExonDistances, bans);
		this.percentCDSDistances = ResultParser.getWithoutBanList(this.percentCDSDistances, bans);
		this.fromStartDistances = ResultParser.getWithoutBanList(this.fromStartDistances, bans);
		this.fromStopDistances = ResultParser.getWithoutBanList(this.fromStopDistances, bans);
		
		this.tAnnotations = filteredAnnotations.substring(1);
	}
}
