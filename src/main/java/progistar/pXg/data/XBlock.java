package progistar.pXg.data;

import java.util.ArrayList;
import java.util.Hashtable;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.utils.ENSTMapper;
import progistar.pXg.utils.Priority;

public class XBlock {
	// Note that only unmapped reads store sequenceID.
	// Mapped reads do not need this information
	public String sequenceID		=	null;
	
	public int targetReadCount		=	0;
	public int mockReadCount		=	0;
	public char strand				=	'+';
	public String genomicLocus		=	null;
	public String mutations 		=	null;
	public String genomicSequence	=	null;
	public String peptideSequence	=	null;
	public String tAnnotations		=	null; // transcript and additional annotations
	public String[] fastaIDs		=	new String[0];
	public String fullReadSequence	=	null; // for unmapped read
	public double bestRegionPriority 	= 	Double.MAX_VALUE;
	
	
	// with the same key value block
	// the block will be null if there is no next sibling.
	public ArrayList<XBlock> siblingXBlocks	=	new ArrayList<XBlock>();
	
	public boolean isCannonical () {
		String events = toEvents().get("key");
		boolean isCannonical = false;
		
		// if there is matched known sequences.
		if(fastaIDs.length != 0) {
			isCannonical = true;
		}
		// wildtype
		if(mutations.equalsIgnoreCase("-")) {
			// proteincoding;sense
			if(events.equalsIgnoreCase(Constants.EVENT_PROTEINCODING+";"+Constants.EVENT_SENSE)) {
				isCannonical = true;
			}
		}
		
		return isCannonical;
	}
	
	public boolean isMapped () {
		if(this.fullReadSequence == null) {
			return true;
		}
		return false;
	}
	
	public String getKey () {
		return this.genomicSequence+"_"+this.genomicLocus;
	}
	
	/**
	 * return record. <br>
	 * 
	 */
	public String toString () {
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
		
		
		return peptideSequence +"\t"+genomicLocus+"\t"
				+strand+"\t"+genomicSequence
				+"\t"+mutations
				+"\t"+tAnnotations+"\t"+transCount
				+"\t"+geneIDs.get("key")+"\t"+geneIDs.get("count")
				+"\t"+geneNames.get("key")+"\t"+geneNames.get("count")
				+"\t"+events.get("key")+"\t"+events.get("count")
				+"\t"+fastaIDs.get("key")+"\t"+fastaIDs.get("count")
				+"\t"+targetReadCount;
	}
	
	/**
	 * {@link Deprecated}
	 * Do not print PSM without read mapping.<br>
	 * -- The original usage is to print out PSMs with no reads.<br>
	 * @param fastaIDs_
	 * @return
	 */
	public static String toNullString (String[] fastaIDs_) {
		
		String fastaIDs = XBlock.toFastaIDs(fastaIDs_);
		int fastaIDCount = fastaIDs.equalsIgnoreCase("-") ? 0 : fastaIDs.split("\\|").length;
		
		return "-\t-\t"
				+ "-\t-\t"
				+ "-\t-\t0\t"
				+ "-\t0\t"
				+ "-\t0\t"
				+ "-\t0\t"
				+ ""+fastaIDs+"\t"+fastaIDCount+"\t"
				+ "0";
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
	
	/**
	 * {@link Deprecated}
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
		
		for(int i=0; i<regions.length; i++) {
			if(penalties[i] == bestRegionPriority) {
				filteredAnnotations.append("|").append(regions[i]);
			}
		}
		
		this.tAnnotations = filteredAnnotations.substring(1);
	}
}
