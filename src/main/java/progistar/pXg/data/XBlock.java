package progistar.pXg.data;

import java.util.Hashtable;

import progistar.pXg.constants.Constants;
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
	public String[] fastaIDs		=	null;
	public String fullReadSequence	=	null; // for unmapped read
	public double bestRegionPriority 	= 	Double.MAX_VALUE;
	
	public boolean isCannonical () {
		String events = toEvents();
		boolean isCannonical = false;
		// wildtype
		if(mutations.equalsIgnoreCase("-")) {
			// proteincoding;sense
			if(events.contains(Constants.EVENT_PROTEINCODING) && events.contains(Constants.EVENT_SENSE) && !events.contains(Constants.EVENT_AS)) {
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
		String geneIDs = toGeneIDs();
		String geneNames = toGeneNames();
		String events = toEvents();
		String fastaIDs = toFastaIDs();
		
		int transCount = tAnnotations.equalsIgnoreCase("-") ? 0 : tAnnotations.split("\\|").length;
		int geneIDCount = geneIDs.equalsIgnoreCase("-") ? 0 : geneIDs.split("\\|").length;
		int geneNameCount = geneNames.equalsIgnoreCase("-") ? 0 : geneNames.split("\\|").length;
		int eventCount = events.equalsIgnoreCase("-") ? 0 : events.split("\\|").length;
		int fastaIDCount = fastaIDs.equalsIgnoreCase("-") ? 0 : fastaIDs.split("\\|").length;
		
		
		return peptideSequence +"\t"+genomicLocus+"\t"
				+strand+"\t"+genomicSequence
				+"\t"+mutations
				+"\t"+tAnnotations+"\t"+transCount
				+"\t"+geneIDs+"\t"+geneIDCount
				+"\t"+geneNames+"\t"+geneNameCount
				+"\t"+events+"\t"+eventCount
				+"\t"+fastaIDs+"\t"+fastaIDCount
				+"\t"+targetReadCount
				+"\t"+mockReadCount;
	}
	
	public static String toNullString () {
		return "-\t-\t"
				+ "-\t-\t"
				+ "-\t-\t0\t"
				+ "-\t0\t"
				+ "-\t0\t"
				+ "-\t0\t"
				+ "-\t0\t"
				+ "0\t0";
	}
	
	public String toFastaIDs () {
		if(this.fastaIDs == null || this.fastaIDs.length == 0) {
			return "-";
		}
		
		StringBuilder fasta= new StringBuilder();
		Hashtable<String, Boolean> isDuplicated = new Hashtable<String, Boolean>();
		
		for(String fastaID : this.fastaIDs) {
			if(isDuplicated.get(fastaID) == null) {
				fasta.append("|").append(fastaID);
				isDuplicated.put(fastaID, true);
			}
		}
		
		return fasta.substring(1).toString();
	}
	
	/**
	 * get gene events from geneIDs. <br>
	 * 
	 * @param transcriptIDs
	 * @return
	 */
	public String toEvents () {
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
		
		return events.substring(1).toString();
	}
	
	/**
	 * Retrieve gene names corresponding to the transcript ID. <br> 
	 * Duplicated gene names are reported once. <br>
	 * 
	 * @return
	 */
	private String toGeneNames () {
		StringBuilder gRegions = new StringBuilder();
		String[] tRegions = tAnnotations.split("\\|");
		
		Hashtable<String, Boolean> isDuplicated = new Hashtable<String, Boolean>();
		for(String tRegion : tRegions) {
			String transcriptID = tRegion.split("\\(")[0];
			String geneName = ENSTMapper.getGeneNamebyENST(transcriptID);
			
			if(isDuplicated.get(geneName) == null) {
				isDuplicated.put(geneName, true);
				gRegions.append("|").append(geneName);
			}
			
		}
		
		return gRegions.substring(1).toString();
	}
	
	/**
	 * Retrieve gene ID corresponding to the transcript ID. <br>
	 * Duplicated gene IDs are reported once. <br>
	 * 
	 * @return
	 */
	private String toGeneIDs () {
		StringBuilder gRegions = new StringBuilder();
		String[] tRegions = tAnnotations.split("\\|");
		
		Hashtable<String, Boolean> isDuplicated = new Hashtable<String, Boolean>();
		for(String tRegion : tRegions) {
			String transcriptID = tRegion.split("\\(")[0];
			String geneID = ENSTMapper.getENSGbyENST(transcriptID);
			
			if(isDuplicated.get(geneID) == null) {
				isDuplicated.put(geneID, true);
				gRegions.append("|").append(geneID);
			}
		}
		
		return gRegions.substring(1).toString();
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
