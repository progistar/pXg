package progistar.pXg.data;

import java.util.Hashtable;

import progistar.pXg.utils.ENSTMapper;
import progistar.pXg.utils.Priority;

public class XBlock {
	public int targetReadCount		=	0;
	public int decoyReadCount		=	0;
	public char strand				=	'+';
	public String genomicLocus		=	null;
	public String mutations 		=	null;
	public String genomicSequence	=	null;
	public String peptideSequence	=	null;
	public String tAnnotations		=	null; // transcript and additional annotations
	public String[] fastaIDs		=	null;
	public double bestRegionPriority 	= 	Double.MAX_VALUE;
	
	
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
		
		int transCount = tAnnotations.split("\\|").length;
		int geneIDCount = geneIDs.split("\\|").length;
		int geneNameCount = geneNames.split("\\|").length;
		int eventCount = events.split("\\|").length;
		int fastaIDCount = this.fastaIDs.length;
		
		return peptideSequence +"\t"+genomicLocus+"\t"
				+strand+"\t"+genomicSequence
				+"\t"+mutations
				+"\t"+tAnnotations+"\t"+transCount
				+"\t"+geneIDs+"\t"+geneIDCount
				+"\t"+geneNames+"\t"+geneNameCount
				+"\t"+events+"\t"+eventCount
				+"\t"+fastaIDs+"\t"+fastaIDCount
				+"\t"+targetReadCount;
	}
	
	public static String toNullString () {
		return "-\t-\t-\t-\t-\t-\t0\t-\t0\t-\t0\t-\t0\t-\t0\t0";
	}
	
	public String toFastaIDs () {
		if(this.fastaIDs.length == 0) {
			return "-";
		}
		
		StringBuilder fasta= new StringBuilder();
		
		for(String fastaID : this.fastaIDs) {
			fasta.append("|").append(fastaID);
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
