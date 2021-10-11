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
	public String tAnnotations		=	null;
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
		int geneCount = geneIDs.split("\\|").length;
		return peptideSequence +"\t"+genomicLocus+"\t"+strand+"\t"+genomicSequence+"\t"+mutations+"\t"+tAnnotations
				+"\t"+geneIDs+"\t"+toGeneNames()+"\t"+toEvents(geneIDs)+"\t"+geneCount+"\t"+targetReadCount+"\t"+decoyReadCount;
	}
	
	public static String toNullString () {
		return "-\t-\t-\t-\t-\t-\t-\t-\t0\t0";
	}
	
	/**
	 * get gene events from geneIDs. <br>
	 * 
	 * @param geneIDs
	 * @return
	 */
	public String toEvents (String geneIDs) {
		StringBuilder events= new StringBuilder();
		String[] genes = geneIDs.split("\\|");
		
		for(String gene : genes) {
			events.append("|").append(Priority.getRegionEvent(gene));
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
			String gRegion = tRegion.replace(transcriptID, geneName);
			
			if(isDuplicated.get(gRegion) == null) {
				isDuplicated.put(gRegion, true);
				gRegions.append("|").append(gRegion);
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
			String gRegion = tRegion.replace(transcriptID, geneID);
			
			if(isDuplicated.get(gRegion) == null) {
				isDuplicated.put(gRegion, true);
				gRegions.append("|").append(gRegion);
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
		double[] scores = new double[regions.length];
		
		for(int i=0; i<regions.length; i++) {
			scores[i] = Priority.getRegionScore(regions[i], mutations);
			bestRegionPriority = Math.min(scores[i], bestRegionPriority);
		}
		
		StringBuilder filteredAnnotations = new StringBuilder();
		
		for(int i=0; i<regions.length; i++) {
			if(scores[i] == bestRegionPriority) {
				filteredAnnotations.append("|").append(regions[i]);
			}
		}
		
		this.tAnnotations = filteredAnnotations.substring(1);
	}
}
