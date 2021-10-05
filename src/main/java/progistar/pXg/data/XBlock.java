package progistar.pXg.data;

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
		return peptideSequence +"\t"+genomicLocus+"\t"+strand+"\t"+genomicSequence+"\t"+mutations+"\t"+tAnnotations+"\t"+targetReadCount+"\t"+decoyReadCount;
	}
	
	public static String toNullString () {
		return "-\t-\t-\t-\t-\t-\t0\t0";
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
