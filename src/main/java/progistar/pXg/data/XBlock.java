package progistar.pXg.data;

public class XBlock {
	public int targetReadCount		=	0;
	public int decoyReadCount		=	0;
	public char strand				=	'+';
	public String genomicLocus		=	null;
	public String mutations 		=	null;
	public String genomicSequence	=	null;
	public String peptideSequence	=	null;
	public String tAnnotations		=	null;
	
	
	public String getKey () {
		return this.peptideSequence+"_"+this.genomicLocus;
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
}
