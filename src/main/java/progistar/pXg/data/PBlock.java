package progistar.pXg.data;

public class PBlock {

	private String[] record;
	private String pSeq;
	
	public PBlock (String[] record, String pSeq) {
		this.record = record;
		this.pSeq = pSeq;
	}
	
	public String getPeptideSequence () {
		return this.pSeq;
	}
}
