package progistar.pXg.data;

import progistar.pXg.constants.Parameters;

public class PBlock {

	private String[] record;
	private String pSeq;
	
	public PBlock (String[] record, String pSeq) {
		this.record = record;
		this.pSeq = pSeq;
	}
	
	public String getPeptideSequence () {
		if(Parameters.leucineIsIsoleucine) {
			return this.pSeq.replace("I", "L");
		} else return this.pSeq;
	}
}
