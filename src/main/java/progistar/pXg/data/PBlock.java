package progistar.pXg.data;

import java.util.Hashtable;

import progistar.pXg.constants.Parameters;

public class PBlock {

	// preprocessed
	private String[] record;
	private String pSeq;
	
	// after mapping
	// key: peptide with I!=L
	public Hashtable<String, XBlock> xBlocks = new Hashtable<String, XBlock>();
	
	public PBlock (String[] record, String pSeq) {
		this.record = record;
		this.pSeq = pSeq;
	}
	
	/**
	 * 
	 * Return peptide sequence with IL-replacement option.<br>
	 * 
	 * @return
	 */
	public String getPeptideSequence () {
		if(Parameters.leucineIsIsoleucine) {
			return this.pSeq.replace("I", "L");
		} else return this.pSeq;
	}
}
