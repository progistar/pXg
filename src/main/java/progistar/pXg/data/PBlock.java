package progistar.pXg.data;

import java.util.Hashtable;

import progistar.pXg.constants.Parameters;

public class PBlock implements Comparable<PBlock> {

	// preprocessed
	private String[] record;
	private String pSeq;
	
	public String[] fastaIDs;
	public double score;
	public boolean isTarget = true;
	// after mapping
	// key: peptide with I!=L
	public Hashtable<String, XBlock> xBlocks = new Hashtable<String, XBlock>();
	
	public PBlock (String[] record, String pSeq) {
		this.record = record;
		this.pSeq = pSeq;
		this.fastaIDs = new String[0]; // zero size, initially.
		this.score = Double.parseDouble(record[Parameters.scoreColumnIndex]);
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
	
	public String getScanID () {
		String scanID = "";
		for(int i=0; i<Parameters.scanColumnIndices.length; i++) {
			if(i!=0) scanID += "|";
			scanID += record[Parameters.scanColumnIndices[i]];
		}
		return scanID;
	}
	
	/**
	 * return record. <br>
	 * 
	 */
	public String toString () {
		StringBuilder recordLine = new StringBuilder();
		
		recordLine.append(this.record[0]);
		for(int i=1; i<this.record.length; i++) {
			recordLine.append("\t").append(this.record[i]);
		}
		
		return recordLine.toString();
	}

	@Override
	public int compareTo(PBlock o) {
		if(this.score > o.score) {
			return -1;
		} else if(this.score < o.score) {
			return 1;
		} else if (this.isTarget && !o.isTarget){
			return -1;
		} else if (!this.isTarget && o.isTarget) {
			return 1;
		}
		return 0;
	}
}
