package progistar.pXg.data;

import java.util.Hashtable;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;

public class PBlock implements Comparable<PBlock> {

	// preprocessed
	private String[] record;
	private String pSeq;
	
	public String[] fastaIDs;
	public double score;
	public double deltaScore;
	public byte psmStatus = Constants.PSM_STATUS_UNDEF;
	public boolean isCannonical = false;
	public double fdrRate;
	// after mapping
	// key: peptide with I!=L
	public int rank = -1;
	
	// pBlock stores both target and decoy xBlocks above the RNA threshold.
	public Hashtable<String, XBlock> targetXBlocks = new Hashtable<String, XBlock>();
	public Hashtable<String, XBlock> decoyXBlocks = new Hashtable<String, XBlock>();
	
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
	
	public String getSpecID () {
		String specID = record[Parameters.fileColumnIndex] +"|"+record[Parameters.scanColumnIndex]+"|"+record[Parameters.chargeColumnIndex];
		return specID;
	}
	
	/**
	 * return record. <br>
	 * 
	 */
	public String toString (int genomicID) {
		StringBuilder recordLine = new StringBuilder();
		
		// add spectrum ID
		recordLine.append(this.getSpecID()).append("\t");
		
		// add genomic ID
		recordLine.append(genomicID).append("\t");
		
		// TD labeling
		if(this.psmStatus == Constants.PSM_STATUS_DECOY) {
			recordLine.append("-1");
		} else {
			recordLine.append("1");
		}
		
		for(int i=0; i<this.record.length; i++) {
			recordLine.append("\t").append(this.record[i]);
		}
		
		return recordLine.toString();
	}

	@Override
	/**
	 * Higher scores and decoy last.<br>
	 * 
	 */
	public int compareTo(PBlock o) {
		if(this.score > o.score) {
			return -1;
		} else if(this.score < o.score) {
			return 1;
		} else if (this.psmStatus > o.psmStatus){
			return -1;
		} else if (this.psmStatus < o.psmStatus) {
			return 1;
		}
		return 0;
	}
	
}
