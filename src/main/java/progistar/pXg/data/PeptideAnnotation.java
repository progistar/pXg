package progistar.pXg.data;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import progistar.pXg.constants.Parameters;

public class PeptideAnnotation {

	
	private String[] fields;
	public ArrayList<PBlock> pBlocks = new ArrayList<PBlock>();
	
	public void setFields (String[] fields) {
		this.fields = fields;
	}
	
	public static String getFastaQueryFilePath () {
		return Parameters.peptideFilePath+".query";
	}
	
	/**
	 * Write non-duplicated peptide sequences. <br>
	 * This file will be an input for nBLASTt. <br> 
	 * 
	 * 
	 */
	public void writeFastaQuery () {
		File file = new File(getFastaQueryFilePath());
		
		try {
			BufferedWriter BW = new BufferedWriter(new FileWriter(file));
			int index = 0;
			
			ArrayList<String> sequences = this.enumerateSequence();
			
			for(String sequence : sequences) {
				index++;
				
				BW.append(">peptide"+index); 
				BW.newLine();
				BW.append(sequence);
				BW.newLine();
			}
			
			BW.close();
		}catch(IOException ioe) {
			
		}
	}
	/**
	 * Get non-duplicated peptide sequences from PeptideAnnotation.<br>
	 * 
	 * 
	 * @return
	 */
	private ArrayList<String> enumerateSequence () {
		Hashtable<String, Boolean> checks = new Hashtable<String, Boolean>();
		ArrayList<String> sequences = new ArrayList<String>();
		
		// put peptide sequences into checks
		this.pBlocks.forEach(pBlock -> checks.put(pBlock.getPeptideSequence(), true));
		
		// add to ArrayList without duplications
		checks.forEach((k, v) -> { sequences.add(k); });
		
		return sequences;
	}
	
}
