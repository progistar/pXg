package progistar.pXg.data;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import progistar.pXg.constants.Parameters;

public class PxGAnnotation {

	
	
	// the first key: peptide sequence without I/L consideration
	// the first value: xBlocks corresponding to the key 
	// and the second key: peptide sequence from nucleotides + "_" + genomic locus
	private Hashtable<String, Hashtable<String, XBlock>> xBlockMapper = new Hashtable<String, Hashtable<String, XBlock>>();
	private int maxNGSreadCount = 0;
	
	public void putXBlock (String pSeq, XBlock xBlock) {
		Hashtable<String, XBlock> xBlocks = this.xBlockMapper.get(pSeq);
		
		if(xBlocks == null) {
			xBlocks = new Hashtable<String, XBlock>();
			this.xBlockMapper.put(pSeq, xBlocks);
		}
		
		String key = xBlock.getKey();
		
		XBlock thisXBlock = xBlocks.get(key);
		if(thisXBlock == null) {
			thisXBlock = xBlock;
			xBlocks.put(key, thisXBlock);
		} else {
			thisXBlock.decoyReadCount += xBlock.decoyReadCount;
			thisXBlock.targetReadCount += xBlock.targetReadCount;
			
			maxNGSreadCount = Math.max(thisXBlock.decoyReadCount, maxNGSreadCount);
			maxNGSreadCount = Math.max(thisXBlock.targetReadCount, maxNGSreadCount);
		}
	}
	
	/**
	 * Filter regions below than user-specific p-value. <br>
	 * 
	 */
	public void filterByPvalueThreshold () {
		double[] pValueTable = getPvalueTable();
		
		int cutoffReads = 0;
		for(int i=1; i<pValueTable.length; i++) {
			if(pValueTable[i] < Parameters.ngsPvalue) {
				cutoffReads = i;
				break;
			}
		}
		
		final int cutoff = cutoffReads;
		this.xBlockMapper.forEach((pSeq, xBlocks) -> {
			ArrayList<String> removeList = new ArrayList<String>();
			xBlocks.forEach((key, xBlock) -> {
				if(xBlock.targetReadCount < cutoff) {
					removeList.add(key);
				}
			});
			
			// filter out
			removeList.forEach(key -> xBlocks.remove(key));
		});
	}
	
	/**
	 * Calculate empirical p-value. <br>
	 * 
	 * @return
	 */
	private double[] getPvalueTable () {
		double[] pValueTable = new double[this.maxNGSreadCount+1];
		double[] decoyTable = new double[this.maxNGSreadCount+1];
		double[] targetTable = new double[this.maxNGSreadCount+1];
		
		this.xBlockMapper.forEach((pSeq, xBlocks) -> {
			xBlocks.forEach((key, xBlock) -> {
				decoyTable[xBlock.decoyReadCount]++;
				targetTable[xBlock.targetReadCount]++;
			});
		});
		
		// calculate cumulative decoy & target distribution
		for(int i=decoyTable.length-2; i>0; i--) {
			decoyTable[i] = decoyTable[i] + decoyTable[i+1];
			targetTable[i] = targetTable[i] + targetTable[i+1];
		}
		
		// calculate empirical p-value
		for(int i=1; i<pValueTable.length; i++) {
			pValueTable[i] = decoyTable[i] / decoyTable[1];
		}
		
		System.out.println("reads\tp-value\tregions");
		for(int i=1; i<pValueTable.length; i++) {
			System.out.println(i+"\t"+pValueTable[i]+"\t"+targetTable[i]);
		}
		
		return pValueTable;
	}
	
	public void write (String fileName) {
		try {
			File file = new File(fileName);
			BufferedWriter BW = new BufferedWriter(new FileWriter(file));
			
			ArrayList<PBlock> pBlocks = PeptideAnnotation.pBlocks;
			
			
			BW.append(PeptideAnnotation.toFields()).append("\t");
			BW.append("tPeptide").append("\t");
			BW.append("Loci").append("\t");
			BW.append("Strand").append("\t");
			BW.append("Nucleotide").append("\t");
			BW.append("Mutations").append("\t");
			BW.append("Transcripts").append("\t");
			BW.append("Reads").append("\t");
			BW.append("rReads");
			BW.newLine();
			
			
			for(PBlock pBlock : pBlocks) {
				// peptide sequence without I/L consideration
				String key = pBlock.getPeptideSequence();
				
				Hashtable<String, XBlock> xBlocks = this.xBlockMapper.get(key);
				
				
				
				// there is no available mapping.
				if(xBlocks == null) {
					// skip
					
					//BW.append(pBlock.toString()).append("\t").append(XBlock.toNullString());
					//BW.newLine();
					
				} else {
					xBlocks.forEach((pSeq, xBlock) -> {
						try {
							BW.append(pBlock.toString()).append("\t").append(xBlock.toString());
							BW.newLine();
						}catch(IOException ioe) {
							
						}
					});
				}
			}
			
			BW.close();
		}catch(IOException ioe) {
			
		}
	}
	
	/**
	 * select top-scored PSM only. <br>
	 * 
	 */
	public void topScoreFilter () {
		Hashtable<String, ArrayList<PBlock>> pBlocksByScan = new Hashtable<String, ArrayList<PBlock>>();
		ArrayList<PBlock> pBlocks = PeptideAnnotation.pBlocks;
		
		// aggregate pBlocks by scanID
		pBlocks.forEach(pBlock -> {
			String scanID = pBlock.getScanID();
			ArrayList<PBlock> scanPBlocks = pBlocksByScan.get(scanID);
			if(scanPBlocks == null) {
				scanPBlocks = new ArrayList<PBlock>();
				pBlocksByScan.put(scanID, scanPBlocks);
			}
			scanPBlocks.add(pBlock);
		});
		
		// remove original pBlocks
		pBlocks.clear();
		
		pBlocksByScan.forEach((scanID, scanPBlocks) -> {
			// sort pBlocks by scores, decreasing order.
			Collections.sort(scanPBlocks);
			
			// cutoff
			double topScore = scanPBlocks.get(0).score;
			int size = scanPBlocks.size();
			for(int i=0; i<size; i++) {
				// add only top-score
				double deltaScore = topScore - scanPBlocks.get(i).score;
				if(deltaScore == 0) {
					pBlocks.add(scanPBlocks.get(i));
				}
			}
		});
	}
}
