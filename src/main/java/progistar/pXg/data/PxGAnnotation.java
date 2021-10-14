package progistar.pXg.data;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;

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
		double[][] pValueTable = getPvalueTable();
		
		int[] cutoffReads = new int[Parameters.maxPeptLen+1];
		for(int peptLen = Parameters.minPeptLen; peptLen <= Parameters.maxPeptLen; peptLen++) {
			for(int i=1; i<pValueTable[peptLen].length; i++) {
				if(pValueTable[peptLen][i] < Parameters.ngsPvalue) {
					cutoffReads[peptLen] = i;
					System.out.println("Cutoff for "+peptLen+" aa: "+i);
					break;
				}
			}
		}
		
		final int[] cutoffs = cutoffReads;
		ArrayList<String> zeroSizeList = new ArrayList<String>();
		this.xBlockMapper.forEach((pSeq, xBlocks) -> {
			ArrayList<String> removeList = new ArrayList<String>();
			xBlocks.forEach((key, xBlock) -> {
				if(xBlock.targetReadCount < cutoffs[pSeq.length()]) {
					removeList.add(key);
				}
			});
			
			// filter out
			removeList.forEach(key -> xBlocks.remove(key));
			
			if(xBlocks.size() == 0) {
				zeroSizeList.add(pSeq);
			}
			
		});
		
		// if all xBlocks are discarded?
		zeroSizeList.forEach(key -> this.xBlockMapper.remove(key));
		
		ArrayList<PBlock> pBlocks = PeptideAnnotation.pBlocks;
		
		// update pBlocks!
		for(int i=pBlocks.size()-1; i>=0; i--) {
			PBlock pBlock = pBlocks.get(i);
			// peptide sequence without I/L consideration
			String key = pBlock.getPeptideSequence();
			
			Hashtable<String, XBlock> xBlocks = this.xBlockMapper.get(key);
			if(xBlocks == null) {
				pBlocks.remove(i);
			}
		}
		
	}
	
	/**
	 * Calculate empirical p-value. <br>
	 * By length. <br>
	 * 
	 * @return
	 */
	private double[][] getPvalueTable () {
		double[][] pValueTable = new double[Parameters.maxPeptLen+1][this.maxNGSreadCount+1];
		double[][] decoyTable = new double[Parameters.maxPeptLen+1][this.maxNGSreadCount+1];
		double[][] targetTable = new double[Parameters.maxPeptLen+1][this.maxNGSreadCount+1];
		
		this.xBlockMapper.forEach((pSeq, xBlocks) -> {
			xBlocks.forEach((key, xBlock) -> {
				decoyTable[pSeq.length()][xBlock.decoyReadCount]++;
				targetTable[pSeq.length()][xBlock.targetReadCount]++;
			});
		});
		
		// calculate cumulative decoy & target distribution
		try {
			BufferedWriter BW = new BufferedWriter(new FileWriter(Parameters.statFilePath));
			BW.append("PeptideLength\tReadCount\tTarget\tDecoy");
			BW.newLine();
			for(int peptLen = Parameters.minPeptLen; peptLen <= Parameters.maxPeptLen; peptLen++) {
				for(int i=decoyTable[peptLen].length-2; i>0; i--) {
					decoyTable[peptLen][i] = decoyTable[peptLen][i] + decoyTable[peptLen][i+1];
					targetTable[peptLen][i] = targetTable[peptLen][i] + targetTable[peptLen][i+1];
					
					BW.append(peptLen+"\t"+i+"\t"+targetTable[peptLen][i+1]+"\t"+decoyTable[peptLen][i+1]);
					BW.newLine();
				}
			}
			BW.close();
		}catch (IOException ioe) {
			
		}
		
		// calculate empirical p-value
		for(int peptLen = Parameters.minPeptLen; peptLen <= Parameters.maxPeptLen; peptLen++) {
			if(decoyTable[peptLen][1] == 0) continue;
			
			for(int i=1; i<pValueTable[peptLen].length; i++) {
				pValueTable[peptLen][i] = decoyTable[peptLen][i] / decoyTable[peptLen][1];
			}
		}
		
		for(int peptLen = Parameters.minPeptLen; peptLen <= Parameters.maxPeptLen; peptLen++) {
			System.out.println(peptLen+"aa peptides: "+targetTable[peptLen][1]+" targets  "+decoyTable[peptLen][1]+" decoys");
		}
		
		
		return pValueTable;
	}
	
	public void write (String fileName) {
		try {
			File file = new File(fileName);
			BufferedWriter BW = new BufferedWriter(new FileWriter(file));
			
			ArrayList<PBlock> pBlocks = PeptideAnnotation.pBlocks;
			
			
			BW.append(PeptideAnnotation.toFields()).append("\t");
			BW.append("MultipleRegions").append("\t");
			BW.append("tPeptide").append("\t");
			BW.append("Loci").append("\t");
			BW.append("Strand").append("\t");
			BW.append("Nucleotide").append("\t");
			BW.append("Mutations").append("\t");
			BW.append("TranscriptIDs").append("\t");
			BW.append("GeneIDs").append("\t");
			BW.append("GeneNames").append("\t");
			BW.append("Events").append("\t");
			BW.append("GeneCount").append("\t");
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
							BW.append(pBlock.toString()).append("\t").append(xBlocks.size()+"\t").append(xBlock.toString());
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
	
	public void regionScoreFilter () {
		// filter regions in the same locus and nucleotides.
		Iterator<String> pSeqs = (Iterator<String>) this.xBlockMapper.keys();
		
		while(pSeqs.hasNext()) {
			String pSeq = pSeqs.next();
			Hashtable<String, XBlock> xBlocks = this.xBlockMapper.get(pSeq);
			if(xBlocks == null || xBlocks.size() == 0) continue;
			
			Iterator<String> keys = (Iterator<String>) xBlocks.keys();
			double minPriority = Double.MAX_VALUE;
			while(keys.hasNext()) {
				String key = keys.next();
				XBlock xBlock = xBlocks.get(key);
				
				// filter-out 
				xBlock.filterRegions();
				
				// find minPriority
				minPriority = Math.min(xBlock.bestRegionPriority, minPriority);
			}
			
			keys = (Iterator<String>) xBlocks.keys();
			ArrayList<String> discardList = new ArrayList<String>();
			while(keys.hasNext()) {
				String key = keys.next();
				XBlock xBlock = xBlocks.get(key);
				if(xBlock.bestRegionPriority > minPriority) {
					discardList.add(key);
				}
			}
			
			// remove higher priority than minPriority
			discardList.forEach(key -> {
				xBlocks.remove(key);
			});
		}
	}
}
