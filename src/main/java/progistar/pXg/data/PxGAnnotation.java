package progistar.pXg.data;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;

import progistar.pXg.constants.Constants;
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
			thisXBlock.mockReadCount += xBlock.mockReadCount;
			thisXBlock.targetReadCount += xBlock.targetReadCount;
			
			maxNGSreadCount = Math.max(thisXBlock.mockReadCount, maxNGSreadCount);
			maxNGSreadCount = Math.max(thisXBlock.targetReadCount, maxNGSreadCount);
		}
	}
	
	/**
	 * Filter regions below than user-specific p-value. <br>
	 * 
	 */
	public void estimatePvalueThreshold () {
		System.out.println("Calculating p-values...");
		long startTime = System.currentTimeMillis();
		
		double[][] pValueTable = getPvalueTable();
		
		int[] cutoffReads = new int[Parameters.maxPeptLen+1];
		for(int peptLen = Parameters.minPeptLen; peptLen <= Parameters.maxPeptLen; peptLen++) {
			for(int i=1; i<pValueTable[peptLen].length; i++) {
				if(pValueTable[peptLen][i] < Parameters.ngsPvalue) {
					cutoffReads[peptLen] = i-1;
					System.out.println("Cutoff for "+peptLen+" aa: "+(i-1));
					break;
				}
			}
		}
		
		final int[] cutoffs = cutoffReads;
		ArrayList<String> zeroSizeList = new ArrayList<String>();
		this.xBlockMapper.forEach((pSeq, xBlocks) -> {
			if(xBlocks.size() != 0) {
				ArrayList<String> removeList = new ArrayList<String>();
				
				Iterator<String> keys = (Iterator<String>) xBlocks.keys();
				
				// single peptide can be assigned to multiple loci
				while(keys.hasNext()) {
					String key = keys.next();
					XBlock xBlock = xBlocks.get(key);
					if(xBlock.targetReadCount <= cutoffs[pSeq.length()]) {
						// if debug mode turns on, do not filter out annotations by reads
						if(!Parameters.debugMode) {
							removeList.add(key);
						}
					}
				}
				
				// filter out peptide with NGS count below than cut-off
				removeList.forEach(key -> xBlocks.remove(key));
				
				if(xBlocks.size() == 0) {
					zeroSizeList.add(pSeq);
				}
			}
			
		});
		
		// if all xBlocks are discarded?
		zeroSizeList.forEach(key -> this.xBlockMapper.remove(key));
		
		// finally,
		// xBlockMapper contains peptide with significantly-mapped RNA reads.
		
		long endTime = System.currentTimeMillis();
		System.out.println("\tElapsed time: "+((endTime-startTime)/1000) + " sec");
		
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
				decoyTable[pSeq.length()][xBlock.mockReadCount]++;
				targetTable[pSeq.length()][xBlock.targetReadCount]++;
			});
		});
		
		// calculate cumulative decoy & target distribution
		try {
			BufferedWriter BW = new BufferedWriter(new FileWriter(Parameters.ngsStatFilePath));
			BW.append("PeptideLength\tReadCount\tExperiment\tMock");
			BW.newLine();
			for(int peptLen = Parameters.minPeptLen; peptLen <= Parameters.maxPeptLen; peptLen++) {
				BW.append(peptLen+"\t"+(decoyTable[peptLen].length-1)+"\t"+targetTable[peptLen][decoyTable[peptLen].length-1]+"\t"+decoyTable[peptLen][decoyTable[peptLen].length-1]);
				BW.newLine();
				
				for(int i=decoyTable[peptLen].length-2; i>=0; i--) {
					BW.append(peptLen+"\t"+i+"\t"+targetTable[peptLen][i]+"\t"+decoyTable[peptLen][i]);
					BW.newLine();
					
					decoyTable[peptLen][i] = decoyTable[peptLen][i] + decoyTable[peptLen][i+1];
					targetTable[peptLen][i] = targetTable[peptLen][i] + targetTable[peptLen][i+1];
				}
			}
			BW.close();
		}catch (IOException ioe) {
			
		}
		
		// calculate empirical p-value
		for(int peptLen = Parameters.minPeptLen; peptLen <= Parameters.maxPeptLen; peptLen++) {
			if(decoyTable[peptLen][0] == 0) continue;
			
			for(int i=0; i<pValueTable[peptLen].length; i++) {
				pValueTable[peptLen][i] = decoyTable[peptLen][i] / decoyTable[peptLen][0];
			}
		}
		/*
		for(int peptLen = Parameters.minPeptLen; peptLen <= Parameters.maxPeptLen; peptLen++) {
			System.out.println(peptLen+"aa peptides: "+targetTable[peptLen][0]+" for experiment,  "+decoyTable[peptLen][0]+" for mock");
		}
		*/
		
		return pValueTable;
	}
	
	/**
	 * 
	 * 
	 */
	public void markFasta () {
		
		System.out.print("Marking fasta ids... ");
		long startTime = System.currentTimeMillis();
		
		Fasta fasta = new Fasta(Parameters.proteinFastaPath);
		
		ArrayList<String> sequences = PeptideAnnotation.enumerateSequence();
		Hashtable<String, ArrayList<String>> matchedList = fasta.findAll(sequences);
		
		// assign fasta IDs to pBlock
		for(PBlock pBlock : PeptideAnnotation.pBlocks) {
			ArrayList<String> ids = matchedList.get(pBlock.getPeptideSequence());
			if(ids != null) {
				pBlock.fastaIDs = new String[ids.size()];
				for(int i=0; i<pBlock.fastaIDs.length; i++) {
					pBlock.fastaIDs[i] = ids.get(i);
				}
			}
		}
		
		long endTime = System.currentTimeMillis();
		System.out.println("\tElapsed time: "+((endTime-startTime)/1000) + " sec");
	}
	
	public void write (String fileName) {
		try {
			File file = new File(fileName);
			BufferedWriter BW = new BufferedWriter(new FileWriter(file));
			
			ArrayList<PBlock> pBlocks = PeptideAnnotation.pBlocks;
			
			
			BW.append(PeptideAnnotation.toFields()).append("\t");
			BW.append("Rank").append("\t");
			BW.append("GenomicLociCount").append("\t");
			BW.append("InferredPeptide").append("\t");
			BW.append("GenomicLoci").append("\t");
			BW.append("Strand").append("\t");
			BW.append("Nucleotide").append("\t");
			BW.append("Mutations").append("\t");
			BW.append("TranscriptIDs").append("\t");
			BW.append("TranscriptIDCount").append("\t");
			BW.append("GeneIDs").append("\t");
			BW.append("GeneIDCount").append("\t");
			BW.append("GeneNames").append("\t");
			BW.append("GeneNameCount").append("\t");
			BW.append("Events").append("\t");
			BW.append("EventCount").append("\t");
			BW.append("FastaIDs").append("\t");
			BW.append("FastaIDCount").append("\t");
			BW.append("Reads").append("\t");
			BW.append("mockReads").append("\t");
			BW.append("FDR");
			BW.newLine();
			
			
			for(PBlock pBlock : pBlocks) {
				// peptide sequence without I/L consideration
				String key = pBlock.getPeptideSequence();
				
				Hashtable<String, XBlock> xBlocks = this.xBlockMapper.get(key);
				
				// there is no available mapping.
				if(xBlocks != null) {
					xBlocks.forEach((pSeq, xBlock) -> {
						try {
							// assign fastaIDs.
							xBlock.fastaIDs = pBlock.fastaIDs;
							
							BW.append(pBlock.toString()).append("\t").append(pBlock.rank+"\t").append(xBlocks.size()+"\t").append(xBlock.toString()).append("\t"+pBlock.fdrRate);
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
	
	public void mockReadAssignPolicy () {
		byte policy = Parameters.mockPolicy;
		
		// default behavior
		if(policy == Constants.MOCK_ALL) {
			return;
		} else {
			
			ArrayList<PBlock> pBlocks = PeptideAnnotation.pBlocks;
			
			// update mock reads following policy
			for(int i=pBlocks.size()-1; i>=0; i--) {
				PBlock pBlock = pBlocks.get(i);
				// peptide sequence without I/L consideration
				String key = pBlock.getPeptideSequence();
				
				float nonZeroMockReadCount = 0;
				float mean = 0;
				int max = 0;
				
				Hashtable<String, XBlock> xBlocks = this.xBlockMapper.get(key);
				if(xBlocks != null) {
					
					if(xBlocks.size() == 0) continue;
					// calculate mean and max
					Iterator<String> xBlockKeys = (Iterator<String>) xBlocks.keys();
					while(xBlockKeys.hasNext()) {
						String xBlockKey = xBlockKeys.next();
						XBlock xBlock = xBlocks.get(xBlockKey);
						
						if(xBlock.mockReadCount > 0) {
							nonZeroMockReadCount++;
							mean += xBlock.mockReadCount;
							max = Math.max(xBlock.mockReadCount, max);
						}
					}
					
					if(nonZeroMockReadCount != 0) {
						mean /= nonZeroMockReadCount;
					}
					
					if(max != 0) {
						// assign max or mean by policy
						// note that:
						// we do not need preserve consistency between read-count and genomic locus
						// this is because mock reads are such a virtual read count.
						xBlockKeys = (Iterator<String>) xBlocks.keys();
						int index = 0;
						while(xBlockKeys.hasNext()) {
							String xBlockKey = xBlockKeys.next();
							XBlock xBlock = xBlocks.get(xBlockKey);
							
							if(xBlock.mockReadCount > 0) {
								if(index == 0) {
									if(policy == Constants.MOCK_MAX_ONE) {
										xBlock.mockReadCount = max;
									} else if(policy == Constants.MOCK_MEAN) {
										xBlock.mockReadCount = Math.round(mean);
									}
								} else {
									xBlock.mockReadCount = 0;
								}
								index++;
							}
						}
					}
				}
			}
		}
	}
	
	/**
	 * select top-scored PSM only. <br>
	 * 
	 */
	public void topScoreFilter () {
		// if debug-mode turns on, print all of PSMs regardless of filtering
		if(Parameters.debugMode) return;
		
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
			
			// topScore is determined by target or decoy status
			int bestIndex = 0;
			byte bestStatus = Constants.PSM_STATUS_RANDOM;
			for(int i=0; i<scanPBlocks.size(); i++) {
				
				byte psmStatus = scanPBlocks.get(i).psmStatus;
				
				if(psmStatus == Constants.PSM_STATUS_TARGET) {
					//select top score from targets
					bestIndex = i;
					bestStatus = psmStatus;
					break;
				} else if(psmStatus == Constants.PSM_STATUS_DECOY && psmStatus > bestStatus) {
					bestStatus = psmStatus;
					bestIndex = i;
				} else if(bestStatus == Constants.PSM_STATUS_RANDOM){
					// worst score is selected to random PSMs
					bestIndex = i;
				}
			}
			
			pBlocks.add(scanPBlocks.get(bestIndex));
		});
	}
	/**
	 * Marking target PSM <br>
	 * 
	 */
	public void markTargetPSM () {
		long startTime = System.currentTimeMillis();
		
		ArrayList<PBlock> pBlocks = PeptideAnnotation.pBlocks;
		
		// update pBlocks!
		for(int i=pBlocks.size()-1; i>=0; i--) {
			PBlock pBlock = pBlocks.get(i);
			// peptide sequence without I/L consideration
			String key = pBlock.getPeptideSequence();
			
			Hashtable<String, XBlock> xBlocks = this.xBlockMapper.get(key);
			if(xBlocks != null) {
				
				xBlocks.forEach((key_, xBlock) -> {
					if(xBlock.targetReadCount > 0) {
						pBlock.psmStatus = Constants.PSM_STATUS_TARGET > pBlock.psmStatus ? Constants.PSM_STATUS_TARGET : pBlock.psmStatus;
					} 
					// entrapment psms
					else if(xBlock.mockReadCount > 0) {
						pBlock.psmStatus = Constants.PSM_STATUS_DECOY > pBlock.psmStatus ? Constants.PSM_STATUS_DECOY : pBlock.psmStatus;
					}
				});
			}
		}
		
		long endTime = System.currentTimeMillis();
		System.out.println("\tElapsed time: "+((endTime-startTime)/1000) + " sec");
		
	}
	
	/**
	 * 
	 * 
	 */
	public void fdrEstimation () {
		double targetCount = 0;
		double randomCount = 0;
		double decoyCount = 0;
		int fdrCutoffIndex = 0;
		
		// Assume that, single PSM per scan.
		// This is because topScoreFiler only selects single PSM per scan.
		ArrayList<PBlock> pBlocks = PeptideAnnotation.pBlocks;
		
		// sort pBlocks by descending order.
		Collections.sort(pBlocks);
		int size = pBlocks.size();
		
		try {
			BufferedWriter BW = new BufferedWriter(new FileWriter(Parameters.psmStatFilePath));
			BW.append("Class\tScore\tFDR\tTargetCount\tDecoyCount\tRandomCount");
			BW.newLine();
			for(int i=0; i<size; i++) {
				double fdrRate = 1.0;
				
				PBlock pBlock = pBlocks.get(i);
				
				if(pBlock.psmStatus == Constants.PSM_STATUS_TARGET) {
					targetCount++;
					BW.append("target");
				}
				else if(pBlock.psmStatus == Constants.PSM_STATUS_RANDOM){
					randomCount++;
					BW.append("random");
				}
				else if(pBlock.psmStatus == Constants.PSM_STATUS_DECOY) {
					decoyCount++;
					BW.append("decoy");
				}
				
				if(targetCount != 0 || decoyCount != 0) {
					fdrRate = (randomCount*0.5 + decoyCount)/(targetCount+randomCount*0.5);
				}
				
				BW.append("\t"+pBlock.score+"\t"+fdrRate+"\t"+targetCount+"\t"+decoyCount+"\t"+randomCount);
				BW.newLine();
				
				if(fdrRate < Parameters.fdrThreshold) {
					fdrCutoffIndex = i;
				}
				
				pBlock.fdrRate = fdrRate;
			}
			BW.close();
		}catch(IOException ioe) {
			
		}
		
		// remove below than fdr cutoff
		for(int i=pBlocks.size()-1; i>=0; i--) {
			PBlock pBlock = pBlocks.get(i);
			
			if(pBlock.psmStatus != Constants.PSM_STATUS_TARGET) {
				if(!Parameters.debugMode) {
					pBlocks.remove(i);
				}
			} else if(i >= fdrCutoffIndex) {
				if(!Parameters.debugMode) {
					pBlocks.remove(i);
				}
			}
			
		}
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
