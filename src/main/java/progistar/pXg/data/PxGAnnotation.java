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
		}
		
		maxNGSreadCount = Math.max(thisXBlock.mockReadCount, maxNGSreadCount);
		maxNGSreadCount = Math.max(thisXBlock.targetReadCount, maxNGSreadCount);
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
		double[][] mockTable = new double[Parameters.maxPeptLen+1][this.maxNGSreadCount+1];
		double[][] expTable = new double[Parameters.maxPeptLen+1][this.maxNGSreadCount+1];
		// all peptides
		ArrayList<String> peptides = PeptideAnnotation.enumerateSequence();
		
		for(String peptide : peptides) {
			Hashtable<String, XBlock> xBlocks = this.xBlockMapper.get(peptide);
			int length = peptide.length();
			int[] mockCounts = new int[1];
			int[] expCounts = new int[1];
			if(xBlocks != null) {
				mockCounts = new int[xBlocks.size()];
				expCounts = new int[xBlocks.size()];
				Iterator<String> keys = (Iterator<String>) xBlocks.keys();
				int idx = 0;
				while(keys.hasNext()) {
					String key =keys.next();
					XBlock xBlock = xBlocks.get(key);
					
					if(xBlock.mockReadCount > 0) {
						mockCounts[idx] = xBlock.mockReadCount;
					}
					if(xBlock.targetReadCount > 0) {
						expCounts[idx] = xBlock.targetReadCount;
					}
					idx++;
				}
				
				
			}
			
			if(Parameters.mockPolicy == Constants.MOCK_ALL) {
				// mock counts for all possible regions
				boolean isZeroCount = true;
				for(int mockCount : mockCounts) {
					if(mockCount > 0) {
						mockTable[length][mockCount]++;
						isZeroCount = false;
					}
				}
				// there is no matched region for the peptide,
				// just count zero at once.
				if(isZeroCount) {
					mockTable[length][0]++;
				}
				
				// exp counts for all possible regions
				isZeroCount = true;
				for(int expCount : expCounts) {
					if(expCount > 0) {
						expTable[length][expCount]++;
						isZeroCount = false;
					}
				}
				if(isZeroCount) {
					expTable[length][0]++;
				}
			} else if(Parameters.mockPolicy == Constants.MOCK_MAX_ONE) {
				// find max one through traversing counts
				int maxCount = 0;
				for(int mockCount : mockCounts) {
					maxCount = Math.max(mockCount, maxCount);
				}
				mockTable[length][maxCount]++;
				
				maxCount = 0;
				for(int expCount : expCounts) {
					maxCount = Math.max(expCount, maxCount);
				}
				expTable[length][maxCount]++;
				
			} else if(Parameters.mockPolicy == Constants.MOCK_MEAN) {
				// sum all counts through traversing counts
				// and get the average values
				double count = 0;
				double size = 0;
				for(int mockCount : mockCounts) {
					if(mockCount > 0) {
						count += mockCount;
						size ++;
					}
				}
				if(size != 0) {
					count = Math.round(count/size);
				}
				mockTable[length][(int)count]++;
				
				count = 0;
				size = 0;
				for(int expCount : expCounts) {
					if(expCount > 0) {
						count += expCount;
						size ++;
					}
				}
				if(size != 0) {
					count = Math.round(count/size);
				}
				expTable[length][(int)count]++;
				
			}
		}
		
		// calculate cumulative decoy & target distribution
		try {
			BufferedWriter BW = new BufferedWriter(new FileWriter(Parameters.ngsStatFilePath));
			BW.append("PeptideLength\tReadCount\tExperiment\tMock");
			BW.newLine();
			for(int peptLen = Parameters.minPeptLen; peptLen <= Parameters.maxPeptLen; peptLen++) {
				BW.append(peptLen+"\t"+(mockTable[peptLen].length-1)+"\t"+expTable[peptLen][mockTable[peptLen].length-1]+"\t"+mockTable[peptLen][mockTable[peptLen].length-1]);
				BW.newLine();
				
				for(int i=mockTable[peptLen].length-2; i>=0; i--) {
					BW.append(peptLen+"\t"+i+"\t"+expTable[peptLen][i]+"\t"+mockTable[peptLen][i]);
					BW.newLine();
					
					mockTable[peptLen][i] = mockTable[peptLen][i] + mockTable[peptLen][i+1];
					expTable[peptLen][i] = expTable[peptLen][i] + expTable[peptLen][i+1];
				}
			}
			BW.close();
		}catch (IOException ioe) {
			
		}
		
		// calculate empirical p-value
		for(int peptLen = Parameters.minPeptLen; peptLen <= Parameters.maxPeptLen; peptLen++) {
			if(mockTable[peptLen][0] == 0) continue;
			
			for(int i=0; i<pValueTable[peptLen].length; i++) {
				pValueTable[peptLen][i] = mockTable[peptLen][i] / mockTable[peptLen][0];
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
				// assign rank
				scanPBlocks.get(i).rank = (i+1);
				
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
						pBlock.isCannonical |= xBlock.isCannonical();
					} 
					// entrapment psms
					else if(xBlock.mockReadCount > 0) {
						pBlock.psmStatus = Constants.PSM_STATUS_DECOY > pBlock.psmStatus ? Constants.PSM_STATUS_DECOY : pBlock.psmStatus;
						pBlock.isCannonical |= xBlock.isCannonical();
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
		double decoyCount = 0;
		int fdrCutoffIndex = 0;
		
		// Assume that, single PSM per scan.
		// This is because topScoreFiler only selects single PSM per scan.
		ArrayList<PBlock> pBlocks = PeptideAnnotation.pBlocks;
		
		// sort pBlocks by descending order.
		Collections.sort(pBlocks);
		int size = pBlocks.size();
		
		try {
			ArrayList<Double> scores = new ArrayList<Double>();
			Hashtable<Double, Boolean> scoreList = new Hashtable<Double, Boolean>();
			Hashtable<Double, Double> cTargetCounts = new Hashtable<Double, Double>();
			Hashtable<Double, Double> cDecoyCounts = new Hashtable<Double, Double>();
			Hashtable<Double, Double> ncTargetCounts = new Hashtable<Double, Double>();
			Hashtable<Double, Double> ncDecoyCounts = new Hashtable<Double, Double>();
			
			for(int i=0; i<size; i++) {
				double fdrRate = 1.0;
				
				PBlock pBlock = pBlocks.get(i);
				
				byte case_ = pBlock.psmStatus;
				double score = pBlock.score;
				if(scoreList.get(score) == null) {
					scores.add(score);
					scoreList.put(score, true);
				}
				
				if(case_ == Constants.PSM_STATUS_TARGET) {
					targetCount++;
					Double count = .0;
					if(pBlock.isCannonical) {
						count = cTargetCounts.get(score);
					} else {
						count = ncTargetCounts.get(score);
					}
					if(count == null) {
						count = 0.0;
					}
					count++;
					if(pBlock.isCannonical) {
						cTargetCounts.put(score, count);
					} else {
						ncTargetCounts.put(score, count);
					}
				}
				else if(case_ == Constants.PSM_STATUS_DECOY) {
					decoyCount++;
					Double count = .0;
					if(pBlock.isCannonical) {
						count = cDecoyCounts.get(score);
					} else {
						count = ncDecoyCounts.get(score);
					}
					if(count == null) {
						count = 0.0;
					}
					count++;
					if(pBlock.isCannonical) {
						cDecoyCounts.put(score, count);
					} else {
						ncDecoyCounts.put(score, count);
					}
					
				}
				
				if(targetCount != 0) {
					fdrRate = decoyCount/targetCount;
				}
				if(fdrRate < Parameters.fdrThreshold) {
					fdrCutoffIndex = i;
				}
				
				pBlock.fdrRate = fdrRate;
			}
			
			BufferedWriter BW = new BufferedWriter(new FileWriter(Parameters.psmStatFilePath));
			BW.append("Score\tcTargetCount\tcDecoyCount\tncTargetCount\tncDecoyCount");
			BW.newLine();
			
			for(int i=0; i<scores.size(); i++) {
				Double score = scores.get(i);
				Double ctCount = cTargetCounts.get(score);
				Double cdCount = cDecoyCounts.get(score);
				Double nctCount = ncTargetCounts.get(score);
				Double ncdCount = ncDecoyCounts.get(score);
				
				
				if(ctCount == null) {
					ctCount = 0.0;
				}
				if(cdCount == null) {
					cdCount = 0.0;
				}
				if(nctCount == null) {
					nctCount = 0.0;
				}
				if(ncdCount == null) {
					ncdCount = 0.0;
				}
				
				BW.append(score+"\t").append(ctCount+"\t").append(cdCount+"\t").append(nctCount+"\t").append(ncdCount+"");
				BW.newLine();
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
			} else if(i > fdrCutoffIndex) {
				if(!Parameters.debugMode) {
					pBlocks.remove(i);
				}
			}
			
		}
	}
	
	public void regionScoreFilter () {
		// prevent zeor-size
		if(this.xBlockMapper.size() == 0) return;
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
