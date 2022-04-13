package progistar.pXg.data;

import java.awt.Point;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Iterator;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.constants.RunInfo;
import progistar.pXg.utils.Logger;

public class PxGAnnotation {

	
	
	// the first key: peptide sequence without I/L consideration
	// the first value: xBlocks corresponding to the key 
	// and the second key: peptide sequence from nucleotides + "_" + genomic locus
	private Hashtable<String, Hashtable<String, XBlock>> xBlockMapper = new Hashtable<String, Hashtable<String, XBlock>>();
	private int maxNGSreadCount = 0;
	
	/**
	 * Do not modify by this return object.<br>
	 * Strongly prevent modifying object! Look-up is okay. <br>
	 * 
	 * @return
	 */
	public Hashtable<String, Hashtable<String, XBlock>> getXBlockMapper () {
		return this.xBlockMapper;
	}
	
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
			thisXBlock.siblingXBlocks.add(xBlock); // 
		}
		
		maxNGSreadCount = Math.max(thisXBlock.mockReadCount, maxNGSreadCount);
		maxNGSreadCount = Math.max(thisXBlock.targetReadCount, maxNGSreadCount);
	}
	
	/**
	 * Filter regions below than user-specific p-value. <br>
	 * As a result, there are only two possible xBlocks are remained. <br>
	 * 1) xBlocks with passing the threshold. <br>
	 * 2) xBlocks with only matched to mock reads. <br>
	 * <br>
	 * Note that the first case will be regarded as target PSMs; <br>
	 * and the second one will be regarded as decoy PSMs<br>
	 * <br>
	 * We enforce xBlocks with passing the threshold to have zero-mock reads. <br>
	 * This is because existence of mock read count will be used to recognize target and decoy PSMs. <br>
	 * And then, pSeq must be determined by exp read or mock read.<br> 
	 */
	public void estimatePvalueThreshold () {
		System.out.println("Calculating p-values...");
		Logger.append("Calculating p-values...");
		Logger.newLine();
		long startTime = System.currentTimeMillis();
		
		double[][] pValueTable = getPvalueTable();
		
		RunInfo.cutoffReads = new int[Parameters.maxPeptLen+1];
		for(int peptLen = Parameters.minPeptLen; peptLen <= Parameters.maxPeptLen; peptLen++) {
			for(int i=1; i<pValueTable[peptLen].length; i++) {
				if(pValueTable[peptLen][i] < Parameters.pvalue) {
					RunInfo.cutoffReads[peptLen] = i;
					System.out.println("Minimum reads threshold to accept at "+peptLen+" aa peptide: "+(i));
					// append to logger
					Logger.append("Minimum reads threshold to accept at "+peptLen+" aa peptide: "+(i));
					Logger.newLine();
					break;
				}
			}
		}
		
		/**
		 * After filter xBlock by threshold, if there is at least one xBlock with exp read,<br>
		 * then the all xBlocks with mock will be discarded in corresponding to pSeq.<br>
		 * 
		 */
		final int[] cutoffs = RunInfo.cutoffReads;
		ArrayList<String> zeroSizeList = new ArrayList<String>();
		
		this.xBlockMapper.forEach((pSeq, xBlocks) -> {
			if(xBlocks.size() != 0) {
				ArrayList<String> removeList = new ArrayList<String>();
				
				Iterator<String> keys = (Iterator<String>) xBlocks.keys();
				
				// single peptide can be assigned to multiple loci
				while(keys.hasNext()) {
					String key = keys.next();
					XBlock xBlock = xBlocks.get(key);
					
					if(xBlock.targetReadCount == 0 && xBlock.mockReadCount > 0) {
						continue;
					}
					
					if(xBlock.targetReadCount < cutoffs[pSeq.length()]) {
						// if debug mode turns on, do not filter out annotations by reads
						if(!Parameters.debugMode) {
							removeList.add(key);
						}
					} else {
						xBlock.mockReadCount = 0;
					}
				}
				
				// filter out peptide with NGS count below than cut-off
				removeList.forEach(key -> xBlocks.remove(key));
				removeList.clear();
				
				// pSeq will be determined by either exp read or mock read.
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
				// 22.01.24
				// if the pBlock has fasta!
				pBlock.isCannonical = true;
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
			BW.append("IsCanonical");
			BW.newLine();
			
			File outFile = new File(Parameters.unmappedFilePath);
			BufferedWriter BWUnmapped = new BufferedWriter(new FileWriter(outFile));
			
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
							
							BW.append(pBlock.toString()).append("\t").append(pBlock.rank+"\t").append(xBlocks.size()+"\t").append(xBlock.toString()).append("\t"+pBlock.isCannonical);
							BW.newLine();
							
							// if this is unmapped, then store.
							if(!xBlock.isMapped()) {
								ArrayList<XBlock> unmappedXBlocks = new ArrayList<XBlock>();
								unmappedXBlocks.add(xBlock);
								unmappedXBlocks.addAll(xBlock.siblingXBlocks);
								
								BWUnmapped.append(">"+xBlock.peptideSequence);
								BWUnmapped.newLine();
								for(XBlock thisXBlock : unmappedXBlocks) {
									BWUnmapped.append(thisXBlock.sequenceID).append("\t").append(thisXBlock.fullReadSequence).append("\t").append(thisXBlock.genomicSequence);
									BWUnmapped.newLine();
								}
								
							}
							
						}catch(IOException ioe) {
							
						}
					});
				} 
			}
			
			BW.close();
			BWUnmapped.close();
			
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
			int bestTargetIndex = -1;
			int bestDecoyIndex = -1;
			
			for(int i=0; i<scanPBlocks.size(); i++) {
				
				byte psmStatus = scanPBlocks.get(i).psmStatus;
				
				if(psmStatus == Constants.PSM_STATUS_TARGET) {
					//select top score from targets
					bestTargetIndex = i;
					break;
				} else if(psmStatus == Constants.PSM_STATUS_DECOY) {
					if(bestDecoyIndex == -1) {
						bestDecoyIndex = i;
					}
				}
			}
			
			// top-scored PSMs have priority
			if(bestTargetIndex != -1 && bestDecoyIndex != -1) {
				if(scanPBlocks.get(bestTargetIndex).score >= scanPBlocks.get(bestDecoyIndex).score) {
					pBlocks.add(scanPBlocks.get(bestTargetIndex));
				} else {
					pBlocks.add(scanPBlocks.get(bestDecoyIndex));
				}
				
			} else if(bestTargetIndex != -1) {
				pBlocks.add(scanPBlocks.get(bestTargetIndex));
			} else if(bestDecoyIndex != -1){
				pBlocks.add(scanPBlocks.get(bestDecoyIndex));
			}
		});
	}
	/**
	 * Marking target PSM <br>
	 * 
	 */
	public void markTargetPSMs () {
		long startTime = System.currentTimeMillis();
		
		ArrayList<PBlock> pBlocks = PeptideAnnotation.pBlocks;
		
		// update pBlocks!
		for(int i=pBlocks.size()-1; i>=0; i--) {
			PBlock pBlock = pBlocks.get(i);
			// peptide sequence without I/L consideration
			String key = pBlock.getPeptideSequence();
			
			Hashtable<String, XBlock> xBlocks = this.xBlockMapper.get(key);
			if(xBlocks != null) {
				
				// check error!
				boolean[] expAndMocks = new boolean[2];
				// only select significantly mapped PSMs
				// this is because we are interested in P( decoy (score > X) | significantly mapped PSMs)
				xBlocks.forEach((key_, xBlock) -> {
					if(xBlock.targetReadCount >= RunInfo.cutoffReads[key.length()]) {
						pBlock.psmStatus = Constants.PSM_STATUS_TARGET > pBlock.psmStatus ? Constants.PSM_STATUS_TARGET : pBlock.psmStatus;
						pBlock.isCannonical |= xBlock.isCannonical();
						expAndMocks[0] = true;
					} 
					// decoy PSMs
					else if(xBlock.mockReadCount >= RunInfo.cutoffReads[key.length()]) {
						pBlock.psmStatus = Constants.PSM_STATUS_DECOY > pBlock.psmStatus ? Constants.PSM_STATUS_DECOY : pBlock.psmStatus;
						pBlock.isCannonical |= xBlock.isCannonical();
						expAndMocks[1] = true;
					}
				});
				
				if(expAndMocks[0] && expAndMocks[1]) {
					System.out.println("markTargetPSMs was being invalid behavior!::");
					System.out.println("-- exp and mock read xBlocks were coincident!");
				}
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
		double cTargetCount = 0;
		double ncTargetCount = 0;
		double cDecoyCount = 0;
		double ncDecoyCount = 0;
		int cFDRCutoffIndex = 0;
		int ncFDRCutoffIndex = 0;
		
		// Assume that, single PSM per scan.
		// This is because topScoreFiler only selects single PSM per scan.
		ArrayList<PBlock> pBlocks = PeptideAnnotation.pBlocks;
		
		// sort pBlocks by descending order.
		class decoyFirstOrder implements Comparator<PBlock> {
			  @Override
			  public int compare(PBlock p1, PBlock p2) {
				  if(p1.score > p2.score) {
						return -1;
					} else if(p1.score < p2.score) {
						return 1;
					} else if (p1.psmStatus > p2.psmStatus){
						return 1;
					} else if (p1.psmStatus < p2.psmStatus) {
						return -1;
					}
					return 0;
			  }
		}
		Collections.sort(pBlocks, new decoyFirstOrder());
		
		int size = pBlocks.size();
		
		try {
			ArrayList<Double> scores = new ArrayList<Double>();
			Hashtable<Double, Boolean> isScored = new Hashtable<Double, Boolean>(); // to print score at once.
			Hashtable<Double, Double> cTargetCounts = new Hashtable<Double, Double>();
			Hashtable<Double, Double> cDecoyCounts = new Hashtable<Double, Double>();
			Hashtable<Double, Double> ncTargetCounts = new Hashtable<Double, Double>();
			Hashtable<Double, Double> ncDecoyCounts = new Hashtable<Double, Double>();
			
			for(int i=0; i<size; i++) {
				double fdrRate = 1.0;
				
				PBlock pBlock = pBlocks.get(i);
				
				byte case_ = pBlock.psmStatus;
				
				if(case_ == Constants.PSM_STATUS_UNDEF) {
					continue;
				}
				
				double score = pBlock.score;
				
				if(case_ == Constants.PSM_STATUS_TARGET) {
					Double count = .0;
					if(pBlock.isCannonical) {
						cTargetCount++;
						count = cTargetCounts.get(score);
					} else {
						ncTargetCount++;
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
					Double count = .0;
					if(pBlock.isCannonical) {
						cDecoyCount++;
						count = cDecoyCounts.get(score);
					} else {
						ncDecoyCount++;
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
				
				if(case_ == Constants.PSM_STATUS_TARGET) {
					if(pBlock.isCannonical) {
						// for canonical cutoff
						if(cTargetCount != 0) {
							fdrRate = cDecoyCount/cTargetCount;
						}
						if(fdrRate < Parameters.fdr) {
							cFDRCutoffIndex = i;
							RunInfo.cPSMScoreTreshold = pBlock.score;
						}
					} else {
						// for noncanonical cutoff
						if(ncTargetCount != 0) {
							fdrRate = ncDecoyCount/ncTargetCount;
						}
						if(fdrRate < Parameters.fdr) {
							ncFDRCutoffIndex = i;
							RunInfo.ncPSMScoreTreshold = pBlock.score;
						}
					}
					
					pBlock.fdrRate = fdrRate;
				}
				
				
				if(isScored.get(score) == null) {
					scores.add(score);
					isScored.put(score, true);
				}
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
		
		System.out.println("cFDR cutoff idx: "+cFDRCutoffIndex);
		System.out.println("ncFDR cutoff idx: "+ncFDRCutoffIndex);
		System.out.println("cTarget counts: "+cTargetCount);
		System.out.println("cDecoy counts: "+cDecoyCount);
		System.out.println("ncTarget counts: "+ncTargetCount);
		System.out.println("ncDecoy counts: "+ncDecoyCount);
		
		// remove below than FDR cutoff
		ArrayList<PBlock> cutoffedPBlocks = new ArrayList<PBlock>();
		for(int i=0; i<pBlocks.size(); i++) {
			PBlock pBlock = pBlocks.get(i);
			
			if(pBlock.psmStatus == Constants.PSM_STATUS_TARGET) {
				if(pBlock.isCannonical) {
					if(pBlock.score >= RunInfo.cPSMScoreTreshold) {
						cutoffedPBlocks.add(pBlock);
					}
				} else {
					if(pBlock.score >= RunInfo.ncPSMScoreTreshold) {
						cutoffedPBlocks.add(pBlock);
					}
				}
			}
		}
		
		// update
		PeptideAnnotation.pBlocks = cutoffedPBlocks;
	}
	
	/**
	 * This method expects that: <br>
	 * 1) following estimatePvalueTreshold.<br>
	 * 2) so, pSeq has exp or mock read only.<br>
	 * 
	 * 
	 */
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
			
			double minExpXBlockPenalty = Double.MAX_VALUE;
			double minMockXBlockPenalty = Double.MAX_VALUE;
			boolean isExpXBlock = false;
			
			while(keys.hasNext()) {
				String key = keys.next();
				XBlock xBlock = xBlocks.get(key);
				
				// filter-out 
				xBlock.filterRegions();
				
				// find min penalty
				if(xBlock.targetReadCount > 0) {
					minExpXBlockPenalty = Math.min(xBlock.bestRegionPriority, minExpXBlockPenalty);
					isExpXBlock = true;
				} else if(xBlock.mockReadCount > 0) {
					minMockXBlockPenalty = Math.min(xBlock.bestRegionPriority, minMockXBlockPenalty);
				}
			}
			
			keys = (Iterator<String>) xBlocks.keys();
			ArrayList<String> discardList = new ArrayList<String>();
			while(keys.hasNext()) {
				String key = keys.next();
				XBlock xBlock = xBlocks.get(key);
				
				// if there is at least one exp read,
				// xBlocks with mock reads will be discarded.
				if(isExpXBlock && xBlock.mockReadCount > 0) {
					discardList.add(key);
					continue;
				}
				
				if(xBlock.bestRegionPriority > minExpXBlockPenalty) {
					discardList.add(key);
				}
			}
			
			// remove higher penalties
			discardList.forEach(key -> {
				xBlocks.remove(key);
			});
		}
	}
}
