package progistar.pXg.processor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.constants.RunInfo;
import progistar.pXg.data.GenomicSequence;
import progistar.pXg.data.Mutation;
import progistar.pXg.data.Output;
import progistar.pXg.data.PeptideAnnotation;
import progistar.pXg.data.TBlock;
import progistar.pXg.data.parser.SamParser;

public class Worker extends Thread {

	private Task task = null;
	private int workerID = 0;
	private long startTime = 0;
	private File tmpOutput = null;
	
	public Worker (int workerID, Task task) {
		super();
		setTask(task);
		this.workerID = workerID;
		this.startTime = System.currentTimeMillis();
		this.tmpOutput = new File(genTmpFilePath());
		
		// enroll tmpOutput file path
		Master.enrollTmpOutputFilePath(this.tmpOutput.getAbsolutePath());
		
		//System.out.println("Worker "+this.workerID+" takes task "+task.taskID);
	}
	
	public void setTask (Task task) {
		this.task = task;
	}
	
	private String genTmpFilePath () {
		return Parameters.sequenceFilePaths[Parameters.CURRENT_FILE_INDEX] + 
				"." + Constants.UNIQUE_RUN_ID + ".worker" + this.workerID+ ".tmp";
	}
	
	/**
	 * 
	 * The decision of task is defined by Task.taskType <br>
	 * 
	 */
	public void run () {
		
		BufferedWriter BW = Master.getOutputBW(this.tmpOutput.getAbsolutePath());
		// task for mapping genomic annotation
		if(this.task.taskType == Constants.TASK_G_MAP) {
			for(String samRead : this.task.samReads) {
				RunInfo.workerProcessedReads[this.workerID] ++; // increase a number of processed reads
				
				GenomicSequence genomicSequence = SamParser.parseSam(samRead);
				ArrayList<Output> matches = PeptideAnnotation.find(genomicSequence);
				
				/**
				 * Only consider matched NGS-reads.
				 * Because we are not interested in unmatched NGS-reads.
				 * 
				 */
				
				if(matches.size() > 0) {
					Mapper.gMap(genomicSequence, task.gIndexStart, task.genomicAnnotationIndex, task.genomicAnnotation);
					
					/**
					 * Mapping output result to genomic annotation.
					 * 
					 */
					ArrayList<Output> targetMatches = new ArrayList<Output>();
					ArrayList<Output> decoyMatches = new ArrayList<Output>();
					
					for(Output output : matches) {
						output.mapGenomicAnnotation();
						
						if(output.isTarget) {
							targetMatches.add(output);
						} else {
							decoyMatches.add(output);
						}
					}
					
					if(targetMatches.size() != 0) {
						writeTmpOutput(BW, targetMatches, genomicSequence, "");
					}
					
					if(decoyMatches.size() != 0) {
						writeTmpOutput(BW, decoyMatches, genomicSequence, Constants.DECOY_PREFIX);
					}
				}
			}
			
			this.task.samReads.clear();
		}
		 
		
		long endTime = System.currentTimeMillis();
		System.out.println("Worker "+this.workerID+" was done with task"+task.taskID +"\telapsed time: "+((endTime-startTime)/1000)+" sec");
	}
	/**
	 * Write temporary output file. <br>
	 * 
	 * @param BW
	 * @param outputs
	 * @param gSeq
	 */
	public void writeTmpOutput (BufferedWriter BW, ArrayList<Output> outputs, GenomicSequence gSeq, String prefixID) {
		try {
			ArrayList<Byte> cdsTypes = new ArrayList<Byte>();
			cdsTypes.add(Constants.CDS);
			
			ArrayList<Byte> exonTypes = new ArrayList<Byte>();
			exonTypes.addAll(cdsTypes); exonTypes.add(Constants.UTR5); exonTypes.add(Constants.UTR3); exonTypes.add(Constants.NCDS);
			
			ArrayList<Byte> fullTypes = new ArrayList<Byte>();
			fullTypes.addAll(exonTypes); fullTypes.add(Constants.INTRON);
			
			BW.append(Constants.OUTPUT_G_UNIQUE_ID+"\t"+prefixID+gSeq.uniqueID+"@"+gSeq.getLocus());
			BW.newLine();
			

			// for unmapped/softclip we need more sequence information
			char mappingStatus = gSeq.getMappingStatus();
			if(mappingStatus == Constants.MARK_UNMAPPED || mappingStatus == Constants.MARK_SOFTCLIP) {
				BW.append(Constants.OUTPUT_G_SEQUENCE).append("\t").append(gSeq.getNucleotideString());
				BW.newLine();
			}
			
			// write average QScore.
			BW.append(Constants.OUTPUT_G_QSCORE+"\t"+gSeq.meanQScore);
			BW.newLine();
			
			for(Output output : outputs) {
				String strand = null;
				if(output.strand) {
					strand = "+";
				} else {
					strand = "-";
				}
				
				BW.append(Constants.OUTPUT_G_PEPTIDE);
				BW.append("\t");
				BW.append(output.getPeptide());
				BW.append("\t");
				BW.append(output.getLocus());
				BW.append("\t");
				BW.append(strand);
				BW.append("\t");
				BW.append(output.getLeftFlankNucleotide());
				BW.append("\t");
				BW.append(output.getMatchedNucleotide());
				BW.append("\t");
				BW.append(output.getRightFlankNucleotide());
				BW.append("\t");
				BW.append(output.getLeftFlankRefNucleotide());
				BW.append("\t");
				BW.append(output.getMatchedRefNucleotide());
				BW.append("\t");
				BW.append(output.getRightFlankRefNucleotide());
				BW.append("\t");
				
				// mutation check
				ArrayList<Mutation> mutations = output.getMutations();
				if(mutations.size() == 0) {
					BW.append("-");
				} else {
					for(int i=0; i<mutations.size(); i++) {
						if(i!=0) BW.append("|");
						BW.append(mutations.get(i).toString());
					}
				}
				BW.append("\t");
				// mutation status
				BW.append(output.getMutationStatus());
				BW.append("\t");
				
				// exon and full distance of identified genomic regions
				ArrayList<String> exonLengths = new ArrayList<String>();
				ArrayList<String> percentFullDistances = new ArrayList<String>();
				ArrayList<String> percentExonDistances = new ArrayList<String>();
				ArrayList<String> percentCDSDistances = new ArrayList<String>();
				ArrayList<String> fromStartDistances = new ArrayList<String>();
				ArrayList<String> fromStopDistances = new ArrayList<String>();
				
				for(int i=0; i<gSeq.matchedTxds; i++) {
					if(i!=0) BW.append("|");
					// intergenic
					String senseMarker = "-";
					
					String exonLength = "-";
					// percent distance for both exons and introns
					String percentFullDistance = "-";
					// percent distance for only exons
					String percentExonDistance = "-";
					// percent distance for only CDSs
					String percentCDSDistance = "-";
					String fromStartDistance = "-";
					String fromStopDistance = "-";
					
					mappingStatus = gSeq.getMappingStatus();
					
					if(gSeq.tBlocks[i] == null) {
						if(mappingStatus == Constants.MARK_MAPPED || mappingStatus == Constants.MARK_SOFTCLIP) {
							BW.append("intergenic");
						} else if (mappingStatus == Constants.MARK_UNMAPPED) {
							BW.append("unmapped");
						}
					} else {
						BW.append(gSeq.tBlocks[i].transcriptID);
						// check sense
						if(gSeq.tBlocks[i].strand == output.strand) {
							senseMarker = "sense";
						} else {
							senseMarker = "antisense";
						}
						
						int startPos = output.startGenomicPositions.get(0);
						int endPos = output.endGenomicPositions.get(output.endGenomicPositions.size()-1);
						int pos = output.strand ? startPos : endPos;
						
						
						// -1 means that there is no matched exons or introns
						percentFullDistance = getPercentDistance(pos, fullTypes, gSeq.tBlocks[i], output.strand);
						percentExonDistance = getPercentDistance(pos, exonTypes, gSeq.tBlocks[i], output.strand);
						
						if(senseMarker.equalsIgnoreCase("sense")) {
							percentCDSDistance = getPercentDistance(pos, cdsTypes, gSeq.tBlocks[i], output.strand);
							fromStartDistance = distFromStartSite(startPos, endPos, gSeq.tBlocks[i], output.strand);
							fromStopDistance = distFromStopSite(startPos, endPos, gSeq.tBlocks[i], output.strand);
						}
						
						int exonLen = gSeq.tBlocks[i].getTranscriptLength(exonTypes);
						exonLength = exonLen == -1 ? "-" : exonLen+"";
						
					}
					// add dist information
					exonLengths.add(exonLength);
					percentFullDistances.add(percentFullDistance);
					percentExonDistances.add(percentExonDistance);
					percentCDSDistances.add(percentCDSDistance);
					fromStartDistances.add(fromStartDistance);
					fromStopDistances.add(fromStopDistance);
					
					
					// if the match contains softclip, then the frame information is useless.
					// for softclip, a user must manually resolve the genomic origin of reliable identifications.
					// char frame = mappingStatus == Constants.MARK_SOFTCLIP ? Constants.NO_FRAME : output.getFrame(i);
					char frame = output.getFrame(i);
					char as = output.getAS(i); // alternative splicing mark
					BW.append("(").append(output.getAARegionAnnotation(i)).append(";")
					.append(senseMarker).append(";")
					.append(frame).append(";")
					.append(as).append(";")
					.append(")");
					
				}
				
				// distance information
				BW.append("\t");
				for(int i=0; i<exonLengths.size(); i++) {
					if(i != 0) BW.append("|");
					BW.append(exonLengths.get(i));
				}
				
				BW.append("\t");
				for(int i=0; i<percentFullDistances.size(); i++) {
					if(i != 0) BW.append("|");
					BW.append(percentFullDistances.get(i));
				}
				
				BW.append("\t");
				for(int i=0; i<percentExonDistances.size(); i++) {
					if(i != 0) BW.append("|");
					BW.append(percentExonDistances.get(i));
				}
				
				BW.append("\t");
				for(int i=0; i<percentCDSDistances.size(); i++) {
					if(i != 0) BW.append("|");
					BW.append(percentCDSDistances.get(i));
				}
				
				BW.append("\t");
				for(int i=0; i<fromStartDistances.size(); i++) {
					if(i != 0) BW.append("|");
					BW.append(fromStartDistances.get(i));
				}
				
				BW.append("\t");
				for(int i=0; i<fromStopDistances.size(); i++) {
					if(i != 0) BW.append("|");
					BW.append(fromStopDistances.get(i));
				}
				
				BW.newLine();
				
			}
		}catch(IOException ioe) {
			
		}
	}
	
	
	
	
	
	// #TODO: more proper Class for those methods.

	private String distFromStartSite (int start, int end, TBlock tBlock, boolean strand) {
		ArrayList<Byte> targets = new ArrayList<Byte>();
		targets.add(Constants.CDS); targets.add(Constants.UTR5); targets.add(Constants.UTR3);
		int peptStart = tBlock.getRelativeLengthOfPosition(start, targets);
		int peptEnd = tBlock.getRelativeLengthOfPosition(end, targets);
		int startSite = tBlock.getStartSite();
		
		// return "-"
		// when 1) intron, 2) noncoding, 3) intergenic, 4) antisense, 5) unmapped
		if(peptStart == -1 || peptEnd == -1 || startSite == -1 || (strand != tBlock.strand)) {
			return "-";
		}
		
		if(strand) {
			peptStart = peptStart - startSite;
			peptEnd = peptEnd - startSite;
		} else {
			peptStart = startSite - peptStart;
			peptEnd = startSite - peptEnd;
			int swap = peptStart;
			peptStart = peptEnd;
			peptEnd = swap;
		}
		

		if(peptStart >= 0) {
			peptStart++;
		}
		if(peptEnd >= 0) {
			peptEnd++;
		}
		
		String out = "";
		if(peptStart > 0) {
			out += "+";
		}
		out += peptStart +"~";
		if(peptEnd > 0) {
			out += "+";
		}
		out += peptEnd;
		
		return out;
	}
	
	private String distFromStopSite (int start, int end, TBlock tBlock, boolean strand) {
		ArrayList<Byte> targets = new ArrayList<Byte>();
		targets.add(Constants.CDS); targets.add(Constants.UTR5); targets.add(Constants.UTR3);
		int peptStart = tBlock.getRelativeLengthOfPosition(start, targets);
		int peptEnd = tBlock.getRelativeLengthOfPosition(end, targets);
		int endSite = tBlock.getStopSite();
		
		// return "-"
		// when 1) intron, 2) noncoding, 3) intergenic, 4) antisense, 5) unmapped
		if(peptStart == -1 || peptEnd == -1 || endSite == -1 || (strand != tBlock.strand)) {
			return "-";
		}
		
		if(strand) {
			peptStart = peptStart - endSite;
			peptEnd = peptEnd - endSite;
		} else {
			peptStart = endSite - peptStart;
			peptEnd = endSite - peptEnd;
			int swap = peptStart;
			peptStart = peptEnd;
			peptEnd = swap;
		}
		
		if(peptStart <= 0) {
			peptStart--;
		}
		if(peptEnd <= 0) {
			peptEnd--;
		}
		
		String out = "";
		if(peptStart > 0) {
			out += "+";
		}
		out += peptStart +"~";
		if(peptEnd > 0) {
			out += "+";
		}
		out += peptEnd;
		
		return out;
	}
	
	private String getPercentDistance (int pos, ArrayList<Byte> targets, TBlock tBlock, boolean strand) {
		int length = tBlock.getTranscriptLength(targets);
		int dist = tBlock.getRelativeLengthOfPosition(pos, targets);
		String percentDistance = "-";
		
		if(dist == -1 || length == -1) {
			percentDistance ="-";
		} else {
			if(!strand) {
				dist = length - dist + 1;
			}
			percentDistance = String.format("%.3f",(double)dist/length);
		}
		
		return percentDistance;
	}
}
