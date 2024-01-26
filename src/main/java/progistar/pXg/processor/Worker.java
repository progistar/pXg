package progistar.pXg.processor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.constants.RunInfo;
import progistar.pXg.data.GenomicSequence;
import progistar.pXg.data.Mutation;
import progistar.pXg.data.Output;
import progistar.pXg.data.PeptideAnnotation;
import progistar.pXg.data.parser.SamParser;
import progistar.pXg.utils.Logger;

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
		return Parameters.sequenceFilePath +"."+ Constants.UNIQUE_RUN_ID + ".worker" + this.workerID+ ".tmp";
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
		System.out.println("Worker "+this.workerID+" was done with task"+task.taskID +"\tElapsed time: "+((endTime-startTime)/1000)+" sec");
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
			
			/*
			 * These information are needed to confirm.
			 * Not informative for the last result.
			BW.append(Constants.OUTPUT_G_SEQUENCE+gSeq.getNucleotideString());
			BW.newLine();
			
			for(int i=0; i<gSeq.matchedTxds; i++) {
				BW.append(Constants.OUTPUT_G_REGION+gSeq.getGenomieRegion(i));
				BW.newLine();
			}
			*/
			
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
				
				int maxLenTxd = -1;
				double maxLen = -1;
				
				for(int i=0; i<gSeq.matchedTxds; i++) {
					if(i!=0) BW.append("|");
					// intergenic
					String senseMarker = "-";
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
						
						if(gSeq.tBlocks[i].getTranscriptFullLength() > maxLen) {
							maxLenTxd = i;
							maxLen = gSeq.tBlocks[i].getTranscriptFullLength();
						}
					}
					
					// if the match contains softclip, then the frame information is useless.
					// for softclip, a user must manually resolve the genomic origin of reliable identifications.
					// char frame = mappingStatus == Constants.MARK_SOFTCLIP ? Constants.NO_FRAME : output.getFrame(i);
					char frame = output.getFrame(i);
					char as = output.getAS(i); // alternative splicing mark
					BW.append("(").append(output.getAARegionAnnotation(i)).append(";").append(senseMarker).append(";").append(frame).append(";").append(as).append(")");
				}
				
				
				// distance proportion
				BW.append("\t");
				if(maxLenTxd == -1) {
					BW.append("-");
				} else {
					int pivotPos = -1;
					if(output.strand) {
						pivotPos = output.startGenomicPositions.get(0);
					} else {
						pivotPos = output.endGenomicPositions.get(output.endGenomicPositions.size()-1);
					}
					double dist = gSeq.tBlocks[maxLenTxd].getRelativeFullLengthRatio(pivotPos, output.strand);
					BW.append(dist+"");
				}
				
				BW.newLine();
				
			}
		}catch(IOException ioe) {
			
		}
	}
	
}
