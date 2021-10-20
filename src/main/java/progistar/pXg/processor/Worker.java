package progistar.pXg.processor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.data.GenomicSequence;
import progistar.pXg.data.Mutation;
import progistar.pXg.data.Output;
import progistar.pXg.data.PeptideAnnotation;

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
		return Parameters.sequenceFilePath + ".worker"+this.workerID+".tmp";
	}
	
	/**
	 * 
	 * The decision of task is defined by Task.taskType <br>
	 * 
	 */
	public void run () {
		
		try {
			BufferedWriter BW = new BufferedWriter(new FileWriter(this.tmpOutput, true));
			// task for mapping genomic annotation
			if(this.task.taskType == Constants.TASK_G_MAP) {
				for(GenomicSequence genomicSequence : this.task.genomicSequences) {
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
						for(Output output : matches) {
							output.mapGenomicAnnotation();
						}
						
						writeTmpOutput(BW, matches, genomicSequence);
					}
				}
				
				this.task.genomicSequences.clear();
			}
			BW.close();
		}catch(IOException ioe) {
			
		}
		 
		
		long endTime = System.currentTimeMillis();
		//System.out.println("Worker "+this.workerID+" was done with task"+task.taskID +"\tElapsed time: "+((endTime-startTime)/1000)+" sec");
	}
	/**
	 * Write temporary output file. <br>
	 * 
	 * @param BW
	 * @param outputs
	 * @param gSeq
	 */
	public void writeTmpOutput (BufferedWriter BW, ArrayList<Output> outputs, GenomicSequence gSeq) {
		try {

			BW.append(Constants.OUTPUT_G_UNIQUE_ID+"\t"+gSeq.uniqueID+"_"+gSeq.getLocus());
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
				BW.append(output.getMatchedNucleotide());
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
				
				for(int i=0; i<gSeq.matchedTxds; i++) {
					if(i!=0) BW.append("|");
					// intergenic
					String senseMarker = "-";
					if(gSeq.tBlocks[i] == null) {
						BW.append("intergenic");
					} else {
						BW.append(gSeq.tBlocks[i].transcriptID);
						// check sense
						if(gSeq.tBlocks[i].strand == output.strand) {
							senseMarker = "sense";
						} else {
							senseMarker = "antisense";
						}
					}
					char frame = output.getFrame(i);
					BW.append("(").append(output.getAARegionAnnotation(i)).append(";").append(senseMarker).append(";").append(frame).append(")");
				}
				BW.newLine();
			}
		}catch(IOException ioe) {
			
		}
	}
	
}
