package progistar.pXg.processor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.data.GenomicSequence;
import progistar.pXg.data.Output;
import progistar.pXg.data.PeptideAnnotation;

public class Worker extends Thread {

	private Task task = null;
	private int workerID = 0;
	private long startTime = 0;
	private File tmpOutput = null;
	
	public Worker (int workerID, Task task) {
		super();
		this.task = task;
		this.workerID = workerID;
		this.startTime = System.currentTimeMillis();
		this.tmpOutput = new File(genTmpFilePath());
		
		System.out.println("Worker "+this.workerID+" takes task "+task.taskID);
	}
	
	private String genTmpFilePath () {
		return Parameters.sequenceFilePath + ".worker"+this.workerID+".tmp";
	}
	
	/**
	 * There are two types of tasks: <br>
	 * 
	 * 1) Mapping genomic sequence to genomic annotation. <br>
	 * 2) Make blast db. <br>
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
					if(matches.size() > 0) {
						Mapper.gMap(genomicSequence, task.gIndexStart, task.genomicAnnotationIndex, task.genomicAnnotation);
						
						/**
						 * Mapping output result to genomic annotation.
						 * 
						 */
						for(Output output : matches) {
							output.mapGenomicAnnotation();
						}
						
						BW.append(genomicSequence.uniqueID+"_"+genomicSequence.getLocus());
						BW.newLine();
						BW.append(genomicSequence.getNucleotideString());
						BW.newLine();
						
						for(int i=0; i<genomicSequence.matchedTxds; i++) {
							BW.append(genomicSequence.getGenomieRegion(i));
							BW.newLine();
						}
						
						for(Output output : matches) {
							String strand = null;
							if(output.strand) {
								strand = "forward";
							} else {
								strand = "reverse";
							}
							for(int i=0; i<genomicSequence.matchedTxds; i++) {
								BW.append(strand+"\t"+output.getPeptide()+"\t"+output.getLocus()+"\t"+output.getMatchedNucleotide()+"\t"+output.getGenomicRegion(i)+"\t"+output.getAARegionAnnotation(i));
								BW.newLine();
							}
						}
					}
				}
			}
			BW.close();
		}catch(IOException ioe) {
			
		}
		 
		
		long endTime = System.currentTimeMillis();
		System.out.println("Worker "+this.workerID+" was done with task"+task.taskID +"\tElapsed time: "+((endTime-startTime)/1000)+" sec");
	}
	
	
}
