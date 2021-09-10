package progistar.pXg.processor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.data.GenomicSequence;
import progistar.pXg.data.PeptideAnnotation;
import progistar.pXg.tset.BLAST;

public class Worker extends Thread {

	private Task task = null;
	private int workerID = 0;
	private long startTime = 0;
	// Result of genomic-annotated NGS-read.
	// fastaFile for BLAST mapping
	// indexFile for genomic-annotation
	private File fastaFile = null;
	private File indexFile = null;
	
	public Worker (int workerID, Task task) {
		super();
		this.task = task;
		this.workerID = workerID;
		this.startTime = System.currentTimeMillis();
		this.fastaFile = new File("");
		
		this.fastaFile = new File(genFastaFilePath());
		this.indexFile = new File(genIndexFilePath());
		
		System.out.println("Worker "+this.workerID+" takes task "+task.taskID);
	}
	
	
	private String genFastaFilePath () {
		return Parameters.sequenceFilePath + ".worker"+this.workerID+".fasta";
	}
	
	private String genIndexFilePath () {
		return Parameters.sequenceFilePath + ".worker"+this.workerID+".idx";
	}
	
	private String getBlastResultPath () {
		return Parameters.sequenceFilePath + ".worker"+this.workerID+".blast.out";
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
		
		// task for mapping genomic annotation
		if(this.task.taskType == Constants.TASK_G_MAP) {
			for(GenomicSequence genomicSequence : this.task.genomicSequences) {
				Mapper.gMap(genomicSequence, task.gIndexStart, task.genomicAnnotationIndex, task.genomicAnnotation);
			}
			
			// write down fasta/index files.
			try {
				BufferedWriter fastaBW = new BufferedWriter(new FileWriter(this.fastaFile, true));
				BufferedWriter indexBW = new BufferedWriter(new FileWriter(this.indexFile, true));
				for(GenomicSequence genomicSequence : this.task.genomicSequences) {
					
					// skip unmapped reads
					if(genomicSequence.getNucleotideString().length() == 0) continue;
					
					fastaBW.append(genomicSequence.getFastaHeader());
					fastaBW.newLine();
					fastaBW.append(genomicSequence.getNucleotideString());
					fastaBW.newLine();
					
					// if you want to change the index table,
					// just fix the GenomicSequence.getGenomicFeature
					indexBW.append(genomicSequence.getGenomicFeature());
				}
				fastaBW.close();
				indexBW.close();
			}catch(IOException ioe) {
				
			}
		} 
		// task for making blast DB
		else if(this.task.taskType == Constants.TASK_MAKE_BLAST_DB) {
			if(this.fastaFile != null && this.fastaFile.exists()) {
				BLAST.makeBlastDB(this.fastaFile.getAbsolutePath(), "nucl");
			}
		}
		// task for mapping blast DB
		else if(this.task.taskType == Constants.TASK_MAP_BLAST) {
			if(this.fastaFile != null && this.fastaFile.exists()) {
				BLAST.tblastn(PeptideAnnotation.getFastaQueryFilePath(), this.fastaFile.getAbsolutePath(), this.getBlastResultPath(), false);
			}
		}
		
		long endTime = System.currentTimeMillis();
		System.out.println("Worker "+this.workerID+" was done with task"+task.taskID +"\tElapsed time: "+((endTime-startTime)/1000)+" sec");
	}
	
	
}
