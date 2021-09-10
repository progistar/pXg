package progistar.pXg.processor;

import progistar.pXg.constants.Constants;
import progistar.pXg.data.GenomicSequence;

public class Worker extends Thread {

	private Task task = null;
	private int workerID = 0;
	private long startTime = 0;
	
	public Worker (int workerID, Task task) {
		super();
		this.task = task;
		this.workerID = workerID;
		this.startTime = System.currentTimeMillis();
		
		System.out.println("Worker "+this.workerID+" takes task "+task.taskID);
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
		} 
		
		long endTime = System.currentTimeMillis();
		System.out.println("Worker "+this.workerID+" was done with task"+task.taskID +"\tElapsed time: "+((endTime-startTime)/1000)+" sec");
	}
	
	
}
