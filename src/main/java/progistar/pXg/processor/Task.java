package progistar.pXg.processor;

import java.util.ArrayList;

import progistar.pXg.data.GenomicAnnotation;
import progistar.pXg.data.GenomicSequence;

public class Task {
	protected int[][] genomicAnnotationIndex;
	protected GenomicAnnotation genomicAnnotation;
	protected ArrayList<GenomicSequence> genomicSequences = new ArrayList<GenomicSequence>();
	protected int gIndexStart;
	protected boolean isAssigned = false;
	protected int taskID;
	protected int taskType;
	
	/**
	 * default values: <br>
	 * isAssigned = false <br>
	 * genomicSequences are allocated <br>
	 * 
	 */
	public Task() {}
	
	public Task(int[][] genomicAnnotationIndex, GenomicAnnotation genomicAnnotation,
			ArrayList<GenomicSequence> genomicSequences, int start, int taskID, int taskType) {
		super();
		this.genomicAnnotationIndex = genomicAnnotationIndex;
		this.genomicAnnotation = genomicAnnotation;
		this.genomicSequences = genomicSequences;
		this.gIndexStart = start;
		this.taskID = taskID;
		this.taskType = taskType;
	}
	
	public int getTaskType () {
		return this.taskType;
	}
	
	public int getTaskID () {
		return taskID;
	}
	
	public String description () {
		return "task " + this.taskID +": " + this.genomicSequences.size() +" sequences";
	}
}
