package progistar.pXg.processor;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.data.GenomicAnnotation;
import progistar.pXg.data.GenomicSequence;
import progistar.pXg.data.PxGAnnotation;
import progistar.pXg.data.parser.GTFParser;
import progistar.pXg.data.parser.PeptideParser;
import progistar.pXg.data.parser.ResultParser;
import progistar.pXg.data.parser.SamParser;
import progistar.pXg.mock.Mock;
import progistar.pXg.utils.Codon;

public class Master {

	private static GenomicAnnotation genomicAnnotation = null;
	private static int taskCount = 0;
	private static Hashtable<String, String> tmpOutputFilePaths = null;
	
	private Master() {}
	
	
	/**
	 * Load GTF and peptide file and ready to read SAM file <br>
	 * 
	 * @param genomicAnnotationFilePath
	 * @param sequenceFilePath
	 */
	public static void ready (String genomicAnnotationFilePath, String sequenceFilePath, String peptideFilePath) {
		// GTF parser and Peptide parser
		genomicAnnotation = GTFParser.parseGTF(genomicAnnotationFilePath);
		PeptideParser.parseResult(peptideFilePath); // static..!
		tmpOutputFilePaths = new Hashtable<String, String>();
				
		// TODO:
		// Make available to BAM file.
		SamParser.ready(sequenceFilePath);
		
		// loading Codon.
		Codon.mapping();
	}
	
	/**
	 * Start to map peptides to NGS-reads <br>
	 * 
	 */
	public static void run () {
		assert genomicAnnotation != null;
		// TODO:
		// Auto detection of already made index files.
		
		try {
			Worker[] workers = new Worker[Parameters.nThreads];
			
			ArrayList<GenomicSequence> genomicSequences = null;
			Vector<Task> taskQueue = new Vector<Task>(); // for synchronized
			Task[] tasks = null;
			
			while(!SamParser.isEndOfFile()) {
				genomicSequences = SamParser.parseSam(Parameters.readSize);
				// the array is initialized as "false"
				boolean[] assignedArray = new boolean[genomicSequences.size()];
				
				while(true) {
					tasks = getTasks(genomicSequences, assignedArray);
					boolean isEmpty = true;
					for(int i=0; i<tasks.length; i++) {
						if(tasks[i].isAssigned) {
							// add task into taskQueue
							taskQueue.add(tasks[i]);
							isEmpty = false;
						} 
					}
					
					while(!taskQueue.isEmpty()) {
						Task task = taskQueue.firstElement();
						taskQueue.remove(0);
						
						boolean isAssigned = false;
						
						while(!isAssigned) {
							for(int i=0; i<workers.length; i++) {
								if(workers[i] == null || !workers[i].isAlive()) {
									workers[i] = new Worker(i+1, task);
									workers[i].start();
									isAssigned = true;
									break;
								}
							}
							
							Thread.yield();
						}
					}
					
					if(isEmpty) break;
				}
			}
			
			SamParser.finish();
			
			// wait for finishing all tasks from workers
			waitUntilAllWorkersDone(workers);
			
			// read tmp output files
			ArrayList<File> tmpOutputFiles = new ArrayList<File>();
			tmpOutputFilePaths.forEach((path, value) ->{
				tmpOutputFiles.add(new File(path));
			});
			
			PxGAnnotation pXgA = ResultParser.parseResult(tmpOutputFiles);
			
			// removing tmpOutputFiles
//			tmpOutputFiles.forEach(file -> {file.delete();});
			
			// balancing mock read count
			pXgA.mockReadAssignPolicy();
			// marking target PSMs
			pXgA.markTargetPSM();
			// filter by pvalue
			pXgA.estimatePvalueThreshold();
			// among them, use highest-scored PSM
			pXgA.topScoreFilter();
			// fdr estimation
			pXgA.fdrEstimation();
			// filter regions
			pXgA.regionScoreFilter();
			// mark fasta result
			pXgA.markFasta();
			
			pXgA.write(Parameters.outputFilePath);
			
		}catch (Exception e) {
			e.printStackTrace();
		}
	}
	/**
	 * waiting for all workers are done with their tasks.<br>
	 * 
	 * @param workers
	 */
	private static void waitUntilAllWorkersDone (Worker[] workers) {
		boolean isProcessing = true;
		while(isProcessing) {
			isProcessing = false;
			for(int i=0; i<workers.length; i++) {
				isProcessing |= workers[i].isAlive();
			}
			Thread.yield();
		}
	}
	
	/**
	 * Partitioning tasks. <br>
	 * @param gSeqs
	 * @param assignedArray
	 * @return
	 */
	private static Task[] getTasks (ArrayList<GenomicSequence> gSeqs, boolean[] assignedArray) {
		
		assert gSeqs.size() != 0;
		
		
		Task[] tasks = new Task[Parameters.nThreads * 2];
		
		// Task class contains "isAssigned" feature.
		// The default value of the feature is "false"
		// The value becomes "true" when task has something to do.
		for(int i=0; i<tasks.length; i++) tasks[i] = new Task();
		
		int gSeqSize	=	gSeqs.size();
		int chrIndex	=	0;
		int start		=	0;
		int end			=	0;
		
		// Setting the pivot start position information
		for(int i=0; i<gSeqSize; i++) {
			// start position of NOT treated sequence
			if(!assignedArray[i]) {
				chrIndex	=	gSeqs.get(i).chrIndex;
				start		=	gSeqs.get(i).startPosition;
				end			=	start + Parameters.partitionSize - 1;
				
				break;
			}
		}
		
		ArrayList<GenomicSequence> gSeqPartitionIn = new ArrayList<GenomicSequence>();
		
		for(int i=0; i<gSeqSize; i++) {
			// already treated sequence
			if(assignedArray[i]) continue;
			
			GenomicSequence gSeq = gSeqs.get(i);
			
			if(gSeq.chrIndex == chrIndex) {
				if(gSeq.startPosition >= start && gSeq.endPosition <= end) {
					assignedArray[i] = true;
					gSeqPartitionIn.add(gSeq);
				}
			}
			
		}
		
		// have a task at least one.
		int partitionInSize = gSeqPartitionIn.size();
		if(partitionInSize != 0) {
			// Annotation index
			int[][] gIndex = genomicAnnotation.getIndexingBlocks(chrIndex, start, end);
			int taskIndex = 0;
			for(int i=0; i<partitionInSize; i++) {
				// target NGS-read
				tasks[taskIndex].genomicSequences.add(gSeqPartitionIn.get(i));
				// mock NGS-read
				if(Parameters.mocks != Constants.MOCK_NONE) {
					tasks[taskIndex].genomicSequences.add(Mock.makeMockRead(gSeqPartitionIn.get(i), Parameters.mocks));
				}
				
				taskIndex++;
				if(taskIndex == tasks.length) taskIndex = 0;
			}
			
			for(int i=0; i<tasks.length; i++) {
				if(tasks[i].genomicSequences.size() != 0) {
					tasks[i].isAssigned = true;
					tasks[i].genomicAnnotationIndex = gIndex;
					tasks[i].genomicAnnotation = genomicAnnotation;
					tasks[i].gIndexStart = start;
					tasks[i].taskID = ++Master.taskCount;
					tasks[i].taskType = Constants.TASK_G_MAP;
					
//					System.out.println(tasks[i].description());
				}
			}
		}
		
		return tasks;
	}
	
	/**
	 * Enroll temporary output file path. <br>
	 * The enrolled file paths will be processed when making the final output file. <br>
	 * 
	 * @param outputFilePath
	 */
	public static void enrollTmpOutputFilePath (String outputFilePath) {
		tmpOutputFilePaths.put(outputFilePath, "");
	}
}
