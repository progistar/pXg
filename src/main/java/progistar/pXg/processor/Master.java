package progistar.pXg.processor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Vector;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.constants.RunInfo;
import progistar.pXg.data.GenomicAnnotation;
import progistar.pXg.data.PIN;
import progistar.pXg.data.PeptideAnnotation;
import progistar.pXg.data.PxGAnnotation;
import progistar.pXg.data.parser.GTFParser;
import progistar.pXg.data.parser.PeptideParser;
import progistar.pXg.data.parser.ResultParser;
import progistar.pXg.data.parser.SamParser;
import progistar.pXg.utils.Codon;
import progistar.pXg.utils.IndexConvertor;

public class Master {

	private static GenomicAnnotation genomicAnnotation = null;
	private static int taskCount = 0;
	private static Hashtable<String, BufferedWriter> tmpOutputFilePaths = null;
	private static File SAM_FILE = null;
	private static int[] chrIndices = null;
	private static int[] startPositions = null;
	private static boolean[] assignedArray = null;
	private static String[] reads = null;
	
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
		tmpOutputFilePaths = new Hashtable<String, BufferedWriter>();
				
		// TODO:
		// Make available to BAM file.
		SAM_FILE = new File(sequenceFilePath);
		
		// for SAM - GTF associated task assignment
		// TASK-related variables
		chrIndices = new int[Parameters.readSize];
		startPositions = new int[Parameters.readSize];
		assignedArray = new boolean[Parameters.readSize];
		reads = new String[Parameters.readSize];
		
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
			RunInfo.workerProcessedReads = new long[Parameters.nThreads+1];
			Worker[] workers = new Worker[Parameters.nThreads];
			
			try (SamReader samReader = SamReaderFactory.makeDefault().open(SAM_FILE)) {
				
				// initialize
				Vector<Task> taskQueue = new Vector<Task>(); // for synchronized
				String lineSeparator = System.getProperty("line.separator");
				int readCount = 0;
				long readPartitionSize = Parameters.readSize;
				
	        	// get records				
	            for (SAMRecord samRecord : samReader) {
	                // Process each SAMRecord as needed
	            	String record = samRecord.getSAMString().replace(lineSeparator, "");
	            	reads[readCount] =record;
	            	
	            	String[] fields = record.split("\\s");
	            	String chr = fields[SamParser.CHR_IDX];
	            	Integer startPosition = Integer.parseInt(fields[SamParser.START_POS_IDX]);
					
					// the index for that chr is automatically assigned by auto-increment key.
					IndexConvertor.putChrIndexer(chr);
					int chrIndex_ = IndexConvertor.chrToIndex(chr);
					// check all chromosomes are well processed.
					RunInfo.processedChromosomes.put(chr, chrIndex_);
					
					// store chr and start positions
					chrIndices[readCount] = chrIndex_;
					startPositions[readCount] = startPosition;
					readCount ++;
					if(readCount == readPartitionSize) {
						// the array is initialized as "false"
						Arrays.fill(assignedArray, false);
						// assign tasks to workers
						assignTasks(taskQueue, workers, readCount);
						// once give the tasks, remove reads
						RunInfo.totalProcessedReads += readCount;
						readCount = 0;
					}
	            }
	            
	            // do last tasks
	            if(readCount > 0) {
	            	Arrays.fill(assignedArray, false);
					// assign tasks to workers
					assignTasks(taskQueue, workers, readCount);
					// once give the tasks, remove reads
					RunInfo.totalProcessedReads += readCount;
					readCount = 0;
	            }
	            
	            
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
				
			// wait for finishing all tasks from workers
			waitUntilAllWorkersDone(workers);
			
			// read tmp output files
			ArrayList<File> tmpOutputFiles = new ArrayList<File>();
			tmpOutputFilePaths.forEach((path, tmpBW) ->{
				tmpOutputFiles.add(new File(path));
				Master.closeOutputBW(path);
			});
			
			PxGAnnotation pXgA = ResultParser.parseResult(tmpOutputFiles);
			
			// removing tmpOutputFiles
			tmpOutputFiles.forEach(file -> {file.delete();});
			
			// count peptides and scans matching to exp.reads
			RunInfo.mappingFilterPeptideNum3 = PeptideAnnotation.getPeptideSizeWithXBlocks(pXgA.getXBlockMapper());
			RunInfo.mappingFilterScanNum3 = PeptideAnnotation.getScanSizeWithXBlocks(pXgA.getXBlockMapper());
			
			// count peptides and scans after p-value
			RunInfo.pvalueFilterPeptideNum4 = PeptideAnnotation.getPeptideSizeWithXBlocks(pXgA.getXBlockMapper());
			RunInfo.pvalueFilterScanNum4 = PeptideAnnotation.getScanSizeWithXBlocks(pXgA.getXBlockMapper());
			
			// filter regions
			pXgA.regionScoreFilter();
			// mark fasta result
			// to distinguish ambiguous interpretation 
			pXgA.markFasta();
			// marking target PSMs
			pXgA.assignXBlocks();
			
			pXgA.write(Parameters.outputFilePath);
			
			// parser to PIN
			PIN.parseOutput();
		}catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void assignTasks (Vector<Task> taskQueue, Worker[] workers, int readCount) {
		Task[] tasks = getTasks(readCount);
		for(int i=0; i<tasks.length; i++) {
			if(tasks[i].isAssigned) {
				// add task into taskQueue
				taskQueue.add(tasks[i]);
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
				if(workers[i] != null)
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
	private static Task[] getTasks (int readCount) {
		
		assert readCount != 0;
		
		
		Task[] tasks = new Task[Parameters.nThreads];
		
		// Task class contains "isAssigned" feature.
		// The default value of the feature is "false"
		// The value becomes "true" when task has something to do.
		for(int i=0; i<tasks.length; i++) tasks[i] = new Task();
		
		int gSeqSize	=	readCount;
		int chrIndex	=	0;
		int start		=	0;
		int end			=	0;
		
		// Setting the pivot start position information
		for(int i=0; i<gSeqSize; i++) {
			// start position of NOT treated sequence
			if(!assignedArray[i]) {
				chrIndex	=	chrIndices[i];
				start		=	startPositions[i];
				end			=	start + Parameters.partitionSize - 1;
				
				break;
			}
		}
		
		ArrayList<String> readPartition = new ArrayList<String>();
		
		for(int i=0; i<gSeqSize; i++) {
			// already treated sequence
			if(assignedArray[i]) continue;
			
			if(chrIndices[i] == chrIndex) {
				if(startPositions[i] >= start && startPositions[i] + Parameters.maxJunctionSize <= end) {
					assignedArray[i] = true;
					readPartition.add(reads[i]);
				}
			}
			
		}
		
		// have a task at least one.
		int partitionInSize = readPartition.size();
		if(partitionInSize != 0) {
			// Annotation index
			int[][] gIndex = genomicAnnotation.getIndexingBlocks(chrIndex, start, end);
			int taskIndex = 0;
			for(int i=0; i<partitionInSize; i++) {
				// target NGS-read
				tasks[taskIndex].samReads.add(readPartition.get(i));
				taskIndex++;
				if(taskIndex == tasks.length) taskIndex = 0;
			}
			
			for(int i=0; i<tasks.length; i++) {
				if(tasks[i].samReads.size() != 0) {
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
		if(tmpOutputFilePaths.get(outputFilePath) == null) {
			try {
				tmpOutputFilePaths.put(outputFilePath, new BufferedWriter(new FileWriter(outputFilePath)));
			}catch(IOException ioe) {
				
			}
		}
	}
	
	public static BufferedWriter getOutputBW (String outputFilePath) {
		BufferedWriter BW = tmpOutputFilePaths.get(outputFilePath);
		if(BW == null) {
			enrollTmpOutputFilePath(outputFilePath);
			BW = tmpOutputFilePaths.get(outputFilePath);
		}
		
		return BW;
	}
	
	public static void closeOutputBW (String outputFilePath) {
		BufferedWriter BW = tmpOutputFilePaths.get(outputFilePath);
		if(BW != null) {
			try {
				BW.close();
			}catch(IOException ioe) {
				
			}
		}
	}
}
