package progistar.pXg.data.parser;

import java.io.File;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;

public class ParameterParser {

	public static void parseParams (String[] args) {
		try {
			System.out.println(Constants.VERSION);
			System.out.println(Constants.INTRODUCE);
			System.out.println();
			
			// print parameter description
			if(args.length == 0) {
				System.out.println("Usage");
				System.out.println();
				System.out.println("Mandatory Fields");
				System.out.println("  -gtf                 : gtf file path. We recommand to use the same gtf corresponding to alignment.");
				System.out.println("  -sam                 : sam file path. The sam file must be sorted by coordinate.");
				System.out.println("  -psm                 : psm file path. It is expected that the psm file is derived from proteomics search by de novo or database search engine.");
				System.out.println("  -pept_col            : Peptide column index in the psm file. One-based!");
				System.out.println("  -scan_cols           : Scan identifier indices in the psm file. Multiple columns are also possible because sometimes it is not enough to distinguish scan by only scan id.");
				System.out.println("                         You can write multiple indices such like that: 1,2,5");
				System.out.println("  -pept_col            : Peptide column index in the psm file. One-based!");
				System.out.println("  -out                 : Output path of pXg.");
				System.out.println();
				System.out.println("Optional Fields");
				System.out.println("  -pval                : p-value cutoff of randomly matched peptide-read pairs. Default is 0.05");
				System.out.println("  -fdr                 : fdr cutoff to discard low-quality peptide-spectrum matches. Default is 0.05");
				System.out.println("  -length              : Range of peptide length to consider. Default is 8-15");
				System.out.println("  -rank                : How many candidates will be considered per a scan. Default is 10");
				System.out.println("                         You can write in this way (min-max, both inclusive) : 8-13");
				System.out.println("  -gtf_partition_size  : The size of treating genomic region at once. Default is 5000000");
				System.out.println("  -sam_partition_size  : The size of treating number of reads at once. Default is 1000000");
				System.out.println("  -threads             : The number of threads. Default is 4");
				System.out.println();
				System.out.println("Example1");
				System.out.println("java -Xmx30G -jar pXg.jar -gtf gencode.gtf -sam aligned.sorted.sam -psm peaks.result -pept_col 4 -score_col 8 -scan_cols 1,2,5 -out peaks.pXg");
				System.out.println("Example2");
				System.out.println("java -Xmx30G -jar pXg.jar -gtf gencode.gtf -sam aligned.sorted.sam -psm peaks.result -pept_col 4 -score_col 8 -scan_cols 1,2,5 -out peaks.pXg -pval 0.05 -fdr 0.01 -length 8-13 -threads 2");
				System.exit(1);
			}
			
			for(int i=0; i<args.length; i+=2) {
				String option = args[i].toLowerCase();
				
				// -gtf (mandatory)
				if(option.equalsIgnoreCase(Parameters.CMD_GENOMIC_ANNOTATION_PATH)) {
					Parameters.genomicAnnotationFilePath = args[i+1];
					if(!isExist(Parameters.genomicAnnotationFilePath)) {
						printNoSuchFileOrDirectory(Parameters.genomicAnnotationFilePath);
						System.exit(1);
					}
				} 
				// -sam (mandatory)
				else if(option.equalsIgnoreCase(Parameters.CMD_GENOMIC_SEQUENCE_PATH)) {
					Parameters.sequenceFilePath = args[i+1];
					if(!isExist(Parameters.sequenceFilePath)) {
						printNoSuchFileOrDirectory(Parameters.sequenceFilePath);
						System.exit(1);
					}
				}
				// -psm (mandatory)
				else if(option.equalsIgnoreCase(Parameters.CMD_PEPTIDE_ANNOTATION_PATH)) {
					Parameters.peptideFilePath = args[i+1];
					if(!isExist(Parameters.peptideFilePath)) {
						printNoSuchFileOrDirectory(Parameters.peptideFilePath);
						System.exit(1);
					}
				}
				// -fasta (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_PROTEIN_SEQUENCE_PATH)) {
					Parameters.proteinFastaPath = args[i+1];
					if(!isExist(Parameters.proteinFastaPath)) {
						printNoSuchFileOrDirectory(Parameters.proteinFastaPath);
						System.exit(1);
					}
				}
				// -out (mandatory)
				else if(option.equalsIgnoreCase(Parameters.CMD_OUTPUT_PATH)) {
					Parameters.outputFilePath = args[i+1];
					Parameters.ngsStatFilePath = Parameters.outputFilePath +".pval.dist";
					Parameters.psmStatFilePath = Parameters.outputFilePath +".fdr.dist";
					Parameters.unmappedFilePath = Parameters.outputFilePath +".unmapped";
					
					if(isExist(Parameters.outputFilePath)) {
						printAlreadyExistFileOrDirectory(Parameters.outputFilePath);
						System.exit(1);
					}
				}
				// -pval (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_P_VALUE)) {
					Parameters.pvalue = Double.parseDouble(args[i+1]);
				}
				// -fdr (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_FDR)) {
					Parameters.fdr = Double.parseDouble(args[i+1]);
				}
				// -pept_col (mandatory)
				else if(option.equalsIgnoreCase(Parameters.CMD_PEPTIDE_COLUMN_INDEX)) {
					Parameters.peptideColumnIndex = Integer.parseInt(args[i+1]);
				}
				// -scan_cols (mandatory)
				else if(option.equalsIgnoreCase(Parameters.CMD_SCAN_COLUMN_INDICES)) {
					String[] indicies = args[i+1].split("\\,");
					Parameters.scanColumnIndices = new int[indicies.length];
					for(int idx=0; idx<indicies.length; idx++) {
						Parameters.scanColumnIndices[idx] = Integer.parseInt(indicies[idx]);
					}
				}
				// -score_col (mandatory)
				else if(option.equalsIgnoreCase(Parameters.CMD_SCORE_COLUMN_INDEX)) {
					Parameters.scoreColumnIndex = Integer.parseInt(args[i+1]);
				}
				// -rank (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_CANDIDATE_RANK)) {
					Parameters.psmRank = Integer.parseInt(args[i+1]);
				}
				// -length (optional)
				// >example: 8-15
				else if(option.equalsIgnoreCase(Parameters.CMD_LENGTH)) {
					String[] range = args[i+1].split("\\-");
					Parameters.minPeptLen = Integer.parseInt(range[0]);
					Parameters.maxPeptLen = Integer.parseInt(range[1]);
				}
				// -gtf_partition_size (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_GENOMIC_ANNOTATION_PARTITION_SIZE)) {
					Parameters.partitionSize = Integer.parseInt(args[i+1]);
				}
				// -sam_partition_size (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_GENOMIC_SEQUENCE_PARTITION_SIZE)) {
					Parameters.readSize = Integer.parseInt(args[i+1]);
				}
				// -sam_partition_size (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_THREADS)) {
					Parameters.nThreads = Integer.parseInt(args[i+1]);
				}

			}
		}catch(Exception e) {
			System.out.println("Wrong parameter was detected. Please check the parameters.");
			System.exit(1);
		}
		
		
		if(!isMandatoryOkay()) {
			System.exit(1);
		}
		
		printSetting();
	}
	
	/**
	 * GTF
	 * SAM
	 * PSM
	 * FASTA
	 * OUT
	 * P-value
	 * FDR
	 * SCAN_COLS
	 * PEPT_COL
	 * SCORE_COL
	 * RANKS TO CONSIDER
	 * LENGTH
	 * GTF_PARTITION_SIZE
	 * SAM_PARTITION_SIZE
	 * THREADS
	 */
	private static void printSetting () {
		System.out.println("Running info");
		System.out.println(" GTF: "+Parameters.genomicAnnotationFilePath);
		System.out.println("  GTF_PARTITION_SIZE: "+Parameters.partitionSize);
		System.out.println(" SAM: "+Parameters.sequenceFilePath);
		System.out.println("  SAM_PARTITION_SIZE: "+Parameters.readSize);
		System.out.println("  READ_CUTOFF_P_VALUE: "+Parameters.pvalue);
		System.out.println(" PSM: "+Parameters.peptideFilePath);
		System.out.println("  PEPT_COL: "+Parameters.peptideColumnIndex);
		System.out.println("  SCORE_COL: "+Parameters.scoreColumnIndex);
		
		// to display array
		String scanCols = "";
		for(int i=0; i<Parameters.scanColumnIndices.length; i++) {
			scanCols += "," + Parameters.scanColumnIndices[i];
		}
		System.out.println("  SCAN_COLS: "+scanCols.substring(1));
		System.out.println("  RANK TO CONSIDER: "+Parameters.psmRank);
		System.out.println("  PEPTIDE_LENGTH: "+Parameters.minPeptLen+"-"+Parameters.maxPeptLen);
		System.out.println("  PSM_FDR: "+Parameters.fdr);
		System.out.println(" OUT_RESULT: "+Parameters.outputFilePath);
		System.out.println(" OUT_READ_DIST.: "+Parameters.ngsStatFilePath);
		System.out.println(" OUT_PSM_DIST.: "+Parameters.psmStatFilePath);
		System.out.println(" OUT_UNMAPPED: "+Parameters.unmappedFilePath);
		System.out.println(" THREADS: "+Parameters.nThreads);
	}
	
	private static boolean isMandatoryOkay () {
		boolean pass = true;
		
		// -gtf
		if(Parameters.genomicAnnotationFilePath == null) {
			System.out.println("mandatory option -gtf is blank...");
			pass = false;
		}
		// -sam
		if(Parameters.sequenceFilePath == null) {
			System.out.println("mandatory option -sam is blank...");
			pass = false;
		}
		// -psm
		if(Parameters.peptideFilePath == null) {
			System.out.println("mandatory option -psm is blank...");
			pass = false;
		}
		// -out
		if(Parameters.outputFilePath == null) {
			System.out.println("mandatory option -out is blank...");
			pass = false;
		}
		// -pept_col
		if(Parameters.peptideColumnIndex == -1) {
			System.out.println("mandatory option -pept_col is blank...");
			pass = false;
		}
		// -score_col
		if(Parameters.scoreColumnIndex == -1) {
			System.out.println("mandatory option -score_col is blank...");
			pass = false;
		}
		// -sacn_cols
		if(Parameters.scanColumnIndices == null) {
			System.out.println("mandatory option -scan_cols is blank...");
			pass = false;
		}
		
		return pass;
	}
	
	private static void printNoSuchFileOrDirectory (String fileName) {
		System.out.println("No such file or directory: "+fileName);
	}
	
	private static void printAlreadyExistFileOrDirectory (String fileName) {
		System.out.println("The file already exists: "+fileName);
	}
	
	private static boolean isExist (String fileName) {
		File file = new File(fileName);
		return file.exists();
	}
}
