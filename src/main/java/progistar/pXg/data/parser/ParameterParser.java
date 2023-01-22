package progistar.pXg.data.parser;

import java.io.File;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.utils.Logger;

public class ParameterParser {

	/**
	 * If parsing successfully, return 0. else -1.
	 * 
	 * @param args
	 * @return
	 */
	public static int parseParams (String[] args) {
		try {
			System.out.println(Constants.VERSION+" "+Constants.RELEASE);
			System.out.println(Constants.INTRODUCE);
			System.out.println();
			
			// print parameter description
			if(args.length == 0) {
				System.out.println("Usage");
				System.out.println();
				System.out.println("Mandatory Fields");
				System.out.println("  -gtf                 : GTF file path. We recommand to use the same gtf corresponding to alignment.");
				System.out.println("  -sam                 : SAM file path. The sam file must be sorted by coordinate.");
				System.out.println("  -psm                 : PSM file path. It is expected that the psm file is derived from proteomics search by de novo or database search engine.");
				System.out.println("  -pept_col            : Peptide column index in the psm file. One-based!");
				System.out.println("  -scan_cols           : Scan identifier indices in the psm file. Multiple columns are also possible because sometimes it is not enough to distinguish scan by only scan id.");
				System.out.println("                         You can write multiple indices such like that: 1,2,5");
				System.out.println("  -out                 : Output path of pXg.");
				System.out.println();
				System.out.println("Optional Fields");
				System.out.println("  -sep                 : Specify the column separator. Possible values are csv or tsv. Default is csv");
				System.out.println("  -pval                : p-value cutoff of randomly matched peptide-read pairs. Default is 0.01");
				System.out.println("  -fdr                 : FDR cutoff to discard low-quality peptide-spectrum matches. Default is 0.1");
				System.out.println("  -ileq                : Controls whether pXg treats isoleucine (I) and leucine (L) as the same/equivalent with respect to a peptide identification. Default is true.");
				System.out.println("  -length              : Range of peptide length to consider. Default is 8-15");
				System.out.println("                         You can write in this way (min-max, both inclusive) : 8-13");
				System.out.println("  -fasta               : Canonical sequence database to report conservative assignment of noncanonical PSMs");
				System.out.println("  -rank                : How many candidates will be considered per a scan. Default is 100 (in other words, use all ranked candidates)");
				System.out.println("  -out_sam             : Report matched reads as SAM format (true or false). Default is true.");
				System.out.println("  -out_gtf             : Report matched peptides as GTF format (true or false). Default is true.");
				System.out.println("  -out_noncanonial     : Report noncaonical peptides for SAM and/or GTF formats (true or false). Default is true.");
				System.out.println("  -out_canonial        : Report caonical peptides for SAM and/or GTF formats (true or false). Default is true.");
				System.out.println("  -pMut                : Penalty per a mutation. Default is 1.");
				System.out.println("  -pAS                 : Penalty for alternative splicing. Default is 10.");
				System.out.println("  -p5UTR               : Penalty for 5`-UTR. Default is 20.");
				System.out.println("  -p3UTR               : Penalty for 3`-UTR. Default is 20.");
				System.out.println("  -pncRNA              : Penalty for noncoding RNA. Default is 20.");
				System.out.println("  -pFS                 : Penalty for frame shift. Default is 20.");
				System.out.println("  -pIR                 : Penalty for intron region. Default is 30.");
				System.out.println("  -pIGR                : Penalty for intergenic region. Default is 30.");
				System.out.println("  -pasRNA              : Penalty for antisense RNA. Default is 30.");
				System.out.println("  -punmap              : Penalty for unmapped reads. Default is 100.");
				System.out.println("  -gtf_partition_size  : The size of treating genomic region at once. Default is 5000000");
				System.out.println("  -sam_partition_size  : The size of treating number of reads at once. Default is 1000000");
				//System.out.println("  -threads             : The number of threads. Default is 4");
				System.out.println();
				System.out.println("Example1");
				System.out.println("java -Xmx30G -jar pXg.jar -gtf gencode.gtf -sam aligned.sorted.sam -psm peaks.result -pept_col 4 -score_col 8 -scan_cols 1,2,5  -pval 0.01 -fdr 0.1 -out_canonical false -out peaks.pXg");
				System.out.println("Example2");
				System.out.println("java -Xmx30G -jar pXg.jar -gtf gencode.gtf -sam aligned.sorted.sam -psm peaks.result -pept_col 4 -score_col 8 -scan_cols 1,2,5  -pval 0.01 -fdr 0.1 -length 8-13 -out peaks.pXg");
				return -1;
			}
			
			for(int i=0; i<args.length; i+=2) {
				String option = args[i].toLowerCase();
				
				// -gtf (mandatory)
				if(option.equalsIgnoreCase(Parameters.CMD_GENOMIC_ANNOTATION_PATH)) {
					Parameters.genomicAnnotationFilePath = args[i+1];
					if(!isExist(Parameters.genomicAnnotationFilePath)) {
						printNoSuchFileOrDirectory(Parameters.genomicAnnotationFilePath);
						return -1;
					}
				} 
				// -sam (mandatory)
				else if(option.equalsIgnoreCase(Parameters.CMD_GENOMIC_SEQUENCE_PATH)) {
					Parameters.sequenceFilePath = args[i+1];
					if(!isExist(Parameters.sequenceFilePath)) {
						printNoSuchFileOrDirectory(Parameters.sequenceFilePath);
						return -1;
					}
				}
				// -psm (mandatory)
				else if(option.equalsIgnoreCase(Parameters.CMD_PEPTIDE_ANNOTATION_PATH)) {
					Parameters.peptideFilePath = args[i+1];
					if(!isExist(Parameters.peptideFilePath)) {
						printNoSuchFileOrDirectory(Parameters.peptideFilePath);
						return -1;
					}
				}
				// -sep (optional)
				else if(option.equalsIgnoreCase(Parameters.SEP_TYPE)) {
					Parameters.sepType = args[i+1];
				}
				// -fasta (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_PROTEIN_SEQUENCE_PATH)) {
					Parameters.proteinFastaPath = args[i+1];
					if(!isExist(Parameters.proteinFastaPath)) {
						printNoSuchFileOrDirectory(Parameters.proteinFastaPath);
						return -1;
					}
				}
				// -out (mandatory)
				else if(option.equalsIgnoreCase(Parameters.CMD_OUTPUT_PATH)) {
					Parameters.outputFilePath = args[i+1];
					Parameters.ngsStatFilePath = Parameters.outputFilePath +".read.dist";
					Parameters.psmStatFilePath = Parameters.outputFilePath +".psm.dist";
					Parameters.unmappedFilePath = Parameters.outputFilePath +".unmapped";
					Parameters.exportGTFPath = Parameters.outputFilePath +".gtf";
					Parameters.exportSAMPath = Parameters.outputFilePath +".sam";
					
					if(isExist(Parameters.outputFilePath)) {
						printAlreadyExistFileOrDirectory(Parameters.outputFilePath);
						return -1;
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
				// -ileq (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_IL)) {
					if(args[i+1].equalsIgnoreCase("false")) {
						Parameters.leucineIsIsoleucine = false;
					}
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
				// -out_sam (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_SAM_FORMAT)) {
					if(args[i+1].equalsIgnoreCase("false")) {
						Parameters.EXPORT_SAM = false;
					}
				}
				// -out_gtf (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_GTF_FORMAT)) {
					if(args[i+1].equalsIgnoreCase("false")) {
						Parameters.EXPORT_GTF = false;
					}
				}
				// -out_canonical (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_CANONICAL)) {
					if(args[i+1].equalsIgnoreCase("false")) {
						Parameters.EXPORT_CANONICAL = false;
					}
				}
				// -out_noncanonical (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_NONCANONICAL)) {
					if(args[i+1].equalsIgnoreCase("false")) {
						Parameters.EXPORT_NONCANONICAL = false;
					}
				}
				// -pMut (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_PENALTY_MUTATION)) {
					Double penalty = Double.parseDouble(args[i+1]);
					Parameters.PENALTY_MUTATION = penalty;
				}
				// -pAS (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_PENALTY_AS)) {
					Double penalty = Double.parseDouble(args[i+1]);
					Parameters.PENALTY_AS = penalty;
				}
				// -p5UTR (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_PENALTY_5UTR)) {
					Double penalty = Double.parseDouble(args[i+1]);
					Parameters.PENALTY_5UTR = penalty;
				}
				// -p3UTR (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_PENALTY_3UTR)) {
					Double penalty = Double.parseDouble(args[i+1]);
					Parameters.PENALTY_3UTR = penalty;
				}
				// -pncRNA (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_PENALTY_ncRNA)) {
					Double penalty = Double.parseDouble(args[i+1]);
					Parameters.PENALTY_ncRNA= penalty;
				}
				// -pFS (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_PENALTY_FS)) {
					Double penalty = Double.parseDouble(args[i+1]);
					Parameters.PENALTY_FS = penalty;
				}
				// -pIR (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_PENALTY_IR)) {
					Double penalty = Double.parseDouble(args[i+1]);
					Parameters.PENALTY_IR= penalty;
				}
				// -pIGR (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_PENALTY_IGR)) {
					Double penalty = Double.parseDouble(args[i+1]);
					Parameters.PENALTY_IGR= penalty;
				}
				// -pasRNA (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_PENALTY_asRNA)) {
					Double penalty = Double.parseDouble(args[i+1]);
					Parameters.PENALTY_asRNA = penalty;
				}
				// -punmap (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_PENALTY_UNMAP)) {
					Double penalty = Double.parseDouble(args[i+1]);
					Parameters.PENALTY_UNMAP= penalty;
				}

			}
		}catch(Exception e) {
			System.out.println("Wrong parameter was detected. Please check the parameters.");
			return -1;
		}
		
		
		if(!isMandatoryOkay()) {
			return -1;
		}
		
		// open logger
    	Logger.create(Parameters.outputFilePath+".log");

		printSetting();
		
		// change one-based to zero-based
		Parameters.peptideColumnIndex--;
		Parameters.scoreColumnIndex--;
		for(int i=0; i<Parameters.scanColumnIndices.length; i++) {
			Parameters.scanColumnIndices[i]--;
		}
		
		return 0;
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
		System.out.println("  OUT_READ_DIST.: "+Parameters.ngsStatFilePath);
		System.out.println("  OUT_PSM_DIST.: "+Parameters.psmStatFilePath);
		System.out.println("  OUT_UNMAPPED: "+Parameters.unmappedFilePath);
		System.out.println("  OUT_SAM: "+Parameters.EXPORT_SAM);
		System.out.println("  OUT_GTF: "+Parameters.EXPORT_GTF);
		System.out.println("  OUT_CANONICAL: "+Parameters.EXPORT_CANONICAL);
		System.out.println("  OUT_NONCANONICAL: "+Parameters.EXPORT_NONCANONICAL);
		System.out.println(" pMut: "+Parameters.PENALTY_MUTATION);
		System.out.println(" pAS: "+Parameters.PENALTY_AS);
		System.out.println(" p5UTR: "+Parameters.PENALTY_5UTR);
		System.out.println(" p3UTR: "+Parameters.PENALTY_3UTR);
		System.out.println(" pncRNA: "+Parameters.PENALTY_ncRNA);
		System.out.println(" pFS: "+Parameters.PENALTY_FS);
		System.out.println(" pIR: "+Parameters.PENALTY_IR);
		System.out.println(" pIGR: "+Parameters.PENALTY_IGR);
		System.out.println(" pasRNA: "+Parameters.PENALTY_asRNA);
		System.out.println(" punmap: "+Parameters.PENALTY_UNMAP);
		//System.out.println(" THREADS: "+Parameters.nThreads);
		
		// append to logger
		Logger.append("Running info");
		Logger.newLine();
		Logger.append(" GTF: "+Parameters.genomicAnnotationFilePath);
		Logger.newLine();
		Logger.append("  GTF_PARTITION_SIZE: "+Parameters.partitionSize);
		Logger.newLine();
		Logger.append(" SAM: "+Parameters.sequenceFilePath);
		Logger.newLine();
		Logger.append("  SAM_PARTITION_SIZE: "+Parameters.readSize);
		Logger.newLine();
		Logger.append("  READ_CUTOFF_P_VALUE: "+Parameters.pvalue);
		Logger.newLine();
		Logger.append(" PSM: "+Parameters.peptideFilePath);
		Logger.newLine();
		Logger.append("  PEPT_COL: "+Parameters.peptideColumnIndex);
		Logger.newLine();
		Logger.append("  SCORE_COL: "+Parameters.scoreColumnIndex);
		Logger.newLine();
		Logger.append("  SCAN_COLS: "+scanCols.substring(1));
		Logger.newLine();
		Logger.append("  RANK TO CONSIDER: "+Parameters.psmRank);
		Logger.newLine();
		Logger.append("  PEPTIDE_LENGTH: "+Parameters.minPeptLen+"-"+Parameters.maxPeptLen);
		Logger.newLine();
		Logger.append("  PSM_FDR: "+Parameters.fdr);
		Logger.newLine();
		Logger.append(" OUT_RESULT: "+Parameters.outputFilePath);
		Logger.newLine();
		Logger.append("  OUT_READ_DIST.: "+Parameters.ngsStatFilePath);
		Logger.newLine();
		Logger.append("  OUT_PSM_DIST.: "+Parameters.psmStatFilePath);
		Logger.newLine();
		Logger.append("  OUT_UNMAPPED: "+Parameters.unmappedFilePath);
		Logger.newLine();
		Logger.append("  OUT_SAM: "+Parameters.EXPORT_SAM);
		Logger.newLine();
		Logger.append("  OUT_GTF: "+Parameters.EXPORT_GTF);
		Logger.newLine();
		Logger.append("  OUT_CANONICAL: "+Parameters.EXPORT_CANONICAL);
		Logger.newLine();
		Logger.append("  OUT_NONCANONICAL: "+Parameters.EXPORT_NONCANONICAL);
		Logger.newLine();
		Logger.append(" pMut: "+Parameters.PENALTY_MUTATION);
		Logger.newLine();
		Logger.append(" pAS: "+Parameters.PENALTY_AS);
		Logger.newLine();
		Logger.append(" p5UTR: "+Parameters.PENALTY_5UTR);
		Logger.newLine();
		Logger.append(" p3UTR: "+Parameters.PENALTY_3UTR);
		Logger.newLine();
		Logger.append(" pncRNA: "+Parameters.PENALTY_ncRNA);
		Logger.newLine();
		Logger.append(" pFS: "+Parameters.PENALTY_FS);
		Logger.newLine();
		Logger.append(" pIR: "+Parameters.PENALTY_IR);
		Logger.newLine();
		Logger.append(" pIGR: "+Parameters.PENALTY_IGR);
		Logger.newLine();
		Logger.append(" pasRNA: "+Parameters.PENALTY_asRNA);
		Logger.newLine();
		Logger.append(" punmap: "+Parameters.PENALTY_UNMAP);
		Logger.newLine();
		//Logger.append(" THREADS: "+Parameters.nThreads);
		//Logger.newLine();
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
