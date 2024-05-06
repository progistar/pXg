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
				System.out.println("  --gtf_file            : GTF file path. We recommand to use the same gtf corresponding to alignment.");
				System.out.println("  --sam_file            : SAM/BAM file path. The sam/bam file must be sorted by coordinate.");
				System.out.println("                          Multiple sam/bam files should be separated by comma (,).");
				System.out.println("  --psm_file            : PSM file path. It is expected that the psm file is derived from proteomics search by de novo or database search engine.");
				System.out.println("  --pept_col            : Peptide column index in the psm file. One-based!");
				System.out.println("  --scan_col            : Scan number index (the value is expected as integer > 0). One-based!");
				System.out.println("  --score_col           : Main score index. One-based!");
				System.out.println("  --file_col            : File name index. One-based!");
				System.out.println("  --charge_col          : Charge state index. One-based!");
				System.out.println("  --output              : Output file name of pXg.");
				System.out.println();
				System.out.println("Optional Fields");
				System.out.println("  --add_feat_cols       : Specify the indices for additional features to generate PIN file. One-based!");
				System.out.println("                          Several features can be added by comma separator. ex> 5,6,7");
				System.out.println("  --sep                 : Specify the column separator. Possible values are csv or tsv. Default is tsv");
				System.out.println("  --mode                : Specify the method of translation nucleotides. 3 for three-frame and 6 for six-frame. Default is 3");
				System.out.println("  --ileq                : Controls whether pXg treats isoleucine (I) and leucine (L) as the same/equivalent with respect to a peptide identification. Default is true.");
				System.out.println("  --lengths             : Range of peptide length to consider. Default is 8-15");
				System.out.println("                          You can write in this way (min-max, both inclusive) : 8-13");
				System.out.println("  --max_flank_size      : Specify to print maximum flank nuleotides from matched sequence. Default is 1,000");
				System.out.println("  --fasta_file          : Canonical sequence database to avoid ambiguous assignment of noncanonical peptides");
				System.out.println("  --rank                : How many candidates will be considered per a scan. Default is 100 (in other words, use all ranked candidates)");
				System.out.println("  --out_sam             : Report matched reads as SAM format (true or false). Default is false.");
				System.out.println("  --out_noncanonial     : Report noncaonical peptides for SAM and/or GTF formats (true or false). Default is true.");
				System.out.println("  --out_canonial        : Report caonical peptides for SAM and/or GTF formats (true or false). Default is true.");
				System.out.println("  --penalty_mutation    : Penalty per a mutation. Default is 1.");
				System.out.println("  --penalty_AS          : Penalty for alternative splicing. Default is 10.");
				System.out.println("  --penalty_5UTR        : Penalty for 5`-UTR. Default is 20.");
				System.out.println("  --penalty_3UTR        : Penalty for 3`-UTR. Default is 20.");
				System.out.println("  --penalty_ncRNA       : Penalty for noncoding RNA. Default is 20.");
				System.out.println("  --penalty_FS          : Penalty for frame shift. Default is 20.");
				System.out.println("  --penalty_IR          : Penalty for intron region. Default is 30.");
				System.out.println("  --penalty_IGR         : Penalty for intergenic region. Default is 30.");
				System.out.println("  --penalty_asRNA       : Penalty for antisense RNA. Default is 30.");
				System.out.println("  --penalty_softclip    : Penalty for softclip reads. Default is 50.");
				System.out.println("  --penalty_unknown     : Penalty for unmapped reads. Default is 100.");
				System.out.println("  --gtf_partition_size  : The size of treating genomic region at once. Default is 10000000");
				System.out.println("  --sam_partition_size  : The size of treating number of reads at once. Default is 1000000");
				System.out.println("  --threads             : The number of threads. Default is 4");
				System.out.println();
				System.out.println("Example1");
				System.out.println("java -Xmx30G -jar pXg.jar -gtf_file gencode.gtf -sam_file aligned.sorted.bam -psm_file peaks.result -scan_col 5 -file_col 2 -pept_col 4 -charge_col 11 -score_col 8 -add_feat_cols 14,15  -out_canonical false -out test");
				System.out.println("Example2");
				System.out.println("java -Xmx30G -jar pXg.jar -gtf_file gencode.gtf -sam_file aligned.sorted.sam -psm_file peaks.result -scan_col 5 -file_col 2 -pept_col 4 -charge_col 11 -score_col 8 -add_feat_cols 14,15 -lengths 8-11 -out test");
				return -1;
			}

			// sam file first
			for(int i=0; i<args.length; i+=2) {
				String option = args[i].toLowerCase();
				// --sam_file (mandatory)
				if(option.equalsIgnoreCase(Parameters.CMD_GENOMIC_SEQUENCE_PATH)) {
					String[] paths = args[i+1].split(",");
					Parameters.NUM_OF_SAM_FILES = paths.length;
					Parameters.sequenceFilePaths = new String[Parameters.NUM_OF_SAM_FILES];
					Parameters.unmappedFilePaths = new String[Parameters.NUM_OF_SAM_FILES];
					Parameters.exportSAMPaths	 = new String[Parameters.NUM_OF_SAM_FILES];
					Parameters.tmpOutputFilePaths= new String[Parameters.NUM_OF_SAM_FILES];

					for(int idx=0; idx < Parameters.NUM_OF_SAM_FILES; idx++) {
						Parameters.sequenceFilePaths[idx] = paths[idx];
						if(!(Parameters.sequenceFilePaths[idx].toLowerCase().endsWith(".bam") ||
								Parameters.sequenceFilePaths[idx].toLowerCase().endsWith(".sam"))) {
							System.out.println(Parameters.sequenceFilePaths[idx] +" must be .sam or .bam file.");
							return -1;
						}

						int lastIdx = Parameters.sequenceFilePaths[idx].lastIndexOf(".");

						Parameters.unmappedFilePaths[idx] = Parameters.sequenceFilePaths[idx].substring(0, lastIdx) + ".unknown.seq";
						Parameters.exportSAMPaths[idx]	  = Parameters.sequenceFilePaths[idx].substring(0, lastIdx) + ".ided.sam";
						Parameters.tmpOutputFilePaths[idx]= Parameters.sequenceFilePaths[idx].substring(0, lastIdx) + "."+Constants.UNIQUE_RUN_ID;
						if(!isExist(Parameters.sequenceFilePaths[idx])) {
							printNoSuchFileOrDirectory(Parameters.sequenceFilePaths[idx]);
							return -1;
						}
					}
				}
			}

			for(int i=0; i<args.length; i+=2) {
				String option = args[i].toLowerCase();

				// --gtf_file (mandatory)
				if(option.equalsIgnoreCase(Parameters.CMD_GENOMIC_ANNOTATION_PATH)) {
					Parameters.genomicAnnotationFilePath = args[i+1];
					if(!isExist(Parameters.genomicAnnotationFilePath)) {
						printNoSuchFileOrDirectory(Parameters.genomicAnnotationFilePath);
						return -1;
					}
				}
				// --psm_file (mandatory)
				else if(option.equalsIgnoreCase(Parameters.CMD_PEPTIDE_ANNOTATION_PATH)) {
					Parameters.peptideFilePath = args[i+1];
					if(!isExist(Parameters.peptideFilePath)) {
						printNoSuchFileOrDirectory(Parameters.peptideFilePath);
						return -1;
					}
				}
				// --sep (optional)
				else if(option.equalsIgnoreCase(Parameters.SEP_TYPE)) {
					Parameters.sepType = args[i+1];
				}
				// --mode (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_TRANSLATION)) {
					String sixFT= String.valueOf(Constants.SIX_FRAME);
					String threeFT= String.valueOf(Constants.THREE_FRAME);
					if(!args[i+1].equalsIgnoreCase(sixFT) && !args[i+1].equalsIgnoreCase(threeFT)) {
						System.out.println(args[i+1] +" is wrong value. Enforce to three-frame translation.");
						args[i+1] = threeFT;
					}

					Parameters.translationMethod = Integer.parseInt(args[i+1]);
				}
				// --fasta_file (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_PROTEIN_SEQUENCE_PATH)) {
					Parameters.proteinFastaPath = args[i+1];
					if(!isExist(Parameters.proteinFastaPath)) {
						printNoSuchFileOrDirectory(Parameters.proteinFastaPath);
						return -1;
					}
				}
				// --output (mandatory)
				else if(option.equalsIgnoreCase(Parameters.CMD_OUTPUT_PATH)) {
					Parameters.outputFilePath = args[i+1] +".pXg";
					Parameters.pinFilePath = Parameters.outputFilePath +".pin";
					if(isExist(Parameters.outputFilePath)) {
						printAlreadyExistFileOrDirectory(Parameters.outputFilePath);
						return -1;
					}
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
				// -scan_col (mandatory)
				else if(option.equalsIgnoreCase(Parameters.CMD_SCAN_COLUMN_INDEX)) {
					Parameters.scanColumnIndex = Integer.parseInt(args[i+1]);
				}
				// -file_col (mandatory)
				else if(option.equalsIgnoreCase(Parameters.CMD_FILE_COLUMN_INDEX)) {
					Parameters.fileColumnIndex = Integer.parseInt(args[i+1]);
				}
				// -charge_col (mandatory)
				else if(option.equalsIgnoreCase(Parameters.CMD_CHARGE_COLUMN_INDEX)) {
					Parameters.chargeColumnIndex = Integer.parseInt(args[i+1]);
				}
				// -add_feat_cols (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_ADD_FEAT_COLUMN_INDICES)) {
					String[] indicies = args[i+1].split("\\,");
					Parameters.additionalFeatureIndices = new int[indicies.length];
					for(int idx=0; idx<indicies.length; idx++) {
						Parameters.additionalFeatureIndices[idx] = Integer.parseInt(indicies[idx]);
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
					if(args[i+1].equalsIgnoreCase("true")) {
						Parameters.EXPORT_SAM = true;
					}
				}
				// -out_sam (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_UNMAPPED)) {
					if(args[i+1].equalsIgnoreCase("false")) {
						Parameters.EXPORT_SAM = false;
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
				// -punmap (optional)
				else if(option.equalsIgnoreCase(Parameters.CMD_PENALTY_SOFTCLIP)) {
					Double penalty = Double.parseDouble(args[i+1]);
					Parameters.PENALTY_SOFTCLIP= penalty;
				}
				else if(option.equalsIgnoreCase(Parameters.CMD_MAX_FLANK_SIZE)) {
					Integer mFlankSize = Integer.parseInt(args[i+1]);
					Parameters.maxFlankNSize= mFlankSize;
				}
				// hidden parameters for revision
				else if(option.equalsIgnoreCase(Parameters.CMD_PHRED_CAL)) {
					Parameters.PHRED_CAL = args[i+1].toLowerCase();
					System.out.println("!!Hidden parameter phred cal: "+Parameters.PHRED_CAL);
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
		Parameters.scanColumnIndex--;
		Parameters.chargeColumnIndex--;
		Parameters.fileColumnIndex--;
		if(Parameters.additionalFeatureIndices != null) {
			for(int i=0; i<Parameters.additionalFeatureIndices.length; i++) {
				Parameters.additionalFeatureIndices[i]--;
			}
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
		String samPaths = Parameters.sequenceFilePaths[0];
		for(int i=1; i<Parameters.sequenceFilePaths.length; i++) {
			samPaths += "," + Parameters.sequenceFilePaths[i];
		}
		System.out.println(" SAM: "+samPaths);
		System.out.println("  SAM_PARTITION_SIZE: "+Parameters.readSize);

		String translation = "NA";
		if(Parameters.translationMethod == Constants.THREE_FRAME) {
			translation = "three-frame translation";
		} else if(Parameters.translationMethod == Constants.SIX_FRAME) {
			translation = "six-frame translation";
		}
		System.out.println("  TRANSLATION_MODE: "+Parameters.translationMethod +" ("+translation+")");
		if(Parameters.proteinFastaPath != null) {
			System.out.println(" PROTEIN_DB: "+Parameters.proteinFastaPath);
		}
		System.out.println(" PSM: "+Parameters.peptideFilePath);
		System.out.println("  PEPT_COL: "+Parameters.peptideColumnIndex);
		System.out.println("  SCORE_COL: "+Parameters.scoreColumnIndex);
		System.out.println("  SCAN_COL: "+Parameters.scanColumnIndex);
		System.out.println("  FILE_COL: "+Parameters.fileColumnIndex);
		System.out.println("  CHARGE_COL: "+Parameters.chargeColumnIndex);
		// to display array
		String addFeatCols = "NA";
		if(Parameters.additionalFeatureIndices != null) {
			addFeatCols = "";
			for (int additionalFeatureIndex : Parameters.additionalFeatureIndices) {
				addFeatCols += "," + additionalFeatureIndex;
			}
			addFeatCols = addFeatCols.substring(1);
		}
		System.out.println("  ADDITIONAL_FEATURE_COLS: "+addFeatCols);
		System.out.println("  RANK TO CONSIDER: "+Parameters.psmRank);
		System.out.println("  PEPTIDE_LENGTHS: "+Parameters.minPeptLen+"-"+Parameters.maxPeptLen);
		System.out.println("  MAX_FLANK_SIZE: "+Parameters.maxFlankNSize);
		System.out.println(" OUT_RESULT: "+Parameters.outputFilePath);
		System.out.println("  OUT_PIN.: "+Parameters.pinFilePath);
		System.out.println("  OUT_UNKNOWN: "+Parameters.EXPORT_UNMAPPED_SEQ);
		System.out.println("  OUT_SAM: "+Parameters.EXPORT_SAM);
		System.out.println("  OUT_CANONICAL: "+Parameters.EXPORT_CANONICAL);
		System.out.println("  OUT_NONCANONICAL: "+Parameters.EXPORT_NONCANONICAL);
		System.out.println(" penalty_mutation: "+Parameters.PENALTY_MUTATION);
		System.out.println(" penalty_AS: "+Parameters.PENALTY_AS);
		System.out.println(" penalty_5UTR: "+Parameters.PENALTY_5UTR);
		System.out.println(" penalty_3UTR: "+Parameters.PENALTY_3UTR);
		System.out.println(" penalty_ncRNA: "+Parameters.PENALTY_ncRNA);
		System.out.println(" penalty_FS: "+Parameters.PENALTY_FS);
		System.out.println(" penalty_IR: "+Parameters.PENALTY_IR);
		System.out.println(" penalty_IGR: "+Parameters.PENALTY_IGR);
		System.out.println(" penalty_asRNA: "+Parameters.PENALTY_asRNA);
		System.out.println(" penalty_softclip: "+Parameters.PENALTY_SOFTCLIP);
		System.out.println(" penalty_unknown: "+Parameters.PENALTY_UNMAP);
		System.out.println(" THREADS: "+Parameters.nThreads);

		// append to logger
		Logger.append("Running info");
		Logger.newLine();
		Logger.append(" GTF: "+Parameters.genomicAnnotationFilePath);
		Logger.newLine();
		Logger.append("  GTF_PARTITION_SIZE: "+Parameters.partitionSize);
		Logger.newLine();
		Logger.append(" SAM: "+samPaths);
		Logger.newLine();
		Logger.append("  SAM_PARTITION_SIZE: "+Parameters.readSize);
		Logger.newLine();
		Logger.append("  TRANSLATION_MODE: "+Parameters.translationMethod +" ("+translation+")");
		Logger.newLine();
		if(Parameters.proteinFastaPath != null) {
			Logger.append(" PROTEIN_DB: "+Parameters.proteinFastaPath);
			Logger.newLine();
		}
		Logger.append(" PSM: "+Parameters.peptideFilePath);
		Logger.newLine();
		Logger.append("  PEPT_COL: "+Parameters.peptideColumnIndex);
		Logger.newLine();
		Logger.append("  SCORE_COL: "+Parameters.scoreColumnIndex);
		Logger.newLine();
		Logger.append("  ADDITIONAL_FEATURE_COLS: "+addFeatCols);
		Logger.newLine();
		Logger.append("  RANK TO CONSIDER: "+Parameters.psmRank);
		Logger.newLine();
		Logger.append("  PEPTIDE_LENGTHS: "+Parameters.minPeptLen+"-"+Parameters.maxPeptLen);
		Logger.newLine();
		Logger.append("  MAX_FLANK_SIZE: "+Parameters.maxFlankNSize);
		Logger.newLine();
		Logger.append(" OUT_RESULT: "+Parameters.outputFilePath);
		Logger.newLine();
		Logger.append("  OUT_PIN: "+Parameters.pinFilePath);
		Logger.newLine();
		Logger.append("  OUT_UNMAPPED: "+Parameters.EXPORT_UNMAPPED_SEQ);
		Logger.newLine();
		Logger.append("  OUT_SAM: "+Parameters.EXPORT_SAM);
		Logger.newLine();
		Logger.append("  OUT_CANONICAL: "+Parameters.EXPORT_CANONICAL);
		Logger.newLine();
		Logger.append("  OUT_NONCANONICAL: "+Parameters.EXPORT_NONCANONICAL);
		Logger.newLine();
		Logger.append(" penalty_mutation: "+Parameters.PENALTY_MUTATION);
		Logger.newLine();
		Logger.append(" penalty_AS: "+Parameters.PENALTY_AS);
		Logger.newLine();
		Logger.append(" penalty_5UTR: "+Parameters.PENALTY_5UTR);
		Logger.newLine();
		Logger.append(" penalty_3UTR: "+Parameters.PENALTY_3UTR);
		Logger.newLine();
		Logger.append(" penalty_ncRNA: "+Parameters.PENALTY_ncRNA);
		Logger.newLine();
		Logger.append(" penalty_FS: "+Parameters.PENALTY_FS);
		Logger.newLine();
		Logger.append(" penalty_IR: "+Parameters.PENALTY_IR);
		Logger.newLine();
		Logger.append(" penalty_IGR: "+Parameters.PENALTY_IGR);
		Logger.newLine();
		Logger.append(" penalty_asRNA: "+Parameters.PENALTY_asRNA);
		Logger.newLine();
		Logger.append(" penalty_unknown: "+Parameters.PENALTY_UNMAP);
		Logger.newLine();
		Logger.append(" THREADS: "+Parameters.nThreads);
		Logger.newLine();
	}

	private static boolean isMandatoryOkay () {
		boolean pass = true;

		// --gtf_file
		if(Parameters.genomicAnnotationFilePath == null) {
			System.out.println("mandatory option --gtf_file is missing...");
			pass = false;
		}
		// --sam_file
		if(Parameters.sequenceFilePaths == null) {
			System.out.println("mandatory option --sam_file is missing...");
			pass = false;
		}
		// --psm_file
		if(Parameters.peptideFilePath == null) {
			System.out.println("mandatory option --psm_file is missing...");
			pass = false;
		}
		// --out_file
		if(Parameters.outputFilePath == null) {
			System.out.println("mandatory option --out_file is missing...");
			pass = false;
		}
		// -pept_col
		if(Parameters.peptideColumnIndex == -1) {
			System.out.println("mandatory option --pept_col is missing...");
			pass = false;
		}
		// -score_col
		if(Parameters.scoreColumnIndex == -1) {
			System.out.println("mandatory option --score_col is missing...");
			pass = false;
		}
		// -scan_col
		if(Parameters.scanColumnIndex == -1) {
			System.out.println("mandatory option --scan_col is missing...");
			pass = false;
		}
		// -file_col
		if(Parameters.fileColumnIndex == -1) {
			System.out.println("mandatory option --file_col is missing...");
			pass = false;
		}
		// -charge_col
		if(Parameters.chargeColumnIndex == -1) {
			System.out.println("mandatory option --charge_col is missing...");
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
