package progistar.pXg.constants;

public class Parameters {
	// Input file paths
	public static final String CMD_GENOMIC_ANNOTATION_PATH				=	"-gtf".toLowerCase();
	public static String genomicAnnotationFilePath	=	null;
	
	public static final String CMD_GENOMIC_SEQUENCE_PATH				=	"-sam".toLowerCase();
	public static String sequenceFilePath			=	null;
	
	public static final String CMD_PEPTIDE_ANNOTATION_PATH				=	"-psm".toLowerCase();
	public static String peptideFilePath			=	null;
	
	public static final String CMD_PROTEIN_SEQUENCE_PATH		=	"-fasta".toLowerCase();
	public static String proteinFastaPath			=	null;
	
	public static final String SEP_TYPE							=	"-sep".toLowerCase();
	public static String sepType					=	"csv".toLowerCase();	
	// Output file path
	public static final String CMD_OUTPUT_PATH				=	"-out".toLowerCase();
	public static String outputFilePath					=	null;
	public static String ngsStatFilePath				=	null;
	public static String psmStatFilePath				=	null;
	public static String unmappedFilePath				=	null;
	public static String exportGTFPath					=	null;
	public static String exportSAMPath					=	null;
	public static String exportVCFPath					=	null;
	

	// Significant NGS-read mapping
	public static final String CMD_P_VALUE			=	"-pval".toLowerCase();
	public static double pvalue						=	0.01;
	
	// FDR estimation at PSM level
	public static final String CMD_FDR				=	"-fdr".toLowerCase();
	public static double fdr						=	0.1;
	
	// Peptide length
	public static final String CMD_LENGTH			=	"-length".toLowerCase();
	public static int minPeptLen					=	8;
	public static int maxPeptLen					=	15;
	
	public static final String CMD_IL				=	"-ileq".toLowerCase();
	public static boolean leucineIsIsoleucine		=	true;
	
	// Those options cannot be accessed by users.
	// It is only used for test which method is more adaptable.
	public static byte	mocks						=	Constants.MOCK_PSD_REVERSE;
	public static byte	mockPolicy				=	Constants.MOCK_MAX_ONE;
	
	// Output format
	public static final String CMD_SAM_FORMAT		=	"-out_sam".toLowerCase();
	public static final String CMD_GTF_FORMAT		=	"-out_gtf".toLowerCase();
	public static final String CMD_NONCANONICAL		=	"-out_noncanonical".toLowerCase();
	public static final String CMD_CANONICAL		=	"-out_canonical".toLowerCase();
	
	// GTF partition size
	public static final String CMD_GENOMIC_ANNOTATION_PARTITION_SIZE		=	"-gtf_partition_size".toLowerCase();
	public static int partitionSize					=	10000000; // 10^7 * 10 * 8 = 0.8 G
	
	public static final String CMD_GENOMIC_SEQUENCE_PARTITION_SIZE			=	"-sam_partition_size".toLowerCase();
	public static int readSize						=	1000000; // 1 * 10^6
	
	// maxJunctionSize MUST not exceed partitionSize in GTF
	public static int maxJunctionSize				=	1000000; // 1 * 10^6
	
	// Peptide file
	// for user-friendly purpose, peptideColumnIndex is taken one-based and converted to zero-based.
	public static final String CMD_PEPTIDE_COLUMN_INDEX	=	"-pept_col".toLowerCase();
	public static int peptideColumnIndex			=	-1; // user-specific peptide index
	
	public static final String CMD_SCAN_COLUMN_INDICES	=	"-scan_cols".toLowerCase();
	public static int[] scanColumnIndices			=	null; // user-specific scan index
	
	public static final String CMD_SCORE_COLUMN_INDEX	=	"-score_col".toLowerCase();
	public static int scoreColumnIndex				=	-1;
	
	public static final String CMD_CANDIDATE_RANK	=	"-rank".toLowerCase();
	public static int psmRank						=	100;
	
	public static final String pParserRegExr		=	"aareg".toLowerCase();
	public static String peptideParserRegExr		=	"[A-Z]"; // read sequence matched to the RegExr.
	// Note that
	// comment is different from field.
	// In pXg definition, field indicates column names and comment shows some meta-data.
	// Comment line must be present on the top of records.
	// Field must follow "comment" when the comment presents in the file.
	public static final String cMarker				=	"cm".toLowerCase();
	public static String commentMarker				=	"#|@|%"; // if line starts with the pattern, the line will be skipped during parsing the file.
	
	public static int maxProteinOut					=	10;
	
	// Penalty
	public static final String CMD_PENALTY_MUTATION	=	"-pMut".toLowerCase();
	public static final String CMD_PENALTY_AS		=	"-pAS".toLowerCase();
	public static final String CMD_PENALTY_5UTR		=	"-p5UTR".toLowerCase();
	public static final String CMD_PENALTY_3UTR		=	"-p3UTR".toLowerCase();
	public static final String CMD_PENALTY_FS		=	"-pFS".toLowerCase();
	public static final String CMD_PENALTY_ncRNA	=	"-pncRNA".toLowerCase();
	public static final String CMD_PENALTY_IR		=	"-pIR".toLowerCase();
	public static final String CMD_PENALTY_IGR		=	"-pIGR".toLowerCase();
	public static final String CMD_PENALTY_asRNA	=	"-pasRNA".toLowerCase();
	public static final String CMD_PENALTY_UNMAP	=	"-punmap".toLowerCase();
	
	public static double PENALTY_MUTATION			=	1;
	public static double PENALTY_AS					=	10;
	public static double PENALTY_5UTR				=	20;
	public static double PENALTY_3UTR				=	20;
	public static double PENALTY_FS					=	20;
	public static double PENALTY_ncRNA				=	20;
	public static double PENALTY_IR					=	30;
	public static double PENALTY_IGR				=	30;
	public static double PENALTY_asRNA				=	30;
	public static double PENALTY_UNMAP				=	100;
	
	// Export option
	public static boolean EXPORT_TSV				=	true;
	public static boolean EXPORT_GTF				=	true;
	public static boolean EXPORT_SAM				=	true;
	public static boolean EXPORT_NONCANONICAL		=	true;
	public static boolean EXPORT_CANONICAL			=	true;
	
	// System Parameters
	public static final String CMD_THREADS			=	"-threads".toLowerCase();
	public static int nThreads						=	1;
	
	public static boolean isDecoyOut 				=	true;
	
	public static final String CMD_TRANSLATION		=	"-translation";

	// Three or Six frame translation
	public static int translationMethod				=	Constants.SIX_FRAME;
}
