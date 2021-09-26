package progistar.pXg.constants;

public class Parameters {
	// Input file paths
	public static final String GENOMIC_ANNOTATION_PATH				=	"genomic_annotation_path".toLowerCase();
	public static String genomicAnnotationFilePath	=	"C:\\Bioinformatics\\0.Databases\\1.HumanRNAs\\gencode.v19.annotation.sorted.gtf";
	
	public static final String GENOMIC_SEQUENCE_PATH				=	"genomic_sequence_path".toLowerCase();
	public static String sequenceFilePath			=	"C:\\Users\\progi\\Desktop\\Projects\\pXg\\chr1toy.sam";
	
	public static final String PEPTIDE_ANNOTATION_PATH				=	"peptide_annotation_path".toLowerCase();
	public static String peptideFilePath			=	"";
	
	// Output file path
	public static final String OUTPUT_PATH				=	"output_path".toLowerCase();
	public static String outputFilePath				=	"";
	
	public static final String DECOY_METHOD			=	"reverse".toLowerCase();
	public static byte	decoys						=	Constants.DECOY_REVERSE;
	
	// GTF partition size
	public static final String GENOMIC_ANNOTATION_PARTITION_SIZE		=	"genomic_annotatino_partition_size".toLowerCase();
	public static int partitionSize					=	5000000; // 5 * 10^6 * 10 * 8 = 0.4 G
	
	public static final String GENOMIC_SEQUENCE_PARTITION_SIZE			=	"genomic_sequence_partition_size".toLowerCase();
	public static int readSize						=	1000000; // 1 * 10^7
	
	// READ sequencing 
	public static byte READ_STRAND					=	Constants.FORWARD_STRAND_READS;
	
	public static int minPeptLen					=	8;
	public static int maxPeptLen					=	15;
	public static boolean leucineIsIsoleucine		=	true;
	
	// Peptide file
	// for user-friendly purpose, peptideColumnIndex is taken one-based and converted to zero-based.
	public static final String PEPTIDE_COLUMN_INDEX			=	"peptide_column_index".toLowerCase();
	public static int peptideColumnIndex			=	0; // user-specific peptide index
	
	public static final String SCAN_COLUMN_INDICES	=	"scan_column_indices".toLowerCase();
	public static int[] scanColumnIndices			=	null; // user-specific scan index
	
	public static final String SCORE_COLUMN_INDEX		=	"score_column_index".toLowerCase();
	public static int scoreColumnIndex				=	0;
	
	public static final String DELTA_SCORE_THRESHOLD		=	"score_threshold".toLowerCase();
	public static double deltaScoreThreshold				=	10;
	
	public static final String pParserRegExr		=	"aareg".toLowerCase();
	public static String peptideParserRegExr		=	"[A-Z]"; // read sequence matched to the RegExr.
	// Note that
	// comment is different from field.
	// In pXg definition, field indicates column names and comment shows some meta-data.
	// Comment line must be present on the top of records.
	// Field must follow "comment" when the comment presents in the file.
	public static final String cMarker				=	"cm".toLowerCase();
	public static String commentMarker				=	"#|@|%"; // if line starts with the pattern, the line will be skipped during parsing the file.
	
	// System Parameters
	public static final String numOfThreads			=	"num_threads".toLowerCase();
	public static int nThreads						=	1;
	
}
